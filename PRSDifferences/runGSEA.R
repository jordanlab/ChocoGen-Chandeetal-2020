library(biomaRt)
library(tidyverse)

ToI = read.delim('./Filt-ImputedMatchedSNP-PRS.tsv')
ToI.clm = ToI %>% filter(Diff > 0)
ToI.cho = ToI %>% filter(Diff < 0)

df = read.delim("D:/Dropbox (ABiL)/JordanLab/ChocoGen/Analysis/ChocoGen3/PRSDifferences/20190416-allele_info-withflips.tsv")


mart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")

snps <- getBM(attributes=c('refsnp_id', 'ensembl_gene_stable_id'), 
              filters = 'snp_filter', 
              values = unique(sort(df$rsid)), 
              mart = mart
              )

snps.filt = snps[snps$ensembl_gene_stable_id != "",]

mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
geneids = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", 'entrezgene_id'), 
                filters = "ensembl_gene_id", 
                values = snps.filt$ensembl_gene_stable_id, 
                mart = mart
                )

gene.snps = merge(snps.filt, 
                  geneids, 
                  by.x = "ensembl_gene_stable_id", 
                  by.y="ensembl_gene_id"
                  )


anno = merge(df,gene.snps, by.y="refsnp_id", by.x="rsid") %>% 
    mutate(CLM_expr = ifelse( 
                      effect_allele == major_allele,
                      CLM_major/CHO_major,
                      CLM_minor/CHO_minor 
                    ),
           CHO_expr = ifelse( 
               effect_allele == major_allele,
               CHO_major/CLM_major,
               CHO_minor/CLM_minor 
               )
           )

anno <- anno[!duplicated(anno) & anno$hgnc_symbol != "" & anno$entrezgene_id != "", ] %>% na.omit()


# Psuedo-dummy values for 0 value or inf valued 
anno[is.infinite(anno$CHO_expr), 'CHO_expr'] <- runif(length(anno[is.infinite(anno$CHO_expr), 'CHO_expr']), 9.100069734290, 11.703201230)
anno[is.nan(anno$CHO_expr), 'CHO_expr'] <- runif(length(anno[is.nan(anno$CHO_expr), 'CHO_expr']), 0.01234, 0.059238490)


CLM.merge = merge(ToI.clm, 
                  anno,
                  by.x="Trait",
                  by.y = "trait"
                  ) %>% 
            filter(CLM_expr > 1, is.finite(CLM_expr)) %>% 
            dplyr::group_by(entrezgene_id) %>%
            filter(CLM_expr == max(CLM_expr)) %>%
            ungroup() %>%
            distinct(entrezgene_id, .keep_all = T) %>%
            arrange(-CLM_expr)

CHO.merge = merge(ToI.cho, 
                  anno,
                  by.x="Trait",
                  by.y = "trait"
            ) %>% 
              filter(CHO_expr > 1, is.finite(CHO_expr)) %>% 
              dplyr::group_by(entrezgene_id) %>%
              filter(CHO_expr == max(CHO_expr)) %>%
              ungroup() %>%
              distinct(entrezgene_id, .keep_all = T) %>%
              arrange(-CHO_expr)




library(clusterProfiler)
library(msigdbr)
m_df = msigdbr(species = "Homo sapiens", ) %>% 
  filter(gs_subcat 
         %in% c("CP:KEGG", "CP:BIOCARTA", "CP:PID") |
           gs_cat == "H"
         ) %>%
  dplyr::select(gs_name, entrez_gene)


geneList.cho <- CHO.merge$CHO_expr
names(geneList.cho) <- CHO.merge$entrezgene_id
gene.cho <- CHO.merge$entrezgene_id
geneList.cho.n <- CHO.merge$CHO_expr
names(geneList.cho.n) <- CHO.merge$hgnc_symbol
gene.cho.n <- CHO.merge$hgnc_symbol

geneList.clm <- CLM.merge$CLM_expr
names(geneList.clm) <- CLM.merge$entrezgene_id
gene.clm <- CLM.merge$entrezgene_id
geneList.clm.n <- CLM.merge$CLM_expr
names(geneList.clm.n) <- CLM.merge$hgnc_symbol
gene.clm.n <- CLM.merge$hgnc_symbol


em.CLM<- enricher(gene.clm, TERM2GENE = m_df, pAdjustMethod = "BH", pvalueCutoff = 10)
em.CHO<- enricher(gene.cho, TERM2GENE = m_df, pAdjustMethod = "BH", pvalueCutoff = 10)

traits.cho = em.CHO@result %>% arrange(ID) %>% dplyr::select(ID) %>% unlist()
traits.clm  = em.CLM@result %>% arrange(ID) %>% dplyr::select(ID) %>% unlist()

traits.common = intersect(traits.cho, traits.clm)


pval.cho = em.CHO@result %>% arrange(ID) %>% 
  filter(ID %in% traits.common) %>%
  dplyr::select(pvalue) %>%
  unlist()

pval.clm = em.CLM@result %>% arrange(ID) %>% 
  filter(ID %in% traits.common) %>%
  dplyr::select(pvalue) %>%
  unlist()

delta.p = log(pval.clm, 2) - log(pval.cho, 2)


pdata = data.frame("clm.logp" = -log(pval.clm, 2),
                   "clm.p" = pval.clm,
                   "cho.logp" = -log(pval.cho, 2),
                   "ho.p" = pval.cho,
                   "delta.p" = -log(pval.clm, 2) - -log(pval.cho, 2),
                   "ID" = traits.common, row.names = traits.common)


d <- ggplot(pdata, aes(x=delta.p, label=ID)) + 
  geom_histogram(binwidth=.1) + 
  theme_minimal() + xlab("Δ log2(P-value)") + ylim(0, 15) + 
  geom_vline(xintercept = c(mean(delta.p), mean(delta.p)-2*sd(delta.p), mean(delta.p)+2*sd(delta.p),  mean(delta.p)-3*sd(delta.p),  mean(delta.p)+3*sd(delta.p))) + 
  geom_text(aes(x=mean(delta.p)-2*sd(delta.p)-.5, label = paste("Δ log2(P-value) < ", round(mean(delta.p)-2*sd(delta.p), 5)), y=4), angle=90) + 
  geom_text(aes(x=mean(delta.p)+2*sd(delta.p)+.5, label = paste("Δ log2(P-value) > ", round(mean(delta.p)+2*sd(delta.p), 5)), y=4), angle=90) + 
  geom_text(aes(x=mean(delta.p)+3*sd(delta.p)+.5, label = paste("Δ log2(P-value) > ", round(mean(delta.p)+3*sd(delta.p), 5)), y=4), angle=90)


d
ggsave("P-value-distribution.pdf")

ggplotly(d)



selected_paths = c("PID_BETA_CATENIN_NUC_PATHWAY",
"KEGG_LONG_TERM_POTENTIATION",
"PID_EPHA_FWDPATHWAY",
"PID_NFAT_3PATHWAY",
"KEGG_STARCH_AND_SUCROSE_METABOLISM",
"HALLMARK_KRAS_SIGNALING_UP",
"KEGG_COLORECTAL_CANCER",
"KEGG_TYPE_II_DIABETES_MELLITUS",
'KEGG_GLYCEROLIPID_METABOLISM',
"PID_AJDISS_2PATHWAY",
"KEGG_NON_SMALL_CELL_LUNG_CANCER",
"KEGG_PROSTATE_CANCER",
"KEGG_PORPHYRIN_AND_CHLOROPHYLL_METABOLISM",
"KEGG_STEROID_HORMONE_BIOSYNTHESIS",
"KEGG_ACUTE_MYELOID_LEUKEMIA",
"KEGG_PENTOSE_AND_GLUCURONATE_INTERCONVERSIONS",
"KEGG_ENDOMETRIAL_CANCER",
"PID_AVB3_INTEGRIN_PATHWAY",
"PID_INTEGRIN3_PATHWAY")


sel.cho = em.CHO@result %>% arrange(ID) %>% 
  filter(ID %in% selected_paths) %>% 
  dplyr::select(Description, GeneRatio, BgRatio, p.adjust) %>%
  mutate(pop="CHO") 


sel.clm = em.CLM@result %>% arrange(ID) %>% 
  filter(ID %in% selected_paths) %>% 
  dplyr::select(Description, GeneRatio, BgRatio, p.adjust) %>%
  mutate(pop="CLM") 

selected <- rbind(sel.cho, sel.clm) %>% group_by(pop)

# selected$GeneRatio <- sapply(selected$GeneRatio, function(x) eval(parse(text=as.character(x))))

selected$GeneRatio <- gsub("/.*", "", selected$GeneRatio)
selected$BgRatio<- gsub("/.*", "", selected$BgRatio)
selected$Ratio = as.numeric(selected$GeneRatio)/as.numeric(selected$BgRatio)


selected = selected %>% arrange(Ratio, pop) %>% mutate(Description = factor(Description, unique(Description)))

ggplot(selected, aes(x=pop, y=Description, size=Ratio, color=-log2(p.adjust))) +
  geom_point() +
  scale_colour_gradient2(name = "-log2(Adjusted P-value)", mid="white", low = "grey", midpoint = 4, high="blue") +
  # scale_size(name="Geneset overlap size") +
  ylab(NULL) + xlab(NULL) + ggtitle("") + DOSE::theme_dose(8) + scale_size(range = c(1,10))


ggsave("PathwayEnrichment-Bubbleplot.pdf", dpi=600, width = 8, height = 4, useDingbats = F)

