library(data.table)
library(openxlsx)
library(ggplot2)
library(ggrepel)
exp_tab = openxlsx::read.xlsx(xlsxFile = "data/Transcript_based_TPM_Farm_data_GeneNames.xlsx",
                              sheet = "All_ESF_with_new_RabbitESF")




motif_counts = data.frame()
motifmats = list()

for (Species in c("Cat", "Cow", "Pig", "GuineaPig", "Horse", "Rat", "Rabbit", "Sheep", "Human", "Dog")) {
  load(paste0("data/countmatrices/", Species, ".Rdata"), verbose = T)
  mmat <- 1e3 * sweep(x = mmat, MARGIN = 1, STATS = rowSums(mmat), FUN = "/") ## NOTE:: this line will nornalize the motif counts, which we used for ELI genes
  motifmats[[Species]] <- mmat
}


for (Species in c("Cat", "Cow", "Pig", "Dog", "GuineaPig", "Horse", "Rat", "Rabbit", "Sheep")) {
  mcount_this = data.frame(t(motifmats[[Species]]))
  mcount_this$Species = Species
  mcount_this$ensembl_gene_id = rownames(mcount_this)
  #mcount_this = data.frame(mcounts = motifmats[[Species]]["FOXF2",], Species = Species, ensembl_gene_id = colnames(motifmats[[Species]]))
  #mcount_this = merge(x=mcount_this, y=eval(parse(text = paste0(tolower(Species),"_name_table"))), by="ensembl_gene_id")
  motif_counts <- rbind(motif_counts, mcount_this)
  #[, c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene", "Species", "mcounts")])
}

load("data/countmatrices/Human.Rdata", verbose = T)
mmat <- 1e3 * sweep(x = mmat, MARGIN = 1, STATS = rowSums(mmat), FUN = "/") ## NOTE:: this line will nornalize the motif counts, which we used for ELI genes
mcount_this = data.frame(t(mmat))
mcount_this$Species = "Human"
#mcount_this$hsapiens_homolog_ensembl_gene = rownames(mcount_this)
mcount_this$ensembl_gene_id = rownames(mcount_this)

motif_counts <- rbind(motif_counts, mcount_this)

rm(motifmats)
rm(mcount_this)

#all_exp <- all_exp[all_exp$hsapiens_homolog_ensembl_gene != "",]
#all_exp <- all_exp[all_exp$ensembl_gene_id != "",]
motif_counts <- motif_counts[motif_counts$ensembl_gene_id != "",]



exp_tab = as.data.table(exp_tab)
names(exp_tab)[names(exp_tab) == "GeneID"] <- "ensembl_gene_id"


all_exp = data.table()

this_exp = melt(exp_tab[, c("ensembl_gene_id", grep(pattern = "Human", x = names(exp_tab), value = T)), with=F], id.vars = "ensembl_gene_id")
this_exp$Species = "Human"
this_exp$hsapiens_gene_id = this_exp$ensembl_gene_id
this_exp$variable <- NULL
all_exp = this_exp


#for (Species in c("Cow", "Horse", "Rat", "Rabbit", "Opossum", "GuineaPig", "Sheep", "Dog")) {
for (Species in c("Cow", "Horse", "Rat", "Rabbit", "GuineaPig", "Sheep", "Dog", "Cat")) {
  this_exp = melt(exp_tab[, c("ensembl_gene_id", grep(pattern = Species, x = names(exp_tab), value = T)), with=F], 
                  id.vars = c("ensembl_gene_id", paste0(Species,"GeneID")))
  this_exp$Species = Species
  names(this_exp)[names(this_exp) == "ensembl_gene_id"] = "hsapiens_gene_id"
  names(this_exp)[names(this_exp) == paste0(Species, "GeneID")] = "ensembl_gene_id"
  
  this_exp$variable <- NULL
  all_exp = rbind(all_exp, this_exp)
}



### Compare Nature Evolution Data


hum_bov <- openxlsx::read.xlsx(xlsxFile = "data/TPM Values_Gunter_human decidual bovine_IMPORTANT.xlsx", sheet = "TPM Values hvsb")
hum_bov <- as.data.table(hum_bov)
old_cow = all_exp[Species == "Cow"]
new_cow = hum_bov[, c("Ensembl.Gene.ID", "BtESF1_TPM", "BtESF2_TPM", "BtESF3_TPM")]



names(new_cow)[names(new_cow)=="Ensembl.Gene.ID"] <- "hsapiens_gene_id"

delme <- as.data.table(old_cow)
delme <- unique(delme[, .(ensembl_gene_id, hsapiens_gene_id)])

new_cow <- merge(x = new_cow, y = delme, by="hsapiens_gene_id")
new_cow$BtESF1_TPM = 1e6 * new_cow$BtESF1_TPM/ sum(new_cow$BtESF1_TPM)
new_cow$BtESF2_TPM = 1e6 * new_cow$BtESF2_TPM/ sum(new_cow$BtESF2_TPM)
new_cow$BtESF2_TPM = 1e6 * new_cow$BtESF3_TPM/ sum(new_cow$BtESF3_TPM)

new_cow <- melt(data = as.data.table(new_cow), id.vars=c("hsapiens_gene_id", "ensembl_gene_id"))
new_cow <- new_cow[, -"variable"]
new_cow$Species = "Cow"
new_cow <- new_cow[, c("ensembl_gene_id", "value", "Species", "hsapiens_gene_id")]


old_human = all_exp[Species == "Human"]
new_human = hum_bov[, c("Ensembl.Gene.ID", "Endometrial1_TPM", "Endometrial2_TPM", "Endometrial3_TPM")]
names(new_human)[names(new_human)=="Ensembl.Gene.ID"] <- "hsapiens_gene_id"
new_human$ensembl_gene_id = new_human$hsapiens_gene_id
new_human <- new_human[hsapiens_gene_id %in% old_human$hsapiens_gene_id]
new_human$Endometrial1_TPM <- 1e6 * new_human$Endometrial1_TPM / sum(new_human$Endometrial1_TPM)
new_human$Endometrial2_TPM <- 1e6 * new_human$Endometrial2_TPM / sum(new_human$Endometrial2_TPM)
new_human$Endometrial3_TPM <- 1e6 * new_human$Endometrial3_TPM / sum(new_human$Endometrial3_TPM)

new_human <- melt(data = as.data.table(new_human), id.vars=c("hsapiens_gene_id", "ensembl_gene_id"))
new_human <- new_human[, -"variable"]
new_human$Species = "Human"
new_human <- new_human[, c("ensembl_gene_id", "value", "Species", "hsapiens_gene_id")]


new_all_exp = rbind(all_exp[! (Species %in% c("Human","Cow"))], new_human, new_cow)




new_exp_counts = merge(new_all_exp, motif_counts, by=c("ensembl_gene_id", "Species"))
new_exp_counts[, corrVal := value - mean(value), by=hsapiens_gene_id]

new_exp_counts[, sqrtcorrVal := sqrt(value) - mean(sqrt(value)), by=hsapiens_gene_id]

#### get names
require(biomaRt)
ensembl = biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl", host = "www.ensembl.org")

namedf2 <- biomaRt::getBM(mart = ensembl, attributes = c("ensembl_gene_id", "hgnc_symbol"))
namedf2 <- as.data.table(namedf2)
setkey(namedf2, "ensembl_gene_id")
namedf2[, nsymbols := c(1:length(hgnc_symbol)), by=ensembl_gene_id]
namedf2[, nids := length(ensembl_gene_id), by=hgnc_symbol]
namedf2 <- namedf2[nsymbols==1]
namedf2 <- namedf2[, .(ensembl_gene_id, hgnc_symbol)]


### get eli genes

Species2Score <- c("Rabbit"=2, "GuineaPig"=2, "Rat"=2, "Human"=3, "Cow"=0, "Sheep"=0, "Horse"=0, "Cat"=1, "Dog"=1)



new_score_exp = new_all_exp
new_score_exp$score = Species2Score[new_score_exp$Species]
getELIScores <- function(scs,vals) {
  ct <- cor.test(x = scs, y = vals)
  ct_sqrt <- cor.test(x = scs, y = sqrt(vals))
  ct_spear <- cor.test(x = scs, y = vals, method = "spearman")
  lmfit = lm(vals ~ scs)
  return(list(ct$p.value, ct$estimate, 
              ct_sqrt$p.value, ct_sqrt$estimate, 
              ct_spear$p.value, ct_spear$estimate,
              lmfit$coefficients["scs"]))
}
new_score_exp[, c("ELIpval", "ELIcor", "ELIpval_sqrt", 
                  "ELIcor_sqrt", "ELIpval_spear", "ELIcor_spear", "ELIslope") := 
                getELIScores(score, value), by=hsapiens_gene_id]
new_score_exp[, heme_epi_ratio := mean(value[score!=0], na.rm = T)/ mean(value[score==0], na.rm = T), by=hsapiens_gene_id]

eli_scores = unique(new_score_exp[, c("hsapiens_gene_id", 
                                      "ELIcor", "ELIpval", 
                                      "ELIcor_sqrt", "ELIpval_sqrt", 
                                      "ELIcor_spear", "ELIpval_spear", 
                                      "heme_epi_ratio", "ELIslope"),with=F])
eli_scores <- eli_scores[order(ELIpval, decreasing = F)]





eli_scores <- merge.data.table(x = eli_scores, y = namedf2, 
                               by.x = "hsapiens_gene_id", by.y = "ensembl_gene_id", all.x = T)

eli_scores$ELIpadj <- p.adjust(p = eli_scores$ELIpval, method = "fdr")
eli_scores <- eli_scores[order(ELIpadj, decreasing = F)]

openxlsx::write.xlsx(x = eli_scores, file = "ELI_scores.xlsx", asTable = T)


new_score_exp <- merge(new_score_exp, namedf2, by.x = "hsapiens_gene_id", by.y = "ensembl_gene_id")



new_all_exp[, sampnum := c(1:length(value)), by="ensembl_gene_id"]

new_all_exp <- merge(new_all_exp, namedf2, by.x = "hsapiens_gene_id", by.y = "ensembl_gene_id")

##define the top 100 genes as ELIup and ELIdn

my_eli_genes = eli_scores[ELIpadj < .01]$hsapiens_gene_id
pos_eli_genes = eli_scores[(ELIpadj < .01) & (ELIcor > 0)]$hsapiens_gene_id
neg_eli_genes = eli_scores[(ELIpadj < .01) & (ELIcor < 0)]$hsapiens_gene_id


# do ELI TFBS tests
getTestDF <- function(test_list) {
  testdf = as.data.frame(do.call(rbind, lapply(X = test_list, FUN = function(x) {
    x = summary(x)
    c(p = x$coefficients[2,4], beta = x$coefficients[2,1], r2 = x$r.squared)
  })))
  testdf$TFBS = rownames(testdf)
  testdf = as.data.table(testdf)
  return(testdf)
}


eli_exp_counts <- new_exp_counts[hsapiens_gene_id %in% my_eli_genes]

eli_tests_simple_sqrt = list()

for (motif in names(motif_counts)[1:572]) {
  #for (motif in names(exp_counts)[5:576]) {
  cat(paste0("doing ", motif, "\n"))
  ff = as.formula(paste0("sqrtcorrVal ~ ",motif))
  eli_tests_simple_sqrt[[motif]] <- lm(formula = ff, data = eli_exp_counts)
}

eli_test_df_sqrt = getTestDF(eli_tests_simple_sqrt)
eli_test_df_sqrt$padj <- p.adjust(p = eli_test_df_sqrt$p, method = "fdr")


pos_eli_exp_counts <- new_exp_counts[hsapiens_gene_id %in% pos_eli_genes]
pos_eli_tests_simple_sqrt = list()
for (motif in names(motif_counts)[1:572]) {
  #for (motif in names(exp_counts)[5:576]) {
  cat(paste0("doing ", motif, "\n"))
  ff = as.formula(paste0("sqrtcorrVal ~ ",motif))
  pos_eli_tests_simple_sqrt[[motif]] <- lm(formula = ff, data = pos_eli_exp_counts)
}
pos_eli_test_df_sqrt = getTestDF(pos_eli_tests_simple_sqrt)
pos_eli_test_df_sqrt$padj <- p.adjust(p = pos_eli_test_df_sqrt$p, method = "fdr")


neg_eli_exp_counts <- new_exp_counts[hsapiens_gene_id %in% neg_eli_genes]
neg_eli_tests_simple_sqrt = list()
for (motif in names(motif_counts)[1:572]) {
  #for (motif in names(exp_counts)[5:576]) {
  cat(paste0("doing ", motif, "\n"))
  ff = as.formula(paste0("sqrtcorrVal ~ ",motif))
  neg_eli_tests_simple_sqrt[[motif]] <- lm(formula = ff, data = neg_eli_exp_counts)
}
neg_eli_test_df_sqrt = getTestDF(neg_eli_tests_simple_sqrt)
neg_eli_test_df_sqrt$padj <- p.adjust(p = neg_eli_test_df_sqrt$p, method = "fdr")



save(eli_test_df_sqrt, pos_eli_test_df_sqrt, neg_eli_test_df_sqrt, file = "eli_tests.Rdata")

openxlsx::write.xlsx(x = eli_test_df_sqrt, file = "ELI_TFsonly_pvals.xlsx", asTable = T)


# Look deeper in some TFs for ELI genes

all_eli_TF_gene_list= list()
## NOTE: this takes a long time
for (mtf in eli_test_df_sqrt[padj < .05]$TFBS) {
  if (!(mtf %in% names(all_eli_TF_gene_list))) {
    all_rows = lapply(X = unique(as.character(eli_exp_counts$hsapiens_gene_id)), FUN = function(gg) {
      ff = as.formula(paste0("sqrtcorrVal ~ ",mtf))
      sfit = summary(lm(data = eli_exp_counts[hsapiens_gene_id==gg], formula = ff))
      list(hsapiens_gene_id=gg, motif=mtf, 
           p = if (nrow(sfit$coefficients)<2) NA else sfit$coefficients[2,4],
           beta = if (nrow(sfit$coefficients)<2) NA else sfit$coefficients[2,1], 
           r2 = sfit$r.squared)
    })
    all_eli_TF_gene_list[[mtf]] <- as.data.frame(do.call(what = rbind, args = all_rows))
    all_eli_TF_gene_list[[mtf]]$hsapiens_gene_id = unlist(all_eli_TF_gene_list[[mtf]]$hsapiens_gene_id)
    all_eli_TF_gene_list[[mtf]]$motif = unlist(all_eli_TF_gene_list[[mtf]]$motif)
    all_eli_TF_gene_list[[mtf]]$p = unlist(all_eli_TF_gene_list[[mtf]]$p)
    all_eli_TF_gene_list[[mtf]]$beta = unlist(all_eli_TF_gene_list[[mtf]]$beta)
    all_eli_TF_gene_list[[mtf]]$r2 = unlist(all_eli_TF_gene_list[[mtf]]$r2)
    all_eli_TF_gene_list[[mtf]] = as.data.table(all_eli_TF_gene_list[[mtf]])
    all_eli_TF_gene_list[[mtf]] <- merge(all_eli_TF_gene_list[[mtf]], namedf2, by.x = "hsapiens_gene_id", by.y="ensembl_gene_id")
    all_eli_TF_gene_list[[mtf]]$padj  <- p.adjust(all_eli_TF_gene_list[[mtf]]$p, method = "fdr")
    cat(paste0("done ", mtf))
  }
}

tfnum=1
all_eli_TF_gene_list[[tfnum]]$label <- ""
all_eli_TF_gene_list[[tfnum]][(padj < exp(-10)) & (abs(beta) > 5), label := hgnc_symbol ]
all_eli_TF_gene_list[[tfnum]][(padj < exp(-7)) & (beta < -3), label := hgnc_symbol ]
all_eli_TF_gene_list[[tfnum]][(padj < exp(-20)) , label := hgnc_symbol ]
tiff(filename = paste0("figures/ELI_tests_", names(all_eli_TF_gene_list)[tfnum], ".tiff"), 
     width = 4, height = 3, units = "in", res = 600, compression = "zip")
ggplot(data = all_eli_TF_gene_list[[tfnum]], aes(x=beta,y=-log(padj))) + geom_point(alpha = .2, shape=20) + 
  geom_hline(yintercept =-log(.05), color="red") + 
  labs(title = names(all_eli_TF_gene_list)[[tfnum]], y = "log(FDR)") + geom_label_repel(aes(label=label), min.segment.length = 0) +
  theme_classic()
dev.off()



# make and write all data
for(tfnum in c(1:length(all_eli_TF_gene_list))) {
  if ("label" %in% names(all_eli_TF_gene_list[[tfnum]]))
    all_eli_TF_gene_list[[tfnum]]$label <- NULL
}
wb = openxlsx::createWorkbook()
all_eli_TFgenecombined <- do.call(what = rbind, args = all_eli_TF_gene_list)
openxlsx::addWorksheet(wb = wb, sheetName = "allTFgenes")
openxlsx::writeDataTable(wb = wb, sheet = "allTFgenes", x = all_eli_TFgenecombined)
eli_TF_gene_summary <- all_eli_TFgenecombined
eli_TF_gene_summary[, numSig := sum(padj < .05, na.rm = T), by=hsapiens_gene_id]
eli_TF_gene_summary <- unique(eli_TF_gene_summary[, .(hgnc_symbol, hsapiens_gene_id, numSig)])
openxlsx::addWorksheet(wb = wb, sheetName = "gene_summry")
openxlsx::writeDataTable(wb = wb, x = eli_TF_gene_summary, sheet = "gene_summry")
eli_TF_gene_TFsumm <- all_eli_TFgenecombined
eli_TF_gene_TFsumm <- merge(eli_TF_gene_TFsumm, unique(eli_scores[, .(hsapiens_gene_id, ELIcor)]), by="hsapiens_gene_id")
eli_TF_gene_TFsumm <- eli_TF_gene_TFsumm[, TFnumSig := sum(padj < .05, na.rm = T), by = motif]
eli_TF_gene_TFsumm <- eli_TF_gene_TFsumm[, TFmeanSigBeta := mean(beta[padj < .05], na.rm = T), by = motif]
eli_TF_gene_TFsumm <- eli_TF_gene_TFsumm[, TFmeanSigELICor := mean(ELIcor[padj < .05], na.rm = T), by = motif]
eli_TF_gene_TFsumm <- unique(eli_TF_gene_TFsumm[, .(motif, TFnumSig, TFmeanSigBeta, TFmeanSigELICor)])
openxlsx::addWorksheet(wb = wb, sheetName = "TF_summry")
openxlsx::writeDataTable(wb = wb, x = eli_TF_gene_TFsumm, sheet = "TF_summry")

openxlsx::saveWorkbook(wb = wb , file = "ELI_TF_gene_combined.xlsx", overwrite = T)


