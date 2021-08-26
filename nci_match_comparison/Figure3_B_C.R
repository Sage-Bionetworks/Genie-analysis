#### PANEL B DATA
datB <- read.delim('freqs.csv', sep=',')

lm_all <- lm(FreqCountsNew ~ FreqNCI, data=datB)
summary(lm_all)$adj.r.squared

datB_filt <- datB[datB$code != 'Z1F',]
lm_filt <- lm(FreqCountsNew ~ FreqNCI, data=datB_filt)
summary(lm_filt)$adj.r.squared


#### PANEL C DATA
## GENIE data
datC_raw <- read.delim('data_clinical_sample.txt', as.is=T, skip=4) # file from synapse

# calculate patient-level cancer type frequency
datC_pat <- unique(datC_raw[,c('PATIENT_ID','CANCER_TYPE')]) # cleans up patients with multiple samples of the same cancer type
datC_pat_freq <- round(100 * table(datC_pat$CANCER_TYPE) / length(datC_pat$CANCER_TYPE), 2)
datC_pat_freq_df <- data.frame(tail(sort(datC_pat_freq), n=8))
names(datC_pat_freq_df) <- c('CANCER_TYPE','GENIE_FREQ')


## NCI-MATCH DATA
# data from Table 2 of Flaherty et al JCO 2020 PMID:33048619
dat_match <- data.frame(
  CANCER_TYPE = c('Colorectal Cancer', 'Breast Cancer', 'Non-Small Cell Lung Cancer', 'Prostate Cancer', 'Ovarian Cancer','Pancreatic Cancer','Glioma','Melanoma'),
  MATCH_COUNT = c(848, 685, 407, 139, 530, 344, 94, 76)
)
match_total = 5540
dat_match$MATCH_FREQ <- round(dat_match$MATCH_COUNT / match_total * 100, 2)

## MERGE GENIE & NCI-MATCH
datC_final <- merge(datC_pat_freq_df, dat_match)
datC_final <- datC_final[order(datC_final$GENIE_FREQ, decreasing=T),]

datC_plot <- datC_final[,c('GENIE_FREQ','MATCH_FREQ')]
colnames(datC_plot) <- c("GENIE", "NCI-MATCH")
rownames(datC_plot) <- datC_final$CANCER_TYPE



#### PLOT
## SET UP PDF
pdf('genie_clinical_trials_BC.pdf', width=12, height=12, useDingbats=FALSE)
par(mfrow=c(2,2)); par(pty='s')

## PANEL B MAIN PLOT
plot(FreqCountsNew ~ FreqNCI, data=datB, xlim=c(0,8), ylim=c(0,8), pch=19, cex=.5, xlab='% Patients Matching in NCI-MATCH', ylab='% Patients Matching in GENIE', las=1)
text(FreqCountsNew ~ FreqNCI, data=datB, labels=datB$code, cex=0.6, pos=4)
abline(lm_all, col='blue')
abline(v=1); abline(h=1)  ## to mark inset region
#abline(0,1, col='red')

## PANEL B INSET
plot(FreqCountsNew ~ FreqNCI, data=datB, xlim=c(0,1), ylim=c(0,1), pch=19, xlab='', ylab='', las=1, xaxp=c(0,1,2), yaxp=c(0,1,2))
text(FreqCountsNew ~ FreqNCI, data=datB, labels=datB$code, cex=0.6, pos=4)
abline(lm_all, col='blue')

## PANEL C PLOT
par(mar=c(10,4,4,2)); barplot(t(datC_plot), beside=T, legend=T, ylim=c(0,16), las=2, cex.names=0.8, col=c("#7F7F7F","#CCCCCC"), ylab="Cancer Type Frequency"); box(bty='l')

#### CLOSE PDF
dev.off()

