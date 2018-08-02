library(ggplot2)
library(grid)
library(formattable)
library(gridExtra)

SUBTYPE = 'SubtypeB'
SHORT_SUBTYPE = 'B'

# SUBTYPE = 'All'
# SHORT_SUBTYPE = ''

DB_PREVALENCE_DATA = read.csv(
  sprintf(
    '~/Dropbox/NGSPcnts/SGS/prevalenceData/dbAminoAcidVariants%s.csv',
    SUBTYPE
  )
)

SGS_PR=read.csv('~/Dropbox/NGSPcnts/SGS/prevalenceData/SGS.PRprevalence.csv')
SGS_RT=read.csv('~/Dropbox/NGSPcnts/SGS/prevalenceData/SGS.RTprevalence.csv')
SGS_IN=read.csv('~/Dropbox/NGSPcnts/SGS/prevalenceData/SGS.INprevalence.csv')

COMP_PR = sprintf(
  '~/Dropbox/NGSPcnts/SGS/prevalenceData/CompPR%s.csv',
  SUBTYPE
)
COMP_RT = sprintf(
  '~/Dropbox/NGSPcnts/SGS/prevalenceData/CompRT%s.csv',
  SUBTYPE
)
COMP_IN = sprintf(
  '~/Dropbox/NGSPcnts/SGS/prevalenceData/CompIN%s.csv',
  SUBTYPE
)

GENES = c('PR', 'RT', 'IN')

plots = list()

breaks = 10**(1:5)/1000

for (gene in GENES) {
  sgsData = switch(
    gene,
    PR = SGS_PR,
    RT = SGS_RT,
    IN = SGS_IN
  )
  sgsData = sgsData[sgsData$Subtype == SHORT_SUBTYPE,]
  sgsData$sgsPcnt = sgsData$Pcnt
  dbData = DB_PREVALENCE_DATA[DB_PREVALENCE_DATA$gene == gene,]
  dbData$dbPcnt = dbData$percent * 100
  data = merge(sgsData, dbData, by.x=c("Pos", "AA"), by.y=c("position", "aa"))
  if (gene == "RT") {
    data = data[data$Pos < 241,]
  }
  data$pcntFold = data$sgsPcnt / data$dbPcnt
  # data = c(data$position, data$aa, data$dbPcnt, data$sgsPcnt)
  df = subset(data, select=c("sgsPcnt", "dbPcnt"))
  model <- lm(dbPcnt ~ sgsPcnt, data=df)
  simx = seq(0, 100, step=0.01)
  confints = predict(model, newdata=data.frame(sgsPcnt=simx),
                      interval='confidence', level=0.95)
  confints = as.data.frame(confints)
  confints$lwr[confints$lwr < 0] = 0
  predints = predict(model, newdata=data.frame(sgsPcnt=simx),
                     interval='prediction', level=0.9)
  predints = as.data.frame(predints)
  predints$lwr[predints$lwr < 0] = 0
  
  #confints <- with(df, data.frame(lwr=))
  plot = ggplot(cbind(cbind(x=simx, confints))) +
    # geom_abline(aes(x=sgsPcnt, y=dbPcnt), color="blue") +
    geom_ribbon(aes(ymin=lwr, ymax=upr, x=x), fill="blue", alpha=0.2) +
    geom_line(aes(y=fit, x=x), color="blue", alpha=0.8) +
    geom_ribbon(aes(ymin=lwr, ymax=upr, x=x), data=cbind(cbind(x=simx, predints)), fill="blue", alpha=0.05) +
    geom_point(aes(x=sgsPcnt, y=dbPcnt), data=df, size=0.3) +
    scale_x_log10("SGSPcnt", breaks=breaks, labels=comma(breaks), limits=c(0.01, 100)) +
    scale_y_log10("DBPcnt", breaks=breaks, labels=comma(breaks), limits=c(0.01, 100)) +
    ggtitle(gene)
  plots = append(plots, list(plot))
  target = switch(
    gene,
    PR = COMP_PR,
    RT = COMP_RT,
    IN = COMP_IN
  )
  out = data[data$sgsPcnt > 0 & data$dbPcnt < 60,]
  write.csv(subset(out, select=c("Pos", "AA", "Count", "PosTotal", "PatientCount",
                                 "PatientPosTotal", "sgsPcnt", "dbPcnt", "pcntFold", "isAPOBEC")
                   ), target, row.names=FALSE)
}
pdf(sprintf('~/Dropbox/NGSPcnts/SGS/prevalenceData/%sCmp.pdf', SUBTYPE), width=8, height=24)
print(grid.arrange(grobs = plots, ncol = 1, nrow = 3))
dev.off()