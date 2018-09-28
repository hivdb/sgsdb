library(ggplot2)
library(grid)
library(formattable)
library(gridExtra)

CATEGORIES = c('All', 'SubtypeB', 'SubtypeC', 'Non-SubtypeBC', 'Naive', 'ART')

SGS_PR=read.csv('data/prevalence/SGS.PRprevalence.csv')
SGS_RT=read.csv('data/prevalence/SGS.RTprevalence.csv')
SGS_IN=read.csv('data/prevalence/SGS.INprevalence.csv')

CONSENSUS=read.csv('data/consensus.csv')

GENES = c('PR', 'RT', 'IN')

breaks = 10**(1:5)/1000


for (cat in CATEGORIES) {

  plots = list()
  DB_PREVALENCE_DATA = read.csv(
    sprintf('data/prevalence/dbAminoAcidVariants%s.csv', cat)
  )

  for (gene in GENES) {
    sgsData = switch(
      gene,
      PR = SGS_PR,
      RT = SGS_RT,
      IN = SGS_IN
    )
    cons = switch(
      gene,
      PR = as.character(CONSENSUS[1,2]),
      RT = as.character(CONSENSUS[2,2]),
      IN = as.character(CONSENSUS[3,2])
    )
    sgsData = sgsData[sgsData$Category == cat,]
    sgsData$sgsPcnt = sgsData$Pcnt
    dbData = DB_PREVALENCE_DATA[DB_PREVALENCE_DATA$gene == gene,]
    dbData$dbPcnt = dbData$percent * 100
    data = merge(sgsData, dbData, by.x=c("Pos", "AA"), by.y=c("position", "aa"))
    data$Pos = as.numeric(as.character(data$Pos))
    data$Cons = substring(cons, data$Pos, data$Pos)
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
      PR = sprintf('data/prevalence/CompPR%s.csv', cat),
      RT = sprintf('data/prevalence/CompRT%s.csv', cat),
      IN = sprintf('data/prevalence/CompIN%s.csv', cat)
    )
    out = data[data$sgsPcnt > 0 & data$AA != data$Cons,]
    # out = data[data$sgsPcnt > 0,]
    write.csv(subset(out, select=c("Pos", "Cons", "AA", "Count", "PosTotal", "PatientCount",
                                   "PatientPosTotal", "sgsPcnt", "dbPcnt", "pcntFold", "isAPOBEC", "isUsual")
                     ), target, row.names=FALSE)
  }
  pdf(sprintf('data/%sPrevalanceCmp.pdf', cat), width=8, height=24)
  print(grid.arrange(grobs = plots, ncol = 1, nrow = 3))
  dev.off()
}
