# This script is used to reproduce table 13, the descriptive statistics for Chilean industries
library(stargazer)
library(ggplot2)
library(reshape)
library(NormalRegPanelMixture)
library(foreign)
library(haven)

df <- read_dta("data/ChileanClean.dta")
ind.code <- na.omit(unique(df$ciiu_3d))

T.cap <- max(df$year)

T = 5
t.start <- T.cap-T+1
t.seq <- seq(from=t.start,to=t.start+T-1)

df = df[df$year >= t.seq[1],]
# df = df[complete.cases(df),]
# m.share <- cast(df,id ~ year,value="si")#Collapse the dataframe into panel form , year against firm id
# row.names(m.share) <- m.share$id 
# m.share <- m.share[,!(colnames(m.share)=="id")] 
# m.share <- m.share[complete.cases(m.share),]
# id.list <- row.names(m.share)
# # ind.code <- ind.code[1:3] #For test purpose, only use the first three industries
# df <- df[df$id %in% id.list,]

df <- df[df$si > -3,]
df <- df[df$si < log(2),]

desc.table = matrix(nc=8,nr=length(ind.code))
count = 0
df.include <- data.frame()
for (each.code in ind.code){
  count = count + 1
  ind.each <- subset(df,ciiu_3d==each.code)
  ind.name <- ind.each$ciiu3d_descr[1]
  ind.each$lny <- log(ind.each$GO)
  ind.each$lnm <- log(ind.each$WI)
  ind.each$lnl <- log(ind.each$L)
  ind.each$lnk <- log(ind.each$K)
  
  
  m.share <- cast(ind.each,id ~ year,value="si")#Collapse the dataframe into panel form , year against firm id
  row.names(m.share) <- m.share$id 
  m.share <- m.share[,!(colnames(m.share)=="id")] 
  m.share <- m.share[complete.cases(m.share),]
  id.list <- row.names(m.share)
  ind.each <- ind.each[ind.each$id %in% id.list,]
  ind.each <- ind.each[ind.each$year >= t.start,]
  df.include <- rbind(df.include,ind.each)
  
  desc.table[count, ] <- c(each.code, ind.name,dim(ind.each)[1],dim(exp(m.share))[1],round(mean(exp(ind.each$si)),2),round(sd(exp(ind.each$si)),2),round(mean(ind.each$lny),2),round(sd(ind.each$lny),2))
  
}
colnames(desc.table) <- c("Code","Industry","NObs", "N","Mean_m_share","Sd_m_share","Mean_lnY","Sd_lnY")

desc.table <- transform(desc.table, N = as.numeric(N))
desc.table <- desc.table[order(desc.table[,'N'],decreasing = TRUE),]  
desc.table <- transform(desc.table, N = as.character(N))

sink("results/Empirical/Chile_desc.table.txt")
stargazer(desc.table,type="latex", summary=FALSE, title="Descriptive statistics for Chilean producer revenue share of intermediate material",rownames=FALSE)
sink()

saveRDS(df.include, file = "data/ChileanClean.rds")
