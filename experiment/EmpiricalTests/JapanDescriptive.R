# This script reproduce table 12 in the appendix
# Descriptive statistics for Japanses producer revenue share of intermediate material
library(stargazer)
library(ggplot2)
library(reshape)
library(NormalRegPanelMixture)
library(foreign)
library(Hmisc)


ind_list <- c("food","textile", "wood","paper", "chemical", 
              "petro","plastic","ceramics","steel","othermetal",
              "metal product","machine","electronics",
              "transportation equipment","precision instrument",
              "other")

df <- read.csv(file="data/data_production_function_missing2zero_ver3.csv")
df <- df[df$lnmY_it > -3,]
df <- df[df$lnmY_it < log(2),]

df$t <- df$year
df <- df[order(df$id,df$t),]
df <- df[df$industry_2!=0,]
df <- df[!is.na(df$industry_2),]
ind.code <- unique(df$industry_2)

T.cap <- max(df$year)

T = 5
t.start <- T.cap-T+1
t.seq <- seq(from=t.start,to=t.start+T-1)

df = df[df$year >= t.seq[1],]
# m.share <- cast(df,id ~ year,value="lnmY_it")#Collapse the dataframe into panel form , year against firm id
# row.names(m.share) <- m.share$id 
# m.share <- m.share[,!(colnames(m.share)=="id")] 
# m.share <- m.share[complete.cases(m.share),]
# id.list <- row.names(m.share)
# # ind.code <- ind.code[1:3] #For test purpose, only use the first three industries
# df <- df[df$id %in% id.list,]

desc.table = matrix(nc=8,nr=length(ind.code))
count = 0
df.include <- data.frame()

for (each.code in ind.code){
  count = count + 1
  ind.each <- subset(df,industry_2==each.code)
  ind.name <- capitalize(ind_list[each.code])
  
  
  m.share <- cast(ind.each,id ~ year,value="lnmY_it")#Collapse the dataframe into panel form , year against firm id
  row.names(m.share) <- m.share$id 
  m.share <- m.share[,!(colnames(m.share)=="id")] 
  m.share <- m.share[complete.cases(m.share),]
  id.list <- row.names(m.share)
  ind.each <- ind.each[ind.each$id %in% id.list,]
  df.include <- rbind(df.include,ind.each)
  
  desc.table[count, ] <- c(each.code, ind.name,dim(ind.each)[1],dim(m.share)[1],round(mean(exp(ind.each$lnmY_it)),2),round(sd(exp(ind.each$lnmY_it)),2),round(mean(ind.each$y_it),2),round(sd(ind.each$y_it),2))
  
}

colnames(desc.table) <- c("code","Industry","NObs","N","Mean_m_share","Sd_m_share","Mean_lnY","Sd_lnY")


desc.table <- transform(desc.table, N = as.numeric(N))
desc.table <- desc.table[order(desc.table[,'N'],decreasing = TRUE),]  
desc.table <- transform(desc.table, N = as.character(N))

sink("results/Japan/desc.table.txt")
stargazer(desc.table,type="latex",summary=FALSE, title="Descriptive statistics for Japanses producer revenue share of intermediate material",rownames=FALSE)
sink()

saveRDS(df.include, file = "data/JapanClean.rds")
