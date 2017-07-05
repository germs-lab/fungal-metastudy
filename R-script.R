
### Prepare CSV file for co-occurence analysis

data <- read.table("FUNGAL_METASTUDY_TABLES_SOIL/ITS_EXTRACTED_SAMPLES_VS_SH_TABLE_NORMALISED_SOIL_NOSINGLETONS.txt", header=1, row.name=1)

tdata <- t(data)
fdata <- cbind(trt="none", tdata)

write.csv(fdata, "tomas_table_for_co.csv")


### Read data after calculation then group
library(igraph)
results <- read.table("tomas_result_p0.05_r0.6.tsv")
names(results)<-c("trt","taxa1","taxa2","rho","p.value")
gra <- graph.edgelist(as.matrix(results[,c(2:3)]), directed=F)
fgc <- cluster_fast_greedy(gra)

#save member
mem <-membership(fgc)
mamem <- sapply(mem, matrix)
write.csv(mamem,file="member.csv")

#save groups
library(plyr)
gr <- groups(fgc)
ma <- ldply(gr, rbind)
write.csv(ma,file="group.csv")
