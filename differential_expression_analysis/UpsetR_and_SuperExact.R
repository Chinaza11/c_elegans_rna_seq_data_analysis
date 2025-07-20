setwd("G:/My Drive/PhD/project/project_with_Dr.Reid/c.elegans_RNAseq_suppl/data_analysis/data_and_result/differential_expression_analysis/")
# =============================================================================
# UpSetR plot
# =============================================================================

library(UpSetR)

tpmg = read.table("deseq2/tpmg_up_and_down.txt", header=TRUE, sep="\t")
tyfd = read.table("deseq2/tyfd_up_and_down.txt", header=TRUE, sep="\t")
tyex = read.table("deseq2/tyex_up_and_down.txt", header=TRUE, sep="\t")

listInput = list(TPMG = rownames(tpmg), TYFD = rownames(tyfd), TYEX = rownames(tyex))

png("UpsetR_and_SuperExact/UpSetR_plot.png", width=12, height=12, units="in", res=500)
upset(fromList(listInput), order.by="freq", mainbar.y.label="DEGs Intersections", sets.x.label="DEGs Per Treatment", text.scale=1.25, sets.bar.color=c('red','blue','green'))
dev.off()

# =============================================================================
# SuperExact test and Circular Venn diagram
# =============================================================================

library(SuperExactTest)
library(ggplot2)

png("UpsetR_and_SuperExact/SuperExactTest_landscape_layout.png", width=8, height=8, units="in", res=500)
obj=supertest(listInput,n=19257)
plot(obj, Layout="landscape", sort.by="size", title='Landscape layout', cex=1, cex.title=2.5)
dev.off()

#png("UpsetR_and_SuperExact/SuperExactTest_circular_layout.png", width=8, height=8, units="in", res=500) #wasn't saving it properly and I had to save it through the plot window, TIFF at 1000 x 1000

plot(obj, Layout="circular", sort.by="size", title='Circular layout', cex.title=2, legend.pos=c(0.9,0.15), track.area.range=0.275, legend.text.cex=1.5)

#dev.off()

summary(obj)

#write.csv(summary(obj)$Table, file="UpsetR_and_SuperExact/SuperExact_test_result.csv", row.names=FALSE)
