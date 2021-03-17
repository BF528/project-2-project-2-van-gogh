df_upreg <- data.frame(Term= "Contractile Fiber", Enrich_s =6.99)
df_upreg <- rbind(df_upreg, data.frame(Term= "Extracellular Organelle", Enrich_s =10.79))
df_upreg <- rbind(df_upreg, data.frame(Term= "Organic Acid Metabolic Process", Enrich_s =14.57))
df_upreg <- rbind(df_upreg, data.frame(Term= "ATP/Nucleotide Metabolic Process", Enrich_s =15.29))
df_upreg <- rbind(df_upreg, data.frame(Term= "Mitocondrion", Enrich_s =21.34))

df_dnreg <- data.frame(Term= "Organ Morphogenesis", Enrich_s =8.08)
df_dnreg <- rbind(df_dnreg, data.frame(Term= "Regulation of Cell Component", Enrich_s =8.63))
df_dnreg <- rbind(df_dnreg, data.frame(Term= "Cell Proliferation", Enrich_s =9.16))
df_dnreg <- rbind(df_dnreg, data.frame(Term= "Extracellular Matrix", Enrich_s =9.69))
df_dnreg <- rbind(df_dnreg, data.frame(Term= "Cell Cycle", Enrich_s =11.85))

png(filename="up.png", width=1000, bg="white")
par(mar=c(5,25,4,10)+.1)
barplot(df_upreg$Enrich_s, 
        main = "Up-Regulated", 
        xlab = "Enrichment Score", 
        names.arg = df_upreg$Term,
        cex.names = 1.5,
        cex.main = 1.5,
        cex.lab = 1.5,
        horiz = TRUE, las=1)
dev.off()

png(filename="down.png", width=1000, bg="white")
par(mar=c(5,25,4,10)+.1)
barplot(df_dnreg$Enrich_s, 
        main = "Down-Regulated", 
        xlab = "Enrichment Score", 
        names.arg = df_dnreg$Term,
        cex.names = 1.5,
        cex.main = 1.5,
        cex.lab = 1.5,
        horiz = TRUE, las=1)
dev.off()

