library(DESeq2)
library(ggplot2)

setwd("/Users/marigutierrez/Desktop/transcriptomics/INPUTS")

# 1) Cargar datos
counts   <- read.delim("E-ENAD-3-raw-counts.tsv")
metadata <- read.delim("E-ENAD-3-experiment-design.tsv")

# 2) Preparar tabla de conteos
rownames(counts) <- counts$Gene.ID
genes <- counts[, c("Gene.ID", "Gene.Name")]
counts <- counts[, -c(1, 2)]

# 3) Preparar metadata
rownames(metadata) <- metadata$Run

metadata <- metadata[, c(
  "Sample.Characteristic.organism.part.",
  "Factor.Value.environmental.stress."
)]

colnames(metadata) <- c("organ", "stress")

# Limpiar y convertir a factor ANTES de filtrar
metadata$organ  <- factor(gsub(" ", "_", metadata$organ))
metadata$stress <- factor(gsub(" ", "_", metadata$stress))

# AÑADIR CAMBIO CRUCIAL: Definir explícitamente 'stress' como el nivel de referencia (denominador).
# Esto invertirá el LFC predeterminado. Ahora: LFC = log2(Expresión_none / Expresión_drought)
metadata$stress <- relevel(metadata$stress, ref = "drought_environment")

# 4) Filtrar SOLO raíces
metadata_root <- metadata[metadata$organ == "root", ]
counts_root   <- counts[, rownames(metadata_root)]

# Verificar emparejamiento
stopifnot(all(colnames(counts_root) == rownames(metadata_root)))

# 5) Crear DESeqDataSet
dds_root <- DESeqDataSetFromMatrix(
  countData = counts_root,
  colData   = metadata_root,
  design    = ~ stress
)

dds_root

# 5b) Ejecutar DESeq2
dds_root <- DESeq(dds_root)

# 5c) Obtener resultados (drought vs none)
# Opcional pero recomendado: especificar el contraste explícitamente para mayor claridad
# Esto fuerza a comparar 'none' (control) contra 'drought' (referencia)
res_root <- results(dds_root, contrast=c("stress", "drought_environment", "none"))

# Convertir a data frame
res_df <- as.data.frame(res_root)
res_df$gene <- rownames(res_df)

# Quitar NA para evitar errores en el plot
res_df <- res_df[!is.na(res_df$padj), ]


# 6) Volcano plot
# ==============================================================================
res_df$threshold <- "NS"
res_df$threshold[res_df$padj < 0.05 & res_df$log2FoldChange > 1]  <- "Up"
res_df$threshold[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down"

# Crear volcano plot y guardarlo en un objeto
p_volcano <- ggplot(res_df,
                    aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano plot — Roots (drought vs none)",
    x = "log2 Fold Change",
    y = "-log10 adjusted p-value"
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed")
print(p_volcano)
ggsave("volcano_root_drought_vs_control.pdf", plot = p_volcano, width = 8, height = 6)
ggsave("volcano_root_drought_vs_control.png", plot = p_volcano, width = 8, height = 6, dpi = 300)


# === GUARDAR RESULTADOS DE DESEQ2 ===

# Guardar tabla completa como CSV
write.csv(res_df, "DESeq2_results_root_drought_vs_control.csv", row.names = FALSE)

# Guardar solo genes significativos
res_sig <- res_df[res_df$padj < 0.05, ]
write.csv(res_sig, "DESeq2_results_root_drought_vs_control_SIGNIFICANT.csv", row.names = FALSE)

