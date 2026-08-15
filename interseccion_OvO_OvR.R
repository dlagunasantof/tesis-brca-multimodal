# =============================================================================
# Intersección de listas de marcadores: one-vs-one (OvO)  ∩  one-vs-others (OvR)
#
# Cada catálogo asigna a cada gen una etiqueta (Subtype, Direction). El cruce
# identifica los genes en los que AMBOS métodos coinciden -> firma consenso de
# marcadores por subtipo, robusta a las dos definiciones de contraste.
#
# Se calculan tres niveles de exigencia para el match, porque el número de genes
# consenso depende de cuál se use (útil para reproducir la firma original):
#   (1) solo Gene
#   (2) Gene + Subtype
#   (3) Gene + Subtype + Direction   <- el más estricto y defendible
# =============================================================================

# ------------------------------- Parámetros ---------------------------------
ovo_raw <- read.csv("C:/Users/danil/OneDrive/Documentos/tesis-brca-multimodal/resultados_552_one_vs_one/genes_unicos_OvO.tsv", sep="\t")
ovr_raw <- read.csv("C:/Users/danil/OneDrive/Documentos/tesis-brca-multimodal/resultados_552_one_vs_others/unique_genes_catalog_TOP.tsv", sep="\t")
#OUT_DIR <- "resultados_interseccion"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


# =============================================================================
# 1) Cargar y estandarizar a Gene / Subtype / Direction
# =============================================================================
prep <- function(df, etiqueta) {
  req <- c("Gene", "Subtype", "Direction")
  stopifnot(all(req %in% colnames(df)))
  df <- unique(df[, req, drop = FALSE])
  df$Metodo     <- etiqueta
  df$key_gene   <- df$Gene
  df$key_subt   <- paste(df$Gene, df$Subtype, sep = "|")
  df$key_triple <- paste(df$Gene, df$Subtype, df$Direction, sep = "|")
  df
}

ovo <- prep(ovo_raw, "OvO")
ovr <- prep(ovr_raw, "OvR")

cat("OvO:", nrow(ovo), "filas,", length(unique(ovo$Gene)), "genes distintos\n")
cat("OvR:", nrow(ovr), "filas,", length(unique(ovr$Gene)), "genes distintos\n")

# =============================================================================
# 2) Intersección a los tres niveles
# =============================================================================
inter_gene   <- intersect(ovo$key_gene,   ovr$key_gene)
inter_subt   <- intersect(ovo$key_subt,   ovr$key_subt)
inter_triple <- intersect(ovo$key_triple, ovr$key_triple)

cat("\n--- Genes consenso (en ambos métodos) ---\n")
cat("(1) solo Gene            :", length(inter_gene),   "\n")
cat("(2) Gene + Subtype       :", length(inter_subt),   "\n")
cat("(3) Gene + Subtype + Dir :", length(inter_triple), "\n")

# =============================================================================
# 3) Firma consenso concordante (nivel 3) — la recomendada
# =============================================================================
consenso <- ovo[ovo$key_triple %in% inter_triple, c("Gene", "Subtype", "Direction")]
consenso <- consenso[order(consenso$Subtype, consenso$Direction, consenso$Gene), ]
cat("\nConsenso concordante por (Subtype, Direction):\n")
print(table(consenso$Subtype, consenso$Direction))

# =============================================================================
# 4) Discordantes: gen presente en ambos, pero SIN triple compartido
#    (los métodos lo asignan a distinto subtipo o dirección)
# =============================================================================
genes_ambos      <- inter_gene
genes_concordan  <- unique(ovo$Gene[ovo$key_triple %in% inter_triple])
genes_discordan  <- setdiff(genes_ambos, genes_concordan)

discordantes <- merge(
  ovo[ovo$Gene %in% genes_discordan, c("Gene", "Subtype", "Direction")],
  ovr[ovr$Gene %in% genes_discordan, c("Gene", "Subtype", "Direction")],
  by = "Gene", suffixes = c("_OvO", "_OvR")
)
cat("\nGenes en ambos pero discordantes (distinto subtipo/dirección):",
    length(genes_discordan), "\n")

# =============================================================================
# 5) Únicos de cada método (a nivel de gen)
# =============================================================================
unicos_OvO <- setdiff(ovo$Gene, ovr$Gene)
unicos_OvR <- setdiff(ovr$Gene, ovo$Gene)
cat("Únicos de OvO:", length(unicos_OvO),
    " | Únicos de OvR:", length(unicos_OvR), "\n")

# =============================================================================
# 6) Guardar salidas
# =============================================================================
write.table(consenso, file.path(OUT_DIR, "firma_consenso_concordante.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(Gene = inter_gene),
            file.path(OUT_DIR, "consenso_nivel_gen.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(discordantes, file.path(OUT_DIR, "genes_discordantes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(Gene = unicos_OvO), file.path(OUT_DIR, "unicos_OvO.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(Gene = unicos_OvR), file.path(OUT_DIR, "unicos_OvR.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nSalidas en", OUT_DIR, ":\n",
    "- firma_consenso_concordante.tsv  (Gene+Subtype+Direction en ambos) <- la firma\n",
    "- consenso_nivel_gen.tsv          (solo coincidencia de gen)\n",
    "- genes_discordantes.tsv          (en ambos, con etiqueta distinta)\n",
    "- unicos_OvO.tsv / unicos_OvR.tsv\n")


# =============================================================================
# 6) Guardar unión de las dos listas sin duplicados
# =============================================================================
genes_union <- unique(c(ovo$Gene, ovr$Gene))
length(genes_union)

write.table(data.frame(Gene = genes_union),
            file.path(OUT_DIR, "genes_union_OvO_OvR.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# =============================================================================
# 6) Guardar salidas de genes de unión con su expresión (proviene de los analisis previos de procesamiento)
# =============================================================================
expr_union <- log2_filtrado[rownames(log2_filtrado) %in% genes_union, , drop = FALSE]
dim(expr_union)


mapa <- read.delim("C:/Users/danil/OneDrive/Documentos/tesis-brca-multimodal/resultados_552_one_vs_one/mapa_muestras.tsv",
                   sep = "\t", stringsAsFactors = FALSE)

barcodes <- mapa$TCGA_Barcode[match(colnames(expr_union), mapa$Sample_ID)]

salida <- rbind(
  c("IDPaciente", barcodes),
  cbind(rownames(expr_union), as.data.frame(expr_union, check.names = FALSE))
)



subtipos <- sub("\\.[0-9]+$", "", colnames(expr_union))          # BRCA_LumA, BRCA_LumB, ...
barcodes <- mapa$TCGA_Barcode[match(colnames(expr_union), mapa$Sample_ID)]

# --- Armar la tabla de salida: 2 filas de encabezado + expresión ---
salida <- rbind(
  c("Subtipo",    subtipos),
  c("IDPaciente", barcodes),
  cbind(rownames(expr_union), as.data.frame(expr_union, check.names = FALSE))
)

# --- Escribir el archivo ---
write.table(
  salida, file.path(OUT_DIR, "Consensus_DEG_signature.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)


# =============================================================================
# 7) Tabla resumen: genes por (Subtipo, Dirección) en OvR, OvO e intersección
# =============================================================================
# Cuenta genes distintos por combinación Subtype+Direction en cada lista
conteo <- function(df, nombre) {
  agg <- aggregate(Gene ~ Subtype + Direction, data = df,
                   FUN = function(g) length(unique(g)))
  names(agg)[3] <- nombre
  agg
}

tab_ovr <- conteo(ovr,      "Genes_OvR")
tab_ovo <- conteo(ovo,      "Genes_OvO")
tab_int <- conteo(consenso, "Genes_Interseccion")   # consenso = intersección nivel triple

# Unir las tres columnas (todas las combinaciones; 0 donde falte)
resumen <- Reduce(function(a, b) merge(a, b, by = c("Subtype", "Direction"), all = TRUE),
                  list(tab_ovr, tab_ovo, tab_int))
resumen[is.na(resumen)] <- 0

# Etiqueta "UP_LumA" y orden (LumA, LumB, Basal) x (UP, DOWN)
resumen$Subtipo    <- paste(resumen$Direction, resumen$Subtype, sep = "_")
resumen$Subtype    <- factor(resumen$Subtype,    levels = c("LumA", "LumB", "Basal"))
resumen$Direction  <- factor(resumen$Direction,  levels = c("UP", "DOWN"))
resumen <- resumen[order(resumen$Subtype, resumen$Direction), ]
resumen <- resumen[, c("Subtipo", "Genes_OvR", "Genes_OvO", "Genes_Interseccion")]

cat("\n--- Tabla resumen por subtipo y dirección ---\n")
print(resumen, row.names = FALSE)

write.table(resumen, file.path(OUT_DIR, "tabla_resumen_subtipo_direccion.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)