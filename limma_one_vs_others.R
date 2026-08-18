# =============================================================================
# Análisis diferencial de subtipos BRCA — LIMMA sin voom, contrastes ONE-vs-OTHERS
# (LumA vs resto, LumB vs resto, Basal vs resto)
#
# Entrada : matriz FPKM gen x muestra, con nombres de columna tipo BRCA_LumA.1
# Salida  : genes ÚNICOS por subtipo y dirección (UP/DOWN) sobre el TOP de genes
#           más variables, más el panel de expresión y el catálogo etiquetado
#           que alimentan el análisis de intersección posterior.
#           + PDFs de gráficos (PCA, Volcano, MA y Heatmap) por análisis.
#
# NOTA METODOLÓGICA:
#  - La DE se corre sobre todos los genes filtrados (referencia) y sobre los
#    TOP más variables (análisis principal). Las listas finales usan el TOP.
#    Seleccionar por varianza y luego testear diferencias sobre los mismos datos
#    puede sesgar; es una decisión consciente que debe documentarse en la tesis.
#  - Aquí la unicidad sale de UN contraste por subtipo (subtipo vs promedio de
#    los otros dos). En el script one-vs-one se exige consistencia frente a los
#    DOS contrastes. Son dos definiciones distintas de "marcador"; por eso el
#    cruce entre ambas listas identifica genes robustos a ambos planteamientos.
# =============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)          # <-- AÑADIDO: PCA con ggplot
  library(pheatmap)         # <-- AÑADIDO: heatmaps
  library(RColorBrewer)     # <-- AÑADIDO: paletas
  library(EnhancedVolcano)  # <-- AÑADIDO: volcano plots
})

# ------------------------------- Parámetros ---------------------------------
dat   <- read.csv("C:/Users/danil/OneDrive/Escritorio/Propuesta de Tesis/Analisis python/552_expresion_genes_filtrado.tsv", sep="")  # 
OUT_DIR      <- "resultados_552_one_vs_others"

PADJ_CUT     <- 0.05     # umbral de FDR (BH)
LFC_CUT      <- 1        # umbral de |log2FC|
ZERO_FRAC    <- 0.50     # se descarta un gen con > 50% de ceros
COR_CUT      <- 0.75     # correlación media para señalar outliers
TOPN_VAR     <- 1000     # nº de genes más variables (análisis principal)
TOPN_HEATMAP <- 100      # genes por contraste en los heatmaps

# baja correlación media con el resto (revisado manualmente). c() = no eliminar.
OUTLIERS_TO_REMOVE <- c()

SUBTYPE_LEVELS <- c("LumA", "LumB", "Basal")
SAVE_PLOTS   <- TRUE     # guarda PDFs en OUT_DIR; FALSE -> a pantalla

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1) Carga y conversión a numérico
# =============================================================================

cat("Matriz cruda:", nrow(dat), "genes x", ncol(dat) - 1, "muestras\n")

genes_all <- as.character(dat[[1]])
raw       <- dat[, -1, drop = FALSE]

# --- Separar filas de metadatos embebidas (no numéricas) al inicio del archivo ---
# El TSV trae filas de metadatos ANTES de la expresión: etiqueta de muestra y
# barcode TCGA del paciente. Se detectan como filas iniciales no numéricas.
frac_num <- function(r) mean(!is.na(suppressWarnings(as.numeric(as.character(r)))))
n_meta <- 0L
while (n_meta < nrow(raw) && frac_num(unlist(raw[n_meta + 1L, ])) < 0.5)
  n_meta <- n_meta + 1L
cat("Filas de metadatos detectadas al inicio:", n_meta, "\n")

# Metadatos por muestra (barcode TCGA -> emparejar con las WSI del mismo paciente)
sample_meta <- data.frame(
  Column  = colnames(raw),
  Subtype = sub("^BRCA_", "", sub("\\..*$", "", colnames(raw))),
  stringsAsFactors = FALSE
)
if (n_meta >= 1) sample_meta$Sample_Label <- as.character(unlist(raw[1, ]))
if (n_meta >= 2) sample_meta$TCGA_Barcode <- as.character(unlist(raw[2, ]))

# Matriz de expresión = filas posteriores a los metadatos
idx   <- if (n_meta > 0) -(seq_len(n_meta)) else TRUE
genes <- genes_all[idx]
mat   <- raw[idx, , drop = FALSE]

# make.unique: MISMO criterio que el script one-vs-one, para que los símbolos de
# gen duplicados (pseudoautosómicos X/Y: CD99, SLC25A6...) reciban el mismo
# identificador en ambos y la intersección cruce bien.
rownames(mat) <- make.unique(genes)
mat[] <- lapply(mat, function(x) suppressWarnings(as.numeric(x)))

# =============================================================================
# 2) Filtrado por QC y transformación log2
# =============================================================================
mat <- mat[rowSums(mat, na.rm = TRUE) > 0, , drop = FALSE]   # quitar todo-cero
samples_original <- colnames(mat)

log2_mat <- log2(mat + 1)

ceros  <- rowSums(mat == 0, na.rm = TRUE)                     # filtrar por ceros
umbral <- ncol(mat) * ZERO_FRAC
log2_filtrado <- log2_mat[ceros < umbral, , drop = FALSE]

keep <- apply(log2_filtrado, 1, var, na.rm = TRUE) > 0        # quitar var = 0
log2_filtrado <- log2_filtrado[keep, , drop = FALSE]
cat("Tras filtrado QC:", nrow(log2_filtrado), "genes x",
    ncol(log2_filtrado), "muestras\n")

# =============================================================================
# 3) Detección y eliminación de outliers de muestra
# =============================================================================
log2_cc  <- log2_filtrado[complete.cases(log2_filtrado), , drop = FALSE]
cor_prom <- apply(cor(log2_cc), 1, mean)
cat("\nCorrelación media por muestra (resumen):\n"); print(summary(cor_prom))

posibles_out <- names(which(cor_prom < COR_CUT))
if (length(posibles_out)) {
  cat("Posibles outliers (corr media <", COR_CUT, "):\n"); print(posibles_out)
} else cat("No se detectaron outliers bajo el umbral", COR_CUT, "\n")


# Eliminación (según OUTLIERS_TO_REMOVE, decisión documentada arriba)

outliers_detectados <- union(posibles_out, OUTLIERS_TO_REMOVE)
if (length(outliers_detectados)) {
  log2_filtrado <- log2_filtrado[, !colnames(log2_filtrado) %in% outliers_detectados,
                                 drop = FALSE]
  cat("Eliminadas por correlación <", COR_CUT, ":",
      paste(outliers_detectados, collapse = ", "), "\n")
} else {
  cat("No se eliminó ninguna muestra.\n")
}

samples_postfilter <- colnames(log2_filtrado)
cat("Matriz final:", nrow(log2_filtrado), "genes x",
    ncol(log2_filtrado), "muestras\n")

# =============================================================================
# 4) Definición de subtipos y mapa de trazabilidad
# =============================================================================
subtype_raw <- sub("^BRCA_", "", sub("\\..*$", "", colnames(log2_filtrado)))
subtypes    <- factor(subtype_raw, levels = SUBTYPE_LEVELS)
if (any(is.na(subtypes)))
  stop("Etiquetas de subtipo no reconocidas. Revisa colnames(log2_filtrado).")
cat("\nMuestras por subtipo:\n"); print(table(subtypes))

meta_map <- data.frame(
  Sample_ID    = samples_postfilter,
  Subtype      = as.character(subtypes),
  TCGA_Barcode = if (!is.null(sample_meta$TCGA_Barcode))
    sample_meta$TCGA_Barcode[match(samples_postfilter, sample_meta$Column)] else NA,
  stringsAsFactors = FALSE
)
write.table(meta_map, file.path(OUT_DIR, "mapa_muestras.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Muestras eliminadas respecto al inicio:",
    paste(setdiff(samples_original, samples_postfilter), collapse = ", "), "\n")

# =============================================================================
# 5) Diseño y contrastes one-vs-others
# =============================================================================
design <- model.matrix(~ 0 + subtypes)
colnames(design) <- levels(subtypes)

contr <- makeContrasts(
  LumA_vs_others  = LumA  - (LumB + Basal) / 2,
  LumB_vs_others  = LumB  - (LumA + Basal) / 2,
  Basal_vs_others = Basal - (LumA + LumB)  / 2,
  levels = design
)

# =============================================================================
# 6) Funciones auxiliares
# =============================================================================
fit_limma <- function(expr_mat, design, contr) {
  fit  <- lmFit(expr_mat, design)
  fit2 <- eBayes(contrasts.fit(fit, contr))
  tt_list <- lapply(colnames(contr), function(cf) {
    tt <- topTable(fit2, coef = cf, number = Inf,
                   adjust.method = "BH", sort.by = "P")
    tt$Gene <- rownames(tt)
    tt
  })
  names(tt_list) <- colnames(contr)
  list(fit2 = fit2, tt = tt_list)
}

deg_summary <- function(tt_list, padj_cut = PADJ_CUT, lfc_cut = LFC_CUT) {
  data.frame(
    Contraste   = names(tt_list),
    FDR         = sapply(tt_list, function(tt)
                    sum(tt$adj.P.Val < padj_cut, na.rm = TRUE)),
    FDR_y_logFC = sapply(tt_list, function(tt)
                    sum(tt$adj.P.Val < padj_cut & abs(tt$logFC) > lfc_cut, na.rm = TRUE)),
    row.names = NULL
  )
}

# ------------------------- Funciones de gráficos ----------------------------
# (copiadas del script one-vs-one; son genéricas: iteran sobre names(res$tt),
#  por lo que sirven igual para los contrastes *_vs_others)

# PCA sobre la unión de DEGs de todos los contrastes
pca_union_degs <- function(expr_mat, subtypes, tt_list,
                           padj_cut = PADJ_CUT, lfc_cut = LFC_CUT,
                           titulo = "PCA - Unión DEGs") {
  degs_union <- unique(unlist(lapply(tt_list, function(tt)
    rownames(tt[tt$adj.P.Val < padj_cut & abs(tt$logFC) > lfc_cut, , drop = FALSE]))))
  if (length(degs_union) < 2) {
    message("No hay suficientes DEGs para PCA: ", length(degs_union)); return(invisible())
  }
  X  <- expr_mat[intersect(degs_union, rownames(expr_mat)), , drop = FALSE]
  p  <- prcomp(t(X), scale. = TRUE)
  ve <- 100 * (p$sdev^2) / sum(p$sdev^2)
  df <- data.frame(PC1 = p$x[, 1], PC2 = p$x[, 2], Subtype = subtypes)
  print(
    ggplot(df, aes(PC1, PC2, color = Subtype)) +
      geom_point(size = 2.2) +
      theme_minimal(base_size = 12) +
      xlab(sprintf("PC1 (%.1f%%)", ve[1])) +
      ylab(sprintf("PC2 (%.1f%%)", ve[2])) +
      ggtitle(titulo)
  )
}

# Volcano + MA por contraste. AMBOS usan FDR (adj.P.Val), coherente con las
# listas de genes: lo que se marca en rojo es exactamente lo que entra a los DEGs.
plot_volcano_y_MA <- function(fit2, tt, coef_name, titulo_prefix = "",
                              padj_cut = PADJ_CUT, lfc_cut = LFC_CUT) {
  # Volcano sobre FDR
  print(
    EnhancedVolcano(tt, lab = rownames(tt), x = "logFC", y = "adj.P.Val",
                    title   = paste0(coef_name, " - ", titulo_prefix),
                    ylab    = bquote(~-Log[10]~ 'FDR'),
                    pCutoff = padj_cut, FCcutoff = lfc_cut)
  )
  # MA sobre FDR (p ajustado BH del coeficiente)
  coef_idx <- which(colnames(fit2$coefficients) == coef_name)
  padj <- p.adjust(fit2$p.value[, coef_idx], method = "BH")
  lfc  <- fit2$coefficients[, coef_idx]
  status <- rep(0L, length(padj))
  status[padj < padj_cut & lfc >  lfc_cut] <-  1L
  status[padj < padj_cut & lfc < -lfc_cut] <- -1L
  plotMA(fit2, coef = coef_name, status = status,
         main = sprintf("MA - %s (%s) FDR<%.2f & |logFC|>%g",
                        coef_name, titulo_prefix, padj_cut, lfc_cut))
  abline(h = c(-lfc_cut, lfc_cut), lty = 2, col = "red")
}

# Heatmap de los top genes por contraste
heatmap_top_by_contrast <- function(expr_mat, subtypes, tt_list,
                                    topN = TOPN_HEATMAP, titulo = "Heatmap top sig") {
  top_genes <- unique(unlist(lapply(tt_list, function(tt)
    head(rownames(tt[order(tt$adj.P.Val), , drop = FALSE]), topN))))
  sel <- intersect(top_genes, rownames(expr_mat))
  if (length(sel) < 2) {
    message("No hay suficientes genes para heatmap (n=", length(sel), ")"); return(invisible())
  }
  ann <- data.frame(Subtype = subtypes); rownames(ann) <- colnames(expr_mat)
  pal <- colorRampPalette(c("navy", "white", "firebrick3"))(255)
  pheatmap(expr_mat[sel, , drop = FALSE], scale = "row", show_rownames = FALSE,
           annotation_col = ann, clustering_method = "ward.D2",
           main = titulo, color = pal)
}

# Corre el bloque de gráficos de un análisis, a pantalla o a un PDF
run_plots <- function(expr_mat, subtypes, res, etiqueta) {
  pdf_path <- file.path(OUT_DIR, paste0("plots_", etiqueta, ".pdf"))
  if (SAVE_PLOTS) pdf(pdf_path, width = 8, height = 6) else devAskNewPage(TRUE)
  pca_union_degs(expr_mat, subtypes, res$tt, titulo = paste0("PCA (Unión DEGs) - ", etiqueta))
  for (cf in names(res$tt))
    plot_volcano_y_MA(res$fit2, res$tt[[cf]], cf, titulo_prefix = etiqueta)
  heatmap_top_by_contrast(expr_mat, subtypes, res$tt,
                          titulo = paste0("Heatmap top sig - ", etiqueta))
  if (SAVE_PLOTS) { dev.off(); cat("Gráficos ->", pdf_path, "\n") } else devAskNewPage(FALSE)
}

# =============================================================================
# 7) Referencia — DE sobre todos los genes filtrados (no alimenta las listas)
# =============================================================================
cat("\n===== A) TODOS LOS GENES FILTRADOS (referencia) =====\n")
res_all <- fit_limma(log2_filtrado, design, contr)
print(deg_summary(res_all$tt))
write.csv(deg_summary(res_all$tt),
          file.path(OUT_DIR, "DEG_summary_ALL.csv"), row.names = FALSE)
run_plots(log2_filtrado, subtypes, res_all, "ALL")   # <-- AÑADIDO: gráficos ALL

# =============================================================================
# 8) Análisis principal — TOP genes más variables
# =============================================================================
cat("\n===== B) TOP", TOPN_VAR, "GENES MÁS VARIABLES =====\n")
vars     <- apply(log2_filtrado, 1, var, na.rm = TRUE)
top_ids  <- names(sort(vars, decreasing = TRUE))[seq_len(min(TOPN_VAR, length(vars)))]
expr_top <- log2_filtrado[top_ids, , drop = FALSE]

res_top  <- fit_limma(expr_top, design, contr)
LumA_top  <- res_top$tt[["LumA_vs_others"]]
LumB_top  <- res_top$tt[["LumB_vs_others"]]
Basal_top <- res_top$tt[["Basal_vs_others"]]

cat("Filas TOP (LumA/LumB/Basal):",
    nrow(LumA_top), nrow(LumB_top), nrow(Basal_top), "\n")
print(deg_summary(res_top$tt))
write.csv(deg_summary(res_top$tt),
          file.path(OUT_DIR, "DEG_summary_TOP.csv"), row.names = FALSE)
run_plots(expr_top, subtypes, res_top, paste0("TOP", TOPN_VAR))  # <-- AÑADIDO: gráficos TOP

# =============================================================================
# 9) Genes por dirección y unicidad por (subtipo, dirección)
#    ÚNICO-UP en X   = UP en X y NO UP en los otros dos
#    ÚNICO-DOWN en X = DOWN en X y NO DOWN en los otros dos
# =============================================================================
up   <- function(tt) rownames(tt[tt$adj.P.Val < PADJ_CUT & tt$logFC >  LFC_CUT, , drop = FALSE])
down <- function(tt) rownames(tt[tt$adj.P.Val < PADJ_CUT & tt$logFC < -LFC_CUT, , drop = FALSE])

up_LumA <- up(LumA_top);   down_LumA <- down(LumA_top)
up_LumB <- up(LumB_top);   down_LumB <- down(LumB_top)
up_Basal <- up(Basal_top); down_Basal <- down(Basal_top)

unique_UP_LumA    <- setdiff(up_LumA,    union(up_LumB,   up_Basal))
unique_DOWN_LumA  <- setdiff(down_LumA,  union(down_LumB, down_Basal))
unique_UP_LumB    <- setdiff(up_LumB,    union(up_LumA,   up_Basal))
unique_DOWN_LumB  <- setdiff(down_LumB,  union(down_LumA, down_Basal))
unique_UP_Basal   <- setdiff(up_Basal,   union(up_LumA,   up_LumB))
unique_DOWN_Basal <- setdiff(down_Basal, union(down_LumA, down_LumB))

unique_by_dir <- data.frame(
  Subtipo     = SUBTYPE_LEVELS,
  Unicos_UP   = c(length(unique_UP_LumA),   length(unique_UP_LumB),   length(unique_UP_Basal)),
  Unicos_DOWN = c(length(unique_DOWN_LumA), length(unique_DOWN_LumB), length(unique_DOWN_Basal))
)
cat("\n--- Genes únicos por (subtipo, dirección) ---\n"); print(unique_by_dir)

# =============================================================================
# 10) Exportar salidas (mismas que la versión original, en OUT_DIR)
# =============================================================================
# a) Listas simples (una columna de genes)
listas <- list(
  unique_UP_LumA = unique_UP_LumA,     unique_DOWN_LumA = unique_DOWN_LumA,
  unique_UP_LumB = unique_UP_LumB,     unique_DOWN_LumB = unique_DOWN_LumB,
  unique_UP_Basal = unique_UP_Basal,   unique_DOWN_Basal = unique_DOWN_Basal
)
for (nm in names(listas))
  write.table(listas[[nm]], file.path(OUT_DIR, paste0(nm, "_TOP.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# b) Panel de expresión (unión de todos los únicos) + metadatos de muestras
panel_genes <- unique(unlist(listas))
expr_panel  <- expr_top[intersect(panel_genes, rownames(expr_top)), , drop = FALSE]
write.table(cbind(Gene = rownames(expr_panel), expr_panel),
            file.path(OUT_DIR, "expr_panel_unique_TOP.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

sample_info <- data.frame(
  Sample  = colnames(expr_panel),
  Subtype = subtypes[match(colnames(expr_panel), colnames(log2_filtrado))],
  stringsAsFactors = FALSE
)
write.table(sample_info, file.path(OUT_DIR, "sample_info_TOP.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# c) Catálogo combinado con etiqueta (Gene | Subtype | Direction | logFC | FDR)
#    Estructura análoga al catálogo del one-vs-one -> entrada a la intersección.
add_tag <- function(genes, df, subtype, direction) {
  if (!length(genes)) return(NULL)
  sub <- df[genes, c("logFC", "adj.P.Val"), drop = FALSE]
  data.frame(Gene = rownames(sub), Subtype = subtype, Direction = direction,
             logFC = sub$logFC, FDR = sub$adj.P.Val, stringsAsFactors = FALSE)
}
catalogo <- do.call(rbind, list(
  add_tag(unique_UP_LumA,    LumA_top,  "LumA",  "UP"),
  add_tag(unique_DOWN_LumA,  LumA_top,  "LumA",  "DOWN"),
  add_tag(unique_UP_LumB,    LumB_top,  "LumB",  "UP"),
  add_tag(unique_DOWN_LumB,  LumB_top,  "LumB",  "DOWN"),
  add_tag(unique_UP_Basal,   Basal_top, "Basal", "UP"),
  add_tag(unique_DOWN_Basal, Basal_top, "Basal", "DOWN")
))
if (!is.null(catalogo))
  write.table(catalogo, file.path(OUT_DIR, "unique_genes_catalog_TOP.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nSalidas escritas en", OUT_DIR, ":\n",
    "- mapa_muestras.tsv\n",
    "- DEG_summary_ALL.csv / DEG_summary_TOP.csv\n",
    "- unique_{UP,DOWN}_{LumA,LumB,Basal}_TOP.tsv (listas simples)\n",
    "- expr_panel_unique_TOP.tsv + sample_info_TOP.tsv\n",
    "- unique_genes_catalog_TOP.tsv (entrada a la intersección)\n",
    "- plots_ALL.pdf / plots_TOP", TOPN_VAR, ".pdf (PCA, Volcano, MA, Heatmap)\n")
