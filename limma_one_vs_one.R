# =============================================================================
# Análisis diferencial de subtipos BRCA — LIMMA sin voom, contrastes POR PARES
# (one-vs-one: LumA vs LumB, LumA vs Basal, LumB vs Basal)
#
# Entrada : matriz FPKM gen x muestra, con nombres de columna tipo BRCA_LumA.1
# Salida  : listas de genes consistentes y ÚNICOS por subtipo y dirección
#           (UP/DOWN), que alimentan el análisis de intersección posterior.
#
# NOTA METODOLÓGICA: la DE se corre dos veces — (A) sobre todos los genes
# filtrados por QC y (B) sobre los TOP 1000 más variables. Las listas finales
# usan (B). Seleccionar por varianza y luego testear diferencias entre grupos
# sobre los mismos datos puede sesgar hacia genes que separan subtipos; esta es
# una decisión consciente y debe documentarse en la tesis. Para comparar, basta
# cambiar TT_USE <- tt_all / EXPR_USE <- expr_all más abajo.
# =============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(EnhancedVolcano)
})

# ------------------------------- Parámetros ---------------------------------
dat   <- read.csv("C:/Users/danil/OneDrive/Escritorio/Propuesta de Tesis/Analisis python/552_expresion_genes_filtrado.tsv", sep="")  # 
OUT_DIR      <- "resultados_552_one_vs_one"

PADJ_CUT     <- 0.05     # umbral de FDR (BH) para significancia
LFC_CUT      <- 1        # umbral de |log2FC|
ZERO_FRAC    <- 0.50     # se descarta un gen si tiene > 50% de ceros
COR_CUT      <- 0.75     # umbral de correlación media para señalar outliers
TOPN_VAR     <- 1000     # nº de genes más variables para el análisis (B)
TOPN_HEATMAP <- 100      # genes por contraste en los heatmaps

# Muestras a eliminar tras la inspección de outliers.
# manualmente). Deja el vector vacío c() si no quieres eliminar ninguna.
OUTLIERS_TO_REMOVE <- c()   

SUBTYPE_LEVELS <- c("LumA", "LumB", "Basal")  # orden de referencia
SAVE_PLOTS   <- TRUE     # TRUE -> guarda PDFs en OUT_DIR; FALSE -> plot a pantalla

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1) Carga y conversión a numérico
# =============================================================================
cat("Matriz cruda:", nrow(dat), "genes x", ncol(dat) - 1, "muestras\n")

genes_all <- as.character(dat[[1]])
raw       <- dat[, -1, drop = FALSE]

# --- Separar filas de metadatos embebidas (no numéricas) al inicio del archivo ---
# El TSV trae filas de metadatos ANTES de la expresión: etiqueta de muestra y
# barcode TCGA del paciente. Se detectan como filas iniciales no numéricas, de
# modo que si cambia su número el script se adapta solo.
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

# make.unique: genes duplicados (pseudoautosómicos X/Y: CD99, SLC25A6, CRLF2...)
# reciben sufijo .1 para poder usarse como rownames. MISMO criterio en OvR.
rownames(mat) <- make.unique(genes)
mat[] <- lapply(mat, function(x) suppressWarnings(as.numeric(x)))

# =============================================================================
# 2) Filtrado por QC y transformación log2
# =============================================================================
# 2a. Quitar genes todo-cero
mat <- mat[rowSums(mat, na.rm = TRUE) > 0, , drop = FALSE]

# Registro de muestras originales (trazabilidad)
samples_original <- colnames(mat)

# 2b. log2(FPKM + 1)
log2_mat <- log2(mat + 1)

# 2c. Filtrar genes con demasiados ceros (> ZERO_FRAC de las muestras)
ceros  <- rowSums(mat == 0, na.rm = TRUE)
umbral <- ncol(mat) * ZERO_FRAC
log2_filtrado <- log2_mat[ceros < umbral, , drop = FALSE]

# 2d. Quitar genes con varianza cero
keep <- apply(log2_filtrado, 1, var, na.rm = TRUE) > 0
log2_filtrado <- log2_filtrado[keep, , drop = FALSE]
cat("Tras filtrado QC:", nrow(log2_filtrado), "genes x",
    ncol(log2_filtrado), "muestras\n")

# =============================================================================
# 3) Detección y eliminación de outliers de muestra
# =============================================================================
log2_cc  <- log2_filtrado[complete.cases(log2_filtrado), , drop = FALSE]
cor_prom <- apply(cor(log2_cc), 1, mean)

cat("\nCorrelación media por muestra (resumen):\n")
print(summary(cor_prom))

posibles_out <- names(which(cor_prom < COR_CUT))
if (length(posibles_out)) {
  cat("Posibles outliers (corr media <", COR_CUT, "):\n")
  print(posibles_out)
} else {
  cat("No se detectaron outliers bajo el umbral", COR_CUT, "\n")
}

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
# Nombres esperados: BRCA_LumA.1, BRCA_LumB.10, ... -> se toma lo previo al
# primer punto y se retira el prefijo "BRCA_".
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
# 5) Diseño y contrastes por pares
# =============================================================================
design <- model.matrix(~ 0 + subtypes)
colnames(design) <- levels(subtypes)   # "LumA","LumB","Basal"

contr <- makeContrasts(
  LumA_vs_LumB  = LumA - LumB,
  LumA_vs_Basal = LumA - Basal,
  LumB_vs_Basal = LumB - Basal,
  levels = design
)

# =============================================================================
# 6) Funciones auxiliares
# =============================================================================

# Ajusta limma y devuelve las TopTables (FDR/BH) por contraste
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

# Conteo de DEGs por contraste (usa FDR, coherente con el resto del script)
deg_summary <- function(tt_list, padj_cut = PADJ_CUT, lfc_cut = LFC_CUT) {
  data.frame(
    Contraste       = names(tt_list),
    FDR             = sapply(tt_list, function(tt)
                        sum(tt$adj.P.Val < padj_cut, na.rm = TRUE)),
    FDR_y_logFC     = sapply(tt_list, function(tt)
                        sum(tt$adj.P.Val < padj_cut &
                            abs(tt$logFC) > lfc_cut, na.rm = TRUE)),
    row.names = NULL
  )
}

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
# 7) Análisis A — todos los genes filtrados por QC
# =============================================================================
cat("\n===== A) TODOS LOS GENES FILTRADOS =====\n")
expr_all <- log2_filtrado
res_all  <- fit_limma(expr_all, design, contr)

sum_all <- deg_summary(res_all$tt)
print(sum_all)
write.csv(sum_all, file.path(OUT_DIR, "DEG_summary_ALL.csv"), row.names = FALSE)
run_plots(expr_all, subtypes, res_all, "ALL")

# =============================================================================
# 8) Análisis B — TOP genes más variables
# =============================================================================
cat("\n===== B) TOP", TOPN_VAR, "GENES MÁS VARIABLES =====\n")
v_all    <- apply(log2_filtrado, 1, var, na.rm = TRUE)
k        <- min(TOPN_VAR, length(v_all))
expr_top <- log2_filtrado[order(v_all, decreasing = TRUE)[seq_len(k)], , drop = FALSE]
cat("Dim TOP:", nrow(expr_top), "x", ncol(expr_top), "\n")

res_top <- fit_limma(expr_top, design, contr)
sum_top <- deg_summary(res_top$tt)
print(sum_top)
write.csv(sum_top, file.path(OUT_DIR, "DEG_summary_TOP.csv"), row.names = FALSE)
run_plots(expr_top, subtypes, res_top, paste0("TOP", TOPN_VAR))

# =============================================================================
# 9) Genes consistentes y únicos por subtipo y dirección
#    (esta es la salida que alimenta el análisis de intersección)
# =============================================================================
# Selecciona qué análisis usar para las listas finales (ver nota metodológica)
TT_USE   <- res_top$tt
EXPR_USE <- expr_top

tt_AvB <- TT_USE[["LumA_vs_LumB"]]    # LumA - LumB
tt_AvC <- TT_USE[["LumA_vs_Basal"]]   # LumA - Basal
tt_BvC <- TT_USE[["LumB_vs_Basal"]]   # LumB - Basal

sig_up   <- function(tt) rownames(tt[tt$adj.P.Val < PADJ_CUT & tt$logFC >  LFC_CUT, , drop = FALSE])
sig_down <- function(tt) rownames(tt[tt$adj.P.Val < PADJ_CUT & tt$logFC < -LFC_CUT, , drop = FALSE])

# Consistencia por subtipo: mismo signo frente a los DOS subtipos restantes
Up_LumA    <- intersect(sig_up(tt_AvB),   sig_up(tt_AvC))     # A > B  y  A > C
Down_LumA  <- intersect(sig_down(tt_AvB), sig_down(tt_AvC))   # A < B  y  A < C
Up_LumB    <- intersect(sig_down(tt_AvB), sig_up(tt_BvC))     # B > A  y  B > C
Down_LumB  <- intersect(sig_up(tt_AvB),   sig_down(tt_BvC))   # B < A  y  B < C
Up_Basal   <- intersect(sig_down(tt_AvC), sig_down(tt_BvC))   # C > A  y  C > B
Down_Basal <- intersect(sig_up(tt_AvC),   sig_up(tt_BvC))     # C < A  y  C < B

# Unicidad POR DIRECCIÓN (UP contra UP, DOWN contra DOWN)
unique_UP_LumA    <- setdiff(Up_LumA,    union(Up_LumB,   Up_Basal))
unique_DOWN_LumA  <- setdiff(Down_LumA,  union(Down_LumB, Down_Basal))
unique_UP_LumB    <- setdiff(Up_LumB,    union(Up_LumA,   Up_Basal))
unique_DOWN_LumB  <- setdiff(Down_LumB,  union(Down_LumA, Down_Basal))
unique_UP_Basal   <- setdiff(Up_Basal,   union(Up_LumA,   Up_LumB))
unique_DOWN_Basal <- setdiff(Down_Basal, union(Down_LumA, Down_LumB))

cat("\n--- Genes únicos por dirección ---\n")
cat("UP   -> LumA:", length(unique_UP_LumA),
    " LumB:", length(unique_UP_LumB),
    " Basal:", length(unique_UP_Basal), "\n")
cat("DOWN -> LumA:", length(unique_DOWN_LumA),
    " LumB:", length(unique_DOWN_LumB),
    " Basal:", length(unique_DOWN_Basal), "\n")

# Tabla larga: Gene | Subtype | Direction
unique_genes_list_dir <- list(
  LumA_UP = unique_UP_LumA,   LumA_DOWN = unique_DOWN_LumA,
  LumB_UP = unique_UP_LumB,   LumB_DOWN = unique_DOWN_LumB,
  Basal_UP = unique_UP_Basal, Basal_DOWN = unique_DOWN_Basal
)
unique_genes_df_dir <- do.call(rbind, lapply(names(unique_genes_list_dir), function(name) {
  g <- unique_genes_list_dir[[name]]
  if (!length(g)) return(NULL)
  data.frame(Gene = g,
             Subtype   = sub("_.*", "", name),
             Direction = sub(".*_", "", name),
             stringsAsFactors = FALSE)
}))
cat("\nTotal genes únicos con dirección:", nrow(unique_genes_df_dir), "\n")

# Anexar la matriz de expresión y guardar (entrada a la intersección)
expr_unique <- EXPR_USE[match(unique_genes_df_dir$Gene, rownames(EXPR_USE)), , drop = FALSE]
expr_unique_full_dir <- cbind(unique_genes_df_dir, expr_unique)

write.table(unique_genes_df_dir,
            file.path(OUT_DIR, "genes_unicos_OvO.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(expr_unique_full_dir,
            file.path(OUT_DIR, "genes_unicos_OvO_con_expresion.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nSalidas escritas en", OUT_DIR, ":\n",
    "- mapa_muestras.tsv\n",
    "- DEG_summary_ALL.csv / DEG_summary_TOP.csv\n",
    "- genes_unicos_OvO.tsv (lista para intersección)\n",
    "- genes_unicos_OvO_con_expresion.tsv\n")
