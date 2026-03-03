########################################################################
##  Mixing Polygon Simulation (2 isotopes)                            ##
##  v1.3, March 2026, R version 4.4.1                                 ##
##  From Smith et al 2013, Methods in Ecology and Evolution           ##
##  https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12048
##  contact: james.smith@unsw.edu.au; james.a.smith@dpird.nsw.gov.au  ##
########################################################################

## This code wraps the point-in-polygon calculation, plotting, and saving in
## one function.


mix_poly_sim <- function(
    sources, mixture, TEF,
    its = 1500,
    min_C = -50, max_C = -20,
    min_N = -2,  max_N = 10,
    res   = 250,
    out_dir = NULL,
    labels = NULL,
    seed = 117
) {
  if (!is.null(out_dir)) dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  set.seed(seed)
  
  # grid
  C_g <- seq(min_C, max_C, length.out = res)
  N_g <- seq(min_N, max_N, length.out = res)
  grid <- expand.grid(C = C_g, N = N_g)
  
  n_src <- nrow(sources)
  n_mix <- nrow(mixture)
  
  # iteration storage
  Par_values <- matrix(0, nrow = its, ncol = (n_src*4 + 3))
  p          <- matrix(0L, nrow = its, ncol = n_mix)
  mix_reg    <- matrix(0L, nrow = res, ncol = res)
  
  # shoelace polygon area (for convex hull polygon)
  # does the same as splancs::areapl()
  poly_area_shoelace <- function(x, y) {
    x <- c(x, x[1]); y <- c(y, y[1])
    0.5 * abs(sum(x[-1]*y[-length(y)] - x[-length(x)]*y[-1]))
  }
  
  # progress bar
  pb <- progress::progress_bar$new(
    format = "  sim [:bar] :current/:total (:percent) eta::eta",
    total = its, clear = FALSE, width = 60 )
  
  for (i in 1:its) {
    # draw source isotopes and TEF enrichments (columns: meanC, sdC, meanN, sdN)
    vC <- rnorm(n_src, mean = sources[,1], sd = sources[,2])
    vN <- rnorm(n_src, mean = sources[,3], sd = sources[,4])
    fC <- rnorm(n_src, mean = TEF[,1],     sd = TEF[,2])
    fN <- rnorm(n_src, mean = TEF[,3],     sd = TEF[,4])
    
    VC <- vC + fC
    VN <- vN + fN
    
    hull <- chull(VC, VN)
    hull_a <- c(hull, hull[1])
    
    # consumers point-in-polygon
    P <- sp::point.in.polygon(mixture[,1], mixture[,2], VC[hull_a], VN[hull_a])
    p[i,] <- as.integer(P > 0)
    
    # polygon area + running variance
    poly_a <- poly_area_shoelace(VC[hull_a], VN[hull_a])
    
    # mixing region (grid) point-in-polygon
    m_r <- sp::point.in.polygon(grid$C, grid$N, VC[hull_a], VN[hull_a])
    m_r_s <- matrix(as.integer(m_r > 0), nrow = res, byrow = FALSE)  # rows=N, cols=C
    mix_reg <- mix_reg + m_r_s
    
    area_col <- n_src*4 + 1
    iter_col <- n_src*4 + 2
    var_col  <- n_src*4 + 3
    
    Par_values[i, 1:(iter_col)] <- c(vC, vN, fC, fN, poly_a, i)
    Par_values[i, var_col]      <- var(Par_values[1:i, area_col])
    
    pb$tick()
  }
  
  # outputs
  probabilities <- colMeans(p)
  mix_reg <- mix_reg/its
  mix_reg[mix_reg == 0] <- NA
  
  # long df for ggplot
  mix_df <- expand.grid(d13C = C_g, d15N = N_g)
  mix_df$prob <- as.vector(mix_reg)
  
  # enriched source means (for plotting)
  sources_TEF <- data.frame(
    d13C = sources[,1] + TEF[,1],
    d15N = sources[,3] + TEF[,3] )
  
  # a theme with a box around the panel
  theme_box <- theme_minimal() +
    theme(panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.6))
  
  # FIG 1: running variance of polygon area
  var_df <- data.frame(
    iteration = Par_values[, ncol(Par_values)-1],
    variance  = Par_values[, ncol(Par_values)] )
  var_df <- var_df[!is.na(var_df$variance), ]  #first iteration has no variance
  
  area_plot <- ggplot(var_df, aes(iteration, variance)) +
    geom_line() +
    labs(x = "Iteration", y = "Running variance of polygon area") +
    theme_box
  
  # FIG 2: mixing region + contours + points
  cont_levels <- c(0.05, seq(0.1, 1, by = 0.1))
  
  mix_plot <- ggplot(mix_df, aes(d13C, d15N)) +
    geom_raster(aes(fill = prob), na.rm = TRUE) +
    geom_contour(aes(z = prob), breaks = cont_levels, linewidth = 0.6, colour = "black", na.rm = TRUE) +
    geom_point(data = sources_TEF, aes(d13C, d15N),
               inherit.aes = FALSE, shape = 4, size = 3, stroke = 1.2, colour = "white") +
    geom_point(data = mixture, aes(x = mixture[,1], y = mixture[,2]),
               inherit.aes = FALSE, size = 2) +
    scale_fill_gradientn(colours = c("blue","lightblue","green","lightgreen","yellow","red")) +
    labs(x = "d13C", y = "d15N", fill = "P(inside)") +
    coord_cartesian(expand = FALSE) +
    theme_box
  
  # FIG 3: bi-plot with single 95% contour (0.05) + SD error bars + labels
  # SDs displayed for enriched means assuming independence of SDs (visual aid)
  sources_TEF$sdC <- sqrt(sources[,2]^2 + TEF[,2]^2)
  sources_TEF$sdN <- sqrt(sources[,4]^2 + TEF[,4]^2)
  
  if (is.null(labels)) labels <- paste0("s", seq_len(n_src))
  sources_TEF$label <- labels
  
  biplot_95 <- ggplot(mix_df, aes(d13C, d15N)) +
    geom_contour(aes(z = prob), breaks = 0.05, linewidth = 0.7, colour = "black", na.rm = TRUE) +
    geom_point(data = sources_TEF, aes(d13C, d15N), inherit.aes = FALSE, shape = 15, size = 3) +
    geom_errorbarh(data = sources_TEF,
                   aes(y = d15N, xmin = d13C - sdC, xmax = d13C + sdC),
                   inherit.aes = FALSE, height = 0) +
    geom_errorbar(data = sources_TEF,
                  aes(x = d13C, ymin = d15N - sdN, ymax = d15N + sdN),
                  inherit.aes = FALSE, width = 0) +
    geom_text(data = sources_TEF, aes(d13C, d15N, label = label),
              inherit.aes = FALSE, vjust = -1) +
    geom_point(data = mixture, aes(x = mixture[,1], y = mixture[,2]),
               inherit.aes = FALSE, shape = 1, size = 2) +
    labs(x = "d13C", y = "d15N") +
    coord_cartesian(
      xlim = c(min_C, max_C),
      ylim = c(min_N, max_N),
      expand = FALSE ) +
    theme_box
  
  # write files (PDF figures only)
  if (!is.null(out_dir)) {
    # data
    write.csv(rbind(p, probabilities), file = file.path(out_dir, "Consumer_Probabilities.csv"), row.names = FALSE)
    
    col_names <- c(rep("d13C", n_src), rep("d15N", n_src),
                   rep("13C_TEF", n_src), rep("15N_TEF", n_src),
                   "Poly_Area","Iteration","Variance")
    col_nums  <- c(rep(1:n_src, 4), 0, 0, 0)
    colnames(Par_values) <- paste(col_names, col_nums)
    write.csv(Par_values, file = file.path(out_dir, "Parameter_Values.csv"), row.names = FALSE)
    
    # figs
    ggsave(file.path(out_dir, "Fig1_area_variance.pdf"), area_plot, width = 7, height = 4)
    ggsave(file.path(out_dir, "Fig3_mixing_region.pdf"), mix_plot,  width = 7, height = 5)
    ggsave(file.path(out_dir, "Fig5_biplot_95.pdf"),     biplot_95, width = 7, height = 5)
  }
  
  list(
    probabilities = probabilities,
    area_plot = area_plot,
    mix_plot = mix_plot,
    biplot_95 = biplot_95,
    Par_values = Par_values,
    pip = p,
    mix_reg = mix_reg,
    grid = list(C_g = C_g, N_g = N_g) )
}
