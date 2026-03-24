
#---------#    AUXILIARY FUNCTIONS     #---------#
#---------#   don't run this script    #---------#



# auxiliary function that finds the mode of a discrete vector
find_mode <- function(x) {
  freq_table <- table(x)
  mode <- as.numeric(names(freq_table[freq_table == max(freq_table)]))
  return(mode)
}



# ==============================================================================
# Function that plots the selected neurons locations in the hippocampus
# ==============================================================================

plot_neuron_locations <- function(neu1, neu2, shapes, cols = c("royalblue4", "gold")) {
  

  if (!exists("loc_neurons")) {
    stop("Global object 'loc_neurons' is missing.")
  }
  
  dataf <- data.frame(x = loc_neurons[,1], y = loc_neurons[,2])
  dataf$colorn <- "no"
  
  # Handle the case where the same neuron is passed for both arguments
  if (neu1 == neu2) {
    dataf$colorn[neu1] <- "both"
    title_str <- paste("Neuron", neu1, "(selected twice)")
    
    fill_vals <- c("gray22", "purple") 
    alpha_vals <- c(0.4, 1.0)
    shape_vals <- c(21, 23) # Circle, Diamond
    size_vals <- c(3, 5)
    
  } else {
    # Normal case: two different neurons
    dataf$colorn[neu1] <- "yes_neu1"
    dataf$colorn[neu2] <- "yes_neu2"
    title_str <- paste("Neurons", neu1, "and", neu2)
    
    # Define scales for distinct neurons
    # "no", "yes_neu1", "yes_neu2"
    fill_vals <- c("gray22", cols)
    alpha_vals <- c(0.4, 0.8, 0.8)
    shape_vals <- c(21, shapes) # Circle, Triangle up, Diamond
    size_vals <- c(3, 4, 4)
  }
  
  level_order <- unique(c("no", dataf$colorn[neu1], dataf$colorn[neu2]))
  level_order <- c("no", setdiff(level_order, "no"))
  dataf$colorn <- factor(dataf$colorn, levels = level_order)


  G13 <- ggplot(dataf) +
    scale_fill_manual(values = fill_vals) +
    scale_alpha_manual(values = alpha_vals) +
    scale_shape_manual(values = shape_vals) +
    scale_size_manual(values = size_vals) +
    
    geom_point(aes(x = x, y = y, 
                   fill = colorn, 
                   pch = colorn, 
                   alpha = colorn, 
                   size = colorn), 
               color = "black") +
    
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      aspect.ratio = 1,
      panel.background = element_rect(fill = 'transparent', color = NA),
      plot.background = element_rect(fill = 'transparent', color = NA),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "darkgray ", fill = NA),
      axis.line.x.bottom = element_line(color = "gray"),
      legend.position = "none", 
      legend.direction = 'vertical',
      legend.background = element_rect(fill = 'transparent', color = 'white'),
      legend.box.background = element_rect(fill = 'transparent'),
      legend.text = element_text(size = 10),
      strip.text = element_text(size = 10),
      text = element_text(size = 13),
      strip.background = element_rect(fill = NA, color = "gray")
    ) +
    ggtitle(title_str) +
    xlab("x coordinate \n (hippocampus)") +
    ylab("y coordinate \n (hippocampus)")
  
  return(G13)
}


# ==============================================================================
# Function that plots Smoothed Spatial Activity for a Single Neuron
# ==============================================================================

# neu: Integer index of the neuron to plot.
# data_df: Processed data frame containing position (pos1, pos2) and window indices (win). Defaults to global 'data_processed'.
# WIND_list: List of activity matrices per window. Defaults to global 'WIND'.
# idx_vec: Vector of window indices to include. Defaults to global 'idx'.

plot_neuron_smooth_spatial <- function(neu, 
                                       data_df = data, 
                                       WIND_list = WIND, 
                                       idx_vec = idx,
                                       npix = 200,
                                       upper = 1, smooth=1) {
  
  # Basic validation
  if (missing(neu)) stop("Neuron index 'neu' must be provided.")
  
  D_list <- lapply(idx_vec, function(window_idx) {
    # Find rows in data corresponding to the current window
    win_rows <- which(data_df$win == window_idx)
    
    # Extract activity for the specific neuron in this window
    act_vals <- WIND_list[[which(idx_vec == window_idx)]][neu, ]
    n_pos <- length(win_rows)
    n_act <- length(act_vals)
    if(n_pos != n_act){
      warning(paste("Window", window_idx, ": Mismatch between position rows (", n_pos, 
                    ") and activity values (", n_act, "). Truncating to shorter length."))
      min_len <- min(n_pos, n_act)
      win_rows <- win_rows[1:min_len]
      act_vals <- act_vals[1:min_len]
    }
    
    data.frame(
      win = window_idx,
      x = data_df$pos1[win_rows],
      y = data_df$pos2[win_rows],
      act = act_vals
    )
  })
  
  D <- dplyr::bind_rows(D_list)

  data_for_mba <- cbind(x = D$x, y = D$y, z = D$act)

  set.seed(123) 
  fit <- MBA::mba.surf(
    data_for_mba,
    no.X = npix,
    no.Y = npix,
    extend = FALSE
  )
  
  # 4. Smoothing and formatting results

  fit$xyz.est$z2 <- fields::image.smooth(fit$xyz.est$z, theta = smooth)$z
  fit$xyz.est$z2[is.na(fit$xyz.est$z)] <- 0
  fit$xyz.est$z2[(fit$xyz.est$z2)<0] <- 0
  fit$xyz.est$z2[(fit$xyz.est$z2)>upper] <- upper
  
  # Create grid for plotting
  new_grid <- tidyr::expand_grid(
    x_coord = seq(from = fit$xyz.est$x[1], to = fit$xyz.est$x[length(fit$xyz.est$x)], length.out = npix),
    y_coord = seq(from = fit$xyz.est$y[1], to = fit$xyz.est$y[length(fit$xyz.est$y)], length.out = npix)
  )
  
  plot_data <- data.frame(
    x_coord = new_grid$x_coord,
    y_coord = new_grid$y_coord,
    # Flatten the matrix to a vector column-wise to match the grid
    fill_val = as.vector(fit$xyz.est$z2) 
  )
  
  V <- ggplot() +
    geom_path(aes(x = data[idcirc,]$pos1, y= data[idcirc,]$pos2), alpha=.5, lwd=0.25)+
    geom_tile(data = plot_data, aes(x = x_coord, y = y_coord, fill = fill_val),alpha=.85) +
    scale_fill_viridis_c("Spike probability   ",
                         breaks = c(0.0, upper/2, upper), 
                         labels = c("0.0", upper/2, upper),
                         option = "magma",
                         limits = c(0,upper),
                         direction = -1) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      aspect.ratio = 1,
      panel.background = element_rect(fill='transparent', color=NA),
      plot.background = element_rect(fill='transparent', color=NA),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "darkgray ", fill=NA),
      axis.line.x.bottom = element_line(color="gray"),
      legend.position = "bottom",
      legend.title = element_text(vjust = 0.5),
      legend.background = element_rect(fill='transparent', color = 'transparent'),
      legend.box.background = element_rect(fill='transparent', color = 'transparent'),
      legend.text = element_text(size=10),
      strip.background = element_rect( fill=NA, color="gray76" ),
      text = element_text(size=13)
    ) +
    xlab("x coordinate \n (environment)") +
    ylab("y coordinate \n (environment)") +
    ggtitle(paste("Neuron", neu)) +
    guides(fill = guide_colourbar(barwidth = 8))
  return(V)
}

