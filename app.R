library(shiny)
library(minpack.lm)

ui <- fluidPage(
  
  # Custom CSS for a slightly cleaner, flat look (optional but recommended)
  tags$head(
    tags$style(HTML("
      /* Use a slightly different default font for readability */
      body {
        font-family: 'Open Sans', sans-serif;
      }
      /* Style for the main action button */
      .btn-primary {
        background-color: #2c3e50;
        border-color: #2c3e50;
        color: white;
      }
      .btn-primary:hover {
        background-color: #1a242f;
        border-color: #1a242f;
      }
      /* Style for the success button (Show Data) */
      .btn-success {
        background-color: #18bc9c;
        border-color: #18bc9c;
      }
      .btn-success:hover {
        background-color: #14a287;
        border-color: #14a287;
      }
      /* Consistent look for well panels */
      .well {
        border-left: 5px solid #2c3e50;
        background-color: #f8f8f8;
        border-radius: 4px;
        box-shadow: none;
      }
      /* Make tabs slightly more prominent */
      .nav-pills > li.active > a, .nav-pills > li.active > a:hover, .nav-pills > li.active > a:focus {
        background-color: #18bc9c;
      }
    "))
  ),
  
  titlePanel(
    div(
      span("Glycylglycine (GG) Kinetics Model", style = "font-weight: bold; color: #2c3e50;"),
      p("Hydrolysis & Cyclization", style = "font-size: 16px; margin-top: 5px; color = #7f8c8d;")
    ),
    windowTitle = "GG Kinetics Model"
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 4,
      h4("1. Enter Your Data", style = "margin-bottom: 10px;"),
      helpText("Paste your experimental data here. Ensure columns are separated by tabs or commas."),
      tags$textarea(
        id = "data_input",
        rows = 10,
        cols = 50,
        placeholder = "time\tGG\tcGG\tG\n0\t1.8428\t0\t0\n2\t1.6288\t0.0708\t0.4237\n..."
      ),
      
      hr(),
      h4("2. Model Options", style = "margin-top: 20px;"),
      
      checkboxInput(
        inputId = "use_weights",
        label = "Use robust weighting",
        value = TRUE
      ),
      helpText("This weights the fit to minimize the relative error, useful for data with varying magnitude."),
      
      hr(),
      div(
        actionButton("fit_btn", "Run Model & Analyze", class = "btn-primary", icon = icon("dna")),
        # Button to show data in a dedicated tab
        actionButton("show_data_btn", "Show Export Data (Copy Manually)", class = "btn-success", icon = icon("table"))
      )
    ),
    
    mainPanel(
      width = 8,
      tabsetPanel(
        id = "main_tabs",
        type = "pills",
        
        tabPanel("Kinetic Fit",
                 value = "fit_tab",
                 h4("Kinetic Fit Overview"),
                 br(),
                 fluidRow(
                   column(6,
                          wellPanel(
                            h5(strong("Fitted Parameters"), style = "margin-top: 0;"),
                            verbatimTextOutput("fit_results")
                          )
                   ),
                   column(6,
                          wellPanel(
                            h5(strong("R-squared & Baselines"), style = "margin-top: 0;"),
                            verbatimTextOutput("metrics_results")
                          )
                   )
                 ),
                 br(),
                 plotOutput("fit_plot", height = "500px"),
                 hr(),
                 h4("Residuals by Species"),
                 p("A good fit should show residuals scattered randomly around zero."),
                 plotOutput("resid_plot", height = "400px")
        ),
        
        tabPanel("Glycine Fit",
                 h4("Alternative Glycine-Only Fit"),
                 p("Compares the full kinetic model (blue line) with a simplified, single-species fit (orange line) for the glycine data."),
                 br(),
                 fluidRow(
                   column(6,
                          wellPanel(
                            h5(strong("Linear Fit Results"), style = "margin-top: 0;"),
                            verbatimTextOutput("g_fit_results")
                          )
                   ),
                   column(6,
                          wellPanel(
                            h5(strong("Model-Based Fit"), style = "margin-top: 0;"),
                            verbatimTextOutput("model_g_results")
                          )
                   )
                 ),
                 plotOutput("g_fit_plot", height = "500px")
        ),
        
        # Tab for manual data export
        tabPanel("Export Data",
                 value = "export_data_tab",
                 h4("Data for Copy & Paste (Tab Separated)"),
                 helpText(strong("Instructions:"), "This file contains three sections: **1) Parameters**, **2) Raw Data + Fitted Values at Data Points**, and **3) High-Resolution Curves for Plotting** (Glycine only). Select all text below and manually copy (Ctrl+C / Cmd+C). Then paste (Ctrl+V / Cmd+V) directly into Excel, Google Sheets, or a CSV file editor."),
                 verbatimTextOutput("raw_export_data")
        )
      )
    )
  )
)

server <- function(input, output, session) { 
  
  # Reactive values to store all calculated data
  results <- reactiveValues(
    fit = NULL,
    pars_est = NULL,
    R2_global = NULL,
    R2_by_species = NULL,
    long_data = NULL,
    exp_data = NULL,
    g_fit_model = NULL,
    model_g_pred = NULL,
    data_export_text = NULL,
    params_export_text = NULL,
    glycine_plot_data_text = NULL # NEW: For high-res plotting data
  )
  
  observeEvent(input$fit_btn, {
    req(input$data_input)
    
    # --- Data parsing & checks ---
    data_lines <- strsplit(input$data_input, "\n")[[1]]
    data_lines <- trimws(data_lines)
    data_lines <- data_lines[data_lines != ""]
    data_split <- lapply(data_lines, function(x) strsplit(x, "\t")[[1]])
    if (length(data_split) < 2) {
      showNotification("Please include at least a header and one row of data.", type = "error")
      return()
    }
    col_names <- data_split[[1]]
    data_matrix <- data_split[-1]
    data_matrix <- data_matrix[sapply(data_matrix, length) == length(col_names)]
    if (length(data_matrix) == 0) {
      showNotification("No valid data rows found. Check formatting.", type = "error")
      return()
    }
    df <- as.data.frame(do.call(rbind, data_matrix), stringsAsFactors = FALSE)
    colnames(df) <- col_names
    df[] <- lapply(df, function(col) as.numeric(gsub(",", ".", col)))
    if (anyNA(df)) {
      showNotification("Some values could not be converted to numbers. Check input.", type = "error")
      return()
    }
    
    exp_data <- df[order(df$time), ] 
    n <- nrow(exp_data)
    results$exp_data <- exp_data
    
    long_data <- data.frame(
      time = rep(exp_data$time, 3),
      conc = c(exp_data$GG, exp_data$cGG, exp_data$G),
      species = factor(rep(c("GG", "cGG", "G"), each = n)) 
    )
    results$long_data <- long_data
    
    use_weights <- isTRUE(input$use_weights)
    small <- 1e-6
    obs_w <- if (use_weights) 1 / sqrt(pmax(long_data$conc, small)) else rep(1, nrow(long_data))
    
    # --- Model function & fitting (kinetic) ---
    model_fun_pars <- function(pars, time, species) {
      GG0 <- exp(pars["logGG0"])
      k1 <- exp(pars["logk1"])
      k2 <- exp(pars["logk2"])
      ksum <- k1 + k2
      GG_pred <- GG0 * exp(-ksum * time)
      cGG_pred <- GG0 * k2 / ksum * (1 - exp(-ksum * time))
      G_pred <- 2 * GG0 * k1 / ksum * (1 - exp(-ksum * time))
      cGG_pred <- cGG_pred + pars["bg_cGG"]
      G_pred <- G_pred + pars["bg_G"]
      ifelse(species == "GG", GG_pred, ifelse(species == "cGG", cGG_pred, G_pred))
    }
    
    resid_fun <- function(pars, time, species, conc, w) {
      pred <- model_fun_pars(pars, time, species)
      (pred - conc) * w
    }
    
    GG0_start <- max(exp_data$GG, na.rm = TRUE)
    t_span <- max(exp_data$time) - min(exp_data$time)
    k_guess <- if (t_span > 0) log(2) / max(1, t_span) else 0.1
    
    start_pars <- c(
      logGG0 = log(ifelse(GG0_start > 0, GG0_start, 0.5)),
      logk1 = log(max(1e-3, k_guess * 0.5)),
      logk2 = log(max(1e-3, k_guess * 0.2)),
      bg_cGG = min(exp_data$cGG, na.rm = TRUE) * 0.5,
      bg_G = min(exp_data$G, na.rm = TRUE) * 0.5
    )
    start_vec <- as.numeric(start_pars)
    names(start_vec) <- names(start_pars)
    
    fit <- tryCatch(
      nls.lm(
        par = start_vec, fn = resid_fun, time = long_data$time, species = long_data$species,
        conc = long_data$conc, w = obs_w, control = nls.lm.control(maxiter = 200, ftol = 1e-10)
      ),
      error = function(e) e
    )
    
    if (inherits(fit, "error")) {
      showNotification(paste("Fit failed:", fit$message), type = "error", duration = 8)
      results$fit <- "fail"
      return()
    }
    
    results$fit <- fit
    results$pars_est <- fit$par
    
    long_data$predicted <- model_fun_pars(results$pars_est, long_data$time, long_data$species)
    long_data$resid <- long_data$conc - long_data$predicted
    results$long_data <- long_data
    
    ss_res_global <- sum((long_data$conc - long_data$predicted)^2)
    ss_tot_global <- sum((long_data$conc - mean(long_data$conc))^2)
    results$R2_global <- 1 - ss_res_global / ss_tot_global
    
    R2_by_species_list <- by(long_data, long_data$species, function(subset_df) {
      ss_res <- sum((subset_df$conc - subset_df$predicted)^2)
      ss_tot <- sum((subset_df$conc - mean(subset_df$conc))^2)
      R2 <- ifelse(ss_tot > 0, 1 - ss_res / ss_tot, NA_real_)
      data.frame(species = subset_df$species[1], R2 = R2)
    })
    results$R2_by_species <- do.call(rbind, R2_by_species_list)
    
    # --- Glycine Fit Calculations & R2s ---
    g_data <- exp_data[, c("time", "G")]
    g_fit_model_success <- FALSE
    
    # 1. Linear Glycine Fit
    tryCatch({
      t_half_guess <- g_data$time[which.min(abs(g_data$G - max(g_data$G) / 2))]
      k_app_guess <- if(length(t_half_guess) > 0 && is.finite(t_half_guess) && t_half_guess > 0) log(2) / t_half_guess else 0.05
      g_fit_model <- nlsLM(
        G ~ G_inf * (1 - exp(-k_app * time)),
        data = g_data,
        start = list(G_inf = max(g_data$G), k_app = k_app_guess),
        lower = c(0, 0),
        upper = c(Inf, Inf)
      )
      results$g_fit_model <- g_fit_model
      g_fit_model_success <- TRUE
      
      # Calculate R2 for Linear Glycine Fit
      g_pred <- predict(results$g_fit_model, g_data)
      ss_res_g_lin <- sum((g_data$G - g_pred)^2)
      ss_tot_g <- sum((g_data$G - mean(g_data$G))^2) # Need ss_tot_g for R2_model_g later
      R2_g_lin <- 1 - ss_res_g_lin / ss_tot_g
      
    }, error = function(e) {
      results$g_fit_model <- "fail"
      ss_tot_g <- sum((g_data$G - mean(g_data$G))^2)
      R2_g_lin <- NA_real_
    })
    
    # 2. Model-based G predictions (at experimental points)
    pars_est <- results$pars_est
    GG0_fit <- exp(pars_est["logGG0"])
    k1_fit <- exp(pars_est["logk1"])
    k2_fit <- exp(pars_est["logk2"])
    ksum_fit <- k1_fit + k2_fit
    bg_G_fit <- pars_est["bg_G"]
    
    G_model_pred_exp <- 2 * GG0_fit * k1_fit / ksum_fit * (1 - exp(-ksum_fit * g_data$time)) + pars_est["bg_G"]
    ss_res_model <- sum((g_data$G - G_model_pred_exp)^2)
    R2_model_g <- 1 - ss_res_model / ss_tot_g
    
    
    # --- Data Export Preparation (3 Sections) ---
    
    # --- Section 1: Parameters (Kinetic + Glycine) ---
    param_names <- c(
      "Initial_GG_GG0_mM", "k1_hydrolysis_h-1", "k2_cyclization_h-1", "k_sum_h-1", 
      "t_half_total_h", "t_half_k1_h", "t_half_k2_h", 
      "Global_R2", "R2_GG", "R2_cGG", "R2_G", 
      "Baseline_cGG_mM", "Baseline_G_mM",
      "R2_Glycine_Model_Fit"
    )
    param_values <- c(
      GG0_fit, k1_fit, k2_fit, ksum_fit, 
      log(2) / ksum_fit, log(2) / k1_fit, log(2) / k2_fit, 
      results$R2_global, results$R2_by_species$R2[1], results$R2_by_species$R2[2], results$R2_by_species$R2[3], 
      pars_est["bg_cGG"], pars_est["bg_G"],
      R2_model_g
    )
    
    if (g_fit_model_success) {
      coefs <- coef(results$g_fit_model)
      param_names <- c(param_names, 
                       "Glycine_Linear_Fit_k_app_h-1", "Glycine_Linear_Fit_G_inf_mM", "R2_Glycine_Linear_Fit"
      )
      param_values <- c(param_values, 
                        coefs["k_app"], coefs["G_inf"], R2_g_lin
      )
    } else {
      param_names <- c(param_names, 
                       "Glycine_Linear_Fit_k_app_h-1", "Glycine_Linear_Fit_G_inf_mM", "R2_Glycine_Linear_Fit"
      )
      param_values <- c(param_values, 
                        NA, NA, NA
      )
    }
    
    export_params_df <- data.frame(Parameter = param_names, Value = param_values)
    
    # Convert parameter data frame to tab-separated string with a header
    params_header <- "#--- Section 1: Fitted Kinetic Parameters and Model Metrics ---"
    params_string <- paste(
      utils::capture.output(
        write.table(export_params_df, sep = "\t", row.names = FALSE, quote = FALSE)
      ), 
      collapse = "\n"
    )
    results$params_export_text <- paste(params_header, params_string, sep = "\n")
    
    
    # --- Section 2: Raw Data and Fitted Data at Experimental Time Points ---
    export_df <- results$exp_data
    export_df$GG_fit <- GG0_fit * exp(-ksum_fit * export_df$time)
    export_df$cGG_fit <- GG0_fit * k2_fit / ksum_fit * (1 - exp(-ksum_fit * export_df$time)) + pars_est["bg_cGG"]
    export_df$G_fit <- 2 * GG0_fit * k1_fit / ksum_fit * (1 - exp(-ksum_fit * export_df$time)) + pars_est["bg_G"]
    
    data_header <- "\n\n#--- Section 2: Raw Data and Fitted Values (at experimental time points) ---"
    data_string <- paste(
      utils::capture.output(
        write.table(export_df, sep = "\t", row.names = FALSE, quote = FALSE)
      ), 
      collapse = "\n"
    )
    results$data_export_text <- paste(data_header, data_string, sep = "\n")
    
    
    # --- Section 3: High-Resolution Glycine Plotting Data (NEW) ---
    t_plot <- seq(min(exp_data$time), max(exp_data$time), length.out = 200)
    glycine_plot_df <- data.frame(time = t_plot)
    
    # Model-based G curve
    glycine_plot_df$G_Model_Curve <- 2 * GG0_fit * k1_fit / ksum_fit * (1 - exp(-ksum_fit * t_plot)) + bg_G_fit
    
    # Linear Fit G curve
    if (g_fit_model_success) {
      coefs <- coef(results$g_fit_model)
      glycine_plot_df$G_Linear_Curve <- coefs["G_inf"] * (1 - exp(-coefs["k_app"] * t_plot))
    } else {
      glycine_plot_df$G_Linear_Curve <- NA_real_
    }
    
    glycine_plot_header <- "\n\n#--- Section 3: High-Resolution Fitted Curves for Glycine Plotting ---"
    glycine_plot_string <- paste(
      utils::capture.output(
        write.table(glycine_plot_df, sep = "\t", row.names = FALSE, quote = FALSE)
      ), 
      collapse = "\n"
    )
    results$glycine_plot_data_text <- paste(glycine_plot_header, glycine_plot_string, sep = "\n")
    # ------------------------------------------------------------------
  })
  
  # --- Observer to handle the Show Export Data button press ---
  observeEvent(input$show_data_btn, {
    req(results$data_export_text)
    
    # Switch to the 'Export Data' tab
    updateTabsetPanel(session, "main_tabs", selected = "export_data_tab")
    
    # Show a notification
    showNotification("Select the text in the table below and copy manually.", type = "message", duration = 8)
  })
  
  # --- Renderer for the raw export data (Combines all three sections) ---
  output$raw_export_data <- renderPrint({
    req(results$params_export_text, results$data_export_text, results$glycine_plot_data_text)
    
    # Combine parameters, experimental data, and plotting curves
    cat(results$params_export_text, results$data_export_text, results$glycine_plot_data_text, sep = "\n")
  })
  
  # --- All other renderOutputs (fit_results, metrics_results, plots, etc.) remain unchanged for brevity. ---
  
  output$fit_results <- renderPrint({
    req(results$pars_est)
    
    if (is.character(results$fit) && results$fit == "fail") {
      cat("Model fit failed. Please check your data.\n")
      return()
    }
    
    pars_est <- results$pars_est
    GG0_fit <- exp(pars_est["logGG0"])
    k1_fit <- exp(pars_est["logk1"])
    k2_fit <- exp(pars_est["logk2"])
    ksum_fit <- k1_fit + k2_fit 
    
    cat("Fitted Parameters\n")
    cat(sprintf("GG0 = %.4f mM\n", GG0_fit))
    cat(sprintf("k1 (hydrolysis) = %.4f h^-1\n", k1_fit))
    cat(sprintf("k2 (cyclization) = %.4f h^-1\n", k2_fit))
    cat(sprintf("k_sum = %.4f h^-1\n", ksum_fit))
    
    cat("\nHalf-lives:\n")
    cat(sprintf("Total (1/(k1+k2)) t1/2 = %.2f h\n", log(2) / ksum_fit))
    cat(sprintf("Hydrolysis alone t1/2 = %.2f h\n", log(2) / k1_fit))
    cat(sprintf("Cyclization alone t1/2 = %.2f h\n", log(2) / k2_fit))
  })
  
  output$metrics_results <- renderPrint({
    req(results$R2_global)
    
    if (is.character(results$fit) && results$fit == "fail") {
      return()
    }
    
    cat("Goodness of Fit\n")
    cat(sprintf("Global R² = %.4f\n", results$R2_global))
    cat("\nPer-species R²:\n")
    print(results$R2_by_species)
    
    pars_est <- results$pars_est
    cat("\nFitted Baselines:\n")
    cat(sprintf("bg_cGG = %.4f mM\n", pars_est["bg_cGG"]))
    cat(sprintf("bg_G = %.4f mM\n", pars_est["bg_G"]))
  })
  
  output$fit_plot <- renderPlot({
    req(results$long_data, results$pars_est)
    
    if (is.character(results$fit) && results$fit == "fail") {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "Model fit failed. Plot not available.")
      return()
    }
    
    pars_est <- results$pars_est
    GG0_fit <- exp(pars_est["logGG0"])
    k1_fit <- exp(pars_est["logk1"])
    k2_fit <- exp(pars_est["logk2"])
    ksum_fit <- k1_fit + k2_fit 
    bg_cGG_fit <- pars_est["bg_cGG"]
    bg_G_fit <- pars_est["bg_G"]
    
    t_fit <- seq(min(results$exp_data$time), max(results$exp_data$time), length.out = 300)
    
    fit_GG <- GG0_fit * exp(-ksum_fit * t_fit)
    fit_cGG <- GG0_fit * k2_fit / ksum_fit * (1 - exp(-ksum_fit * t_fit)) + bg_cGG_fit
    fit_G <- 2 * GG0_fit * k1_fit / ksum_fit * (1 - exp(-ksum_fit * t_fit)) + bg_G_fit
    
    # Base R Plotting
    cols <- c("#1f78b4", "#e31a1c", "#33a02c") # GG, cGG, G colors
    
    # Set up plot limits
    y_range <- range(results$long_data$conc, fit_GG, fit_cGG, fit_G, na.rm = TRUE)
    
    plot(NA, xlim = range(results$long_data$time), ylim = y_range,
         xlab = "Time (h)", ylab = "Concentration (mM)",
         main = "Kinetic Model Fit (Base R Graphics)", type = "n",
         cex.main = 1.2, font.main = 2)
    
    # Add fit lines
    lines(t_fit, fit_GG, col = cols[1], lwd = 2)
    lines(t_fit, fit_cGG, col = cols[2], lwd = 2)
    lines(t_fit, fit_G, col = cols[3], lwd = 2)
    
    # Add data points
    points(results$exp_data$time, results$exp_data$GG, col = cols[1], pch = 16, cex = 1.5)
    points(results$exp_data$time, results$exp_data$cGG, col = cols[2], pch = 16, cex = 1.5)
    points(results$exp_data$time, results$exp_data$G, col = cols[3], pch = 16, cex = 1.5)
    
    # Add legend
    legend("topright", 
           legend = c("GG (data)", "cGG (data)", "G (data)", "GG (fit)", "cGG (fit)", "G (fit)"),
           col = c(cols, cols),
           pch = c(16, 16, 16, NA, NA, NA),
           lty = c(NA, NA, NA, 1, 1, 1),
           lwd = c(NA, NA, NA, 2, 2, 2),
           cex = 1)
  })
  
  output$resid_plot <- renderPlot({
    req(results$long_data)
    
    if (is.character(results$fit) && results$fit == "fail") {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", main = "Model fit failed. Plot not available.")
      return()
    }
    
    rd <- results$long_data
    species_levels <- levels(rd$species)
    cols <- c("#1f78b4", "#e31a1c", "#33a02c")
    
    # Base R multi-panel plot
    par(mfrow = c(1, length(species_levels)), mar = c(4, 4, 3, 1))
    
    for (i in seq_along(species_levels)) {
      species_data <- rd[rd$species == species_levels[i], ]
      
      plot(species_data$time, species_data$resid,
           xlab = "Time (h)", ylab = "Residual (mM)",
           main = paste("Residuals for", species_levels[i]),
           col = cols[i], pch = 16, 
           cex.main = 1, cex.lab = 1,
           ylim = range(rd$resid)) 
      
      # Add the zero line
      abline(h = 0, lty = 2, col = "gray")
    }
    
    # Reset plotting layout
    par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
  })
  
  output$g_fit_results <- renderPrint({
    if (is.character(results$g_fit_model) && results$g_fit_model == "fail") {
      cat("⚠️ Linear Glycine fit failed.\n")
      cat("Tip: Data may not follow simple exponential rise.")
      return()
    }
    
    coefs <- coef(results$g_fit_model)
    cat("Linear Glycine Fit (apparent kinetics):\n")
    cat(sprintf("k_app = %.4f h^-1\n", coefs["k_app"]))
    cat(sprintf("t1/2 = %.2f h\n", log(2) / coefs["k_app"]))
    cat(sprintf("G_inf = %.2f mM\n", coefs["G_inf"]))
    
    g_data <- results$exp_data[, c("time", "G")]
    g_pred <- predict(results$g_fit_model, g_data)
    ss_res <- sum((g_data$G - g_pred)^2)
    ss_tot <- sum((g_data$G - mean(g_data$G))^2)
    R2_g <- 1 - ss_res / ss_tot
    cat(sprintf("R² = %.4f\n", R2_g))
  })
  
  output$model_g_results <- renderPrint({
    req(results$pars_est, results$exp_data)
    
    if (is.character(results$fit) && results$fit == "fail") {
      cat("⚠️ Model-based Glycine predictions not available due to main model fit failure.\n")
      return()
    }
    
    pars_est <- results$pars_est
    GG0_fit <- exp(pars_est["logGG0"])
    k1_fit <- exp(pars_est["logk1"])
    k2_fit <- exp(pars_est["logk2"])
    ksum_fit <- k1_fit + k2_fit 
    
    cat("Model-Based Glycine Predictions:\n")
    cat(sprintf("k_sum (k1+k2) = %.4f h^-1\n", ksum_fit))
    cat(sprintf("G_inf (as 2*GG0*k1/k_sum) = %.4f mM\n", 2 * GG0_fit * k1_fit / ksum_fit))
    
    g_data <- results$exp_data[, c("time", "G")]
    G_model_pred <- 2 * GG0_fit * k1_fit / ksum_fit * (1 - exp(-ksum_fit * g_data$time)) + pars_est["bg_G"]
    ss_res_model <- sum((g_data$G - G_model_pred)^2)
    ss_tot_g <- sum((g_data$G - mean(g_data$G))^2)
    R2_model_g <- 1 - ss_res_model / ss_tot_g
    cat(sprintf("R² = %.4f\n", R2_model_g))
  })
  
  output$g_fit_plot <- renderPlot({
    req(results$exp_data, results$glycine_plot_data_text) # Req the high-res data here
    
    g_data <- results$exp_data[, c("time", "G")]
    t_vals <- seq(min(g_data$time), max(g_data$time), length.out = 200)
    
    # Extract the high-res data used for export
    # NOTE: In a true shiny app, you would use the reactive `results$glycine_plot_df` if it were stored.
    # Since we only store the text string for export, we'll re-run the calculation or rely on the t_vals vector.
    
    # Recalculate/use existing for plot clarity:
    G_model_pred <- results$model_g_pred
    
    # Base R Plot Setup
    plot(g_data$time, g_data$G,
         xlab = "Time (h)", ylab = "Glycine (mM)",
         main = "Glycine Fit Comparison",
         col = "#33a02c", pch = 16, cex = 1.5,
         ylim = range(0, g_data$G, G_model_pred, na.rm = TRUE))
    
    # Add linear fit line if model is available
    if (!(is.character(results$g_fit_model) && results$g_fit_model == "fail")) {
      coefs <- coef(results$g_fit_model)
      g_lin_pred <- coefs["G_inf"] * (1 - exp(-coefs["k_app"] * t_vals))
      lines(t_vals, g_lin_pred, col = "#e31a1c", lwd = 2, lty = 2) # Dashed line for linear fit
    }
    
    # Add model-based fit line if main fit is successful
    if (!(is.character(results$fit) && results$fit == "fail")) {
      lines(t_vals, G_model_pred, col = "#1f78b4", lwd = 2, lty = 1) # Solid line for model fit
    }
    
    # Add Legend
    legend_text <- c("Data")
    legend_col <- c("#33a02c")
    legend_pch <- c(16)
    legend_lty <- c(NA)
    legend_lwd <- c(NA)
    
    if (!(is.character(results$g_fit_model) && results$g_fit_model == "fail")) {
      legend_text <- c(legend_text, "Linear Fit")
      legend_col <- c(legend_col, "#e31a1c")
      legend_pch <- c(legend_pch, NA)
      legend_lty <- c(legend_lty, 2)
      legend_lwd <- c(legend_lwd, 2)
    }
    
    if (!(is.character(results$fit) && results$fit == "fail")) {
      legend_text <- c(legend_text, "Model Fit")
      legend_col <- c(legend_col, "#1f78b4")
      legend_pch <- c(legend_pch, NA)
      legend_lty <- c(legend_lty, 1)
      legend_lwd <- c(legend_lwd, 2)
    }
    
    legend("bottomright", 
           legend = legend_text,
           col = legend_col,
           pch = legend_pch,
           lty = legend_lty,
           lwd = legend_lwd,
           cex = 1)
  })
  
  output$download_csv <- downloadHandler(
    filename = function() paste0("GG_kinetic_fit_", Sys.Date(), ".csv"),
    content = function(file) {
      showNotification("Download functionality is replaced by the 'Show Export Data' button.", type = "warning")
      return(NULL)
    }
  )
}

shinyApp(ui = ui, server = server)