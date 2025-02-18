# Written using `o3-mini-high` on 4 February 2025

# Literature review on AUC inference:
#
# https://chatgpt.com/share/67b09816-cf7c-800c-bf54-5538836039fb

# Load necessary libraries
library(shiny)
library(ggplot2)
library(dplyr)
library(glue)
library(future)
library(future.apply)
library(bslib)
library(pROC)
library(plotly)
library(MASS)
library(cubature)  # For numerical integration
library(mvtnorm)  # For multivariate normal density

# Set up parallel processing
plan(multisession)

# --- Function to simulate the AR difference for two correlated ROC curves ---
simulate_AR_diff <- function(n, d, sims, rho,
                             alpha_def1, beta_def1, alpha_def2, beta_def2,
                             alpha_non1, beta_non1, alpha_non2, beta_non2) {
  outcomes <- c(rep(1, d), rep(0, n - d))

  results <- future_lapply(1:sims, function(i) {
    # For defaults (positives): generate paired scores via bivariate normal,
    # then transform to uniforms and finally to beta-distributed scores.
    Sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
    defaults_norm <- MASS::mvrnorm(n = d, mu = c(0, 0), Sigma = Sigma)
    defaults_unif <- pnorm(defaults_norm)
    scores_def1 <- qbeta(defaults_unif[, 1], shape1 = alpha_def1, shape2 = beta_def1)
    scores_def2 <- qbeta(defaults_unif[, 2], shape1 = alpha_def2, shape2 = beta_def2)

    # For nondefaults (negatives)
    n_non <- n - d
    nondefaults_norm <- MASS::mvrnorm(n = n_non, mu = c(0, 0), Sigma = Sigma)
    nondefaults_unif <- pnorm(nondefaults_norm)
    scores_non1 <- qbeta(nondefaults_unif[, 1], shape1 = alpha_non1, shape2 = beta_non1)
    scores_non2 <- qbeta(nondefaults_unif[, 2], shape1 = alpha_non2, shape2 = beta_non2)

    # Combine scores for each model
    scores1 <- c(scores_def1, scores_non1)
    scores2 <- c(scores_def2, scores_non2)

    # Compute ROC curves and AR for each model using pROC.
    roc1 <- roc(response = outcomes,
                predictor = scores1,
                levels = c(0, 1),
                direction = "<",
                quiet = TRUE)
    roc2 <- roc(response = outcomes,
                predictor = scores2,
                levels = c(0, 1),
                direction = "<",
                quiet = TRUE)
    AR1 <- 2 * as.numeric(auc(roc1)) - 1
    AR2 <- 2 * as.numeric(auc(roc2)) - 1
    AR_diff <- AR1 - AR2

    # Compute CI for the difference in AUC using the paired DeLong method.
    test_res <- roc.test(roc1, roc2, paired = TRUE, method = "delong")
    if (!is.null(test_res$conf.int)) {
      ci_lower_auc <- test_res$conf.int[1]
      ci_upper_auc <- test_res$conf.int[2]
      AR_ci_lower <- 2 * ci_lower_auc
      AR_ci_upper <- 2 * ci_upper_auc
    } else {
      AR_ci_lower <- NA
      AR_ci_upper <- NA
    }

    c(AR_diff = AR_diff, lower = AR_ci_lower, upper = AR_ci_upper)
  }, future.seed = TRUE)

  res_mat <- do.call(rbind, results)
  return(as.data.frame(res_mat))
}

# --- Shiny UI using bslib for theming ---
ui <- fluidPage(
  # Use a Bootswatch theme (here: "minty")
  theme = bs_theme(bootswatch = "minty"),

  titlePanel("Monte Carlo Simulation of AR Difference from Correlated ROC Curves"),

  sidebarLayout(
    sidebarPanel(
      numericInput("n", "Total Sample Size (n):", value = 500, min = 1, step = 1),
      numericInput("d", "Number of Defaults (d):", value = 50, min = 1, step = 1),
      numericInput("sims", "Number of Simulations:", value = 500, min = 1, step = 1),
      numericInput("rho", "Correlation (rho):", value = 0.7, min = -1, max = 1, step = 0.05),
      hr(),
      h4("Defaults (Positives)"),
      fluidRow(
        column(6,
               numericInput("alpha_def1", "Model 1 Alpha:", value = 7, min = 0, step = 5)
        ),
        column(6,
               numericInput("beta_def1", "Model 1 Beta:", value = 100, min = 0, step = 5)
        )
      ),
      fluidRow(
        column(6,
               numericInput("alpha_def2", "Model 2 Alpha:", value = 7, min = 0, step = 5)
        ),
        column(6,
               numericInput("beta_def2", "Model 2 Beta:", value = 105, min = 0, step = 5)
        )
      ),
      hr(),
      h4("Nondefaults (Negatives)"),
      fluidRow(
        column(6,
               numericInput("alpha_non1", "Model 1 Alpha:", value = 4, min = 0, step = 5)
        ),
        column(6,
               numericInput("beta_non1", "Model 1 Beta:", value = 100, min = 0, step = 5)
        )
      ),
      fluidRow(
        column(6,
               numericInput("alpha_non2", "Model 2 Alpha:", value = 4, min = 0, step = 5)
        ),
        column(6,
               numericInput("beta_non2", "Model 2 Beta:", value = 100, min = 0, step = 5)
        )
      )
    ),

    mainPanel(
      # Top row: Density plot on left and CAP curves on right
      fluidRow(
        column(width = 6, plotlyOutput("densityPlot")),
        column(width = 6, plotOutput("capPlot"))
      ),
      br(),
      # Bottom row: AR difference histogram and CI plot side-by-side
      fluidRow(
        column(width = 6, plotOutput("histPlot")),
        column(width = 6, plotlyOutput("ciPlot"))
      )
    )
  )
)

# --- Shiny Server ---
server <- function(input, output, session) {

  simResults <- reactive({
    n_val    <- input$n
    d_val    <- input$d
    sims_val <- input$sims
    rho_val  <- input$rho

    # Parameters for Defaults for Model 1 and Model 2
    a_def1 <- input$alpha_def1
    b_def1 <- input$beta_def1
    a_def2 <- input$alpha_def2
    b_def2 <- input$beta_def2

    # Parameters for Nondefaults for Model 1 and Model 2
    a_non1 <- input$alpha_non1
    b_non1 <- input$beta_non1
    a_non2 <- input$alpha_non2
    b_non2 <- input$beta_non2

    if(d_val > n_val) {
      showNotification("Number of defaults cannot exceed total sample size.", type = "error")
      return(NULL)
    }
    if(d_val <= 0) {
      showNotification("Number of defaults must be greater than zero.", type = "error")
      return(NULL)
    }

    # Run the Monte Carlo simulation for AR difference and its CI.
    sim_data <- simulate_AR_diff(n = n_val, d = d_val, sims = sims_val, rho = rho_val,
                                 alpha_def1 = a_def1, beta_def1 = b_def1,
                                 alpha_def2 = a_def2, beta_def2 = b_def2,
                                 alpha_non1 = a_non1, beta_non1 = b_non1,
                                 alpha_non2 = a_non2, beta_non2 = b_non2)

    # Compute the implied (population) AUC and AR for each model by integration.
    auc1_implied <- integrate(function(x) {
      dbeta(x, a_def1, b_def1) * pbeta(x, a_non1, b_non1)
    }, lower = 0, upper = 1)$value
    AR1_implied <- 2 * auc1_implied - 1

    auc2_implied <- integrate(function(x) {
      dbeta(x, a_def2, b_def2) * pbeta(x, a_non2, b_non2)
    }, lower = 0, upper = 1)$value
    AR2_implied <- 2 * auc2_implied - 1

    impliedAR_diff <- AR1_implied - AR2_implied

    # --- Compute the "exact" standard error for AR difference using numerical integration ---
    n1 <- d_val
    n0 <- n_val - d_val

    # Variance for the V-components (defaults)
    f_V_diff_sq <- function(u) {
      u1 <- u[1]
      u2 <- u[2]
      x1 <- qbeta(u1, a_def1, b_def1)
      x2 <- qbeta(u2, a_def2, b_def2)
      V1 <- pbeta(x1, a_non1, b_non1)
      V2 <- pbeta(x2, a_non2, b_non2)
      diff <- V1 - V2
      z1 <- qnorm(u1)
      z2 <- qnorm(u2)
      Sigma <- matrix(c(1, rho_val, rho_val, 1), nrow = 2)
      joint_dens <- dmvnorm(c(z1, z2), mean = c(0, 0), sigma = Sigma)
      copula_density <- joint_dens / (dnorm(z1) * dnorm(z2))
      diff^2 * copula_density
    }
    val_V <- adaptIntegrate(f_V_diff_sq, lowerLimit = c(0, 0), upperLimit = c(1, 1))
    E_V_diff_sq <- val_V$integral
    var_V_diff <- E_V_diff_sq - (auc1_implied - auc2_implied)^2

    # Variance for the W-components (nondefaults)
    f_W_diff_sq <- function(u) {
      u1 <- u[1]
      u2 <- u[2]
      x1 <- qbeta(u1, a_non1, b_non1)
      x2 <- qbeta(u2, a_non2, b_non2)
      W1 <- pbeta(x1, a_def1, b_def1)
      W2 <- pbeta(x2, a_def2, b_def2)
      diff <- W1 - W2
      z1 <- qnorm(u1)
      z2 <- qnorm(u2)
      Sigma <- matrix(c(1, rho_val, rho_val, 1), nrow = 2)
      joint_dens <- dmvnorm(c(z1, z2), mean = c(0, 0), sigma = Sigma)
      copula_density <- joint_dens / (dnorm(z1) * dnorm(z2))
      diff^2 * copula_density
    }
    val_W <- adaptIntegrate(f_W_diff_sq, lowerLimit = c(0, 0), upperLimit = c(1, 1))
    E_W_diff_sq <- val_W$integral
    # Compute expected W values (note: E[W] = ∫ F_def(x) f_nondef(x) dx)
    E_W1 <- integrate(function(x) dbeta(x, a_non1, b_non1) * pbeta(x, a_def1, b_def1),
                      lower = 0, upper = 1)$value
    E_W2 <- integrate(function(x) dbeta(x, a_non2, b_non2) * pbeta(x, a_def2, b_def2),
                      lower = 0, upper = 1)$value
    var_W_diff <- E_W_diff_sq - (E_W1 - E_W2)^2

    var_auc_diff <- var_V_diff / n1 + var_W_diff / n0
    se_auc_diff_exact <- sqrt(var_auc_diff)
    se_AR_diff_exact <- 2 * se_auc_diff_exact  # because AR = 2*AUC - 1

    list(sim_data = sim_data,
         impliedAR_diff = impliedAR_diff,
         AR1_implied = AR1_implied,
         AR2_implied = AR2_implied,
         se_AR_diff_exact = se_AR_diff_exact,
         a_def1 = a_def1, b_def1 = b_def1,
         a_def2 = a_def2, b_def2 = b_def2,
         a_non1 = a_non1, b_non1 = b_non1,
         a_non2 = a_non2, b_non2 = b_non2)
  })

  # --- Density Plot: Show population densities for Model 1 and Model 2 ---
  output$densityPlot <- renderPlotly({
    res <- simResults()
    if(is.null(res)) return(NULL)

    x_grid <- seq(0, 1, length.out = 2001)

    # For Model 1: Defaults and Nondefaults
    dens_def1 <- dbeta(x_grid, res$a_def1, res$b_def1)
    dens_non1 <- dbeta(x_grid, res$a_non1, res$b_non1)
    df1 <- data.frame(
      x = rep(x_grid, 2),
      density = c(dens_def1, dens_non1),
      Group = rep(c("Default", "Nondefault"), each = length(x_grid)),
      Model = "Model 1"
    )

    # For Model 2: Defaults and Nondefaults
    dens_def2 <- dbeta(x_grid, res$a_def2, res$b_def2)
    dens_non2 <- dbeta(x_grid, res$a_non2, res$b_non2)
    df2 <- data.frame(
      x = rep(x_grid, 2),
      density = c(dens_def2, dens_non2),
      Group = rep(c("Default", "Nondefault"), each = length(x_grid)),
      Model = "Model 2"
    )

    df_all <- rbind(df1, df2)

    p <- ggplot(df_all, aes(x = x, y = density, color = Group, linetype = Model)) +
      geom_line(size = 0.5) +
      labs(title = "Population Densities of Scores",
           x = "Score", y = "Density") +
      theme_minimal() +
      theme(
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6)
      ) +
      scale_color_manual(values = c("Default" = "red", "Nondefault" = "blue"))

    ggplotly(p) %>% layout(
      legend = list(
        x = 0.60,
        y = 1.00,
        xanchor = "left",
        yanchor = "top",
        bgcolor = 'rgba(255,255,255,1)'
      )
    )
  })

  # --- CAP Curve Plot ---
  output$capPlot <- renderPlot({
    res <- simResults()
    if(is.null(res)) return(NULL)

    p_default <- input$d / input$n
    t_grid <- seq(1, 0, length.out = 1000)

    # Model 1 CAP
    x1 <- 100 * (p_default * (1 - pbeta(t_grid, res$a_def1, res$b_def1)) +
                   (1 - p_default) * (1 - pbeta(t_grid, res$a_non1, res$b_non1)))
    y1 <- 100 * (1 - pbeta(t_grid, res$a_def1, res$b_def1))
    df_cap1 <- data.frame(x = x1, y = y1, Model = "Model 1")

    # Model 2 CAP
    x2 <- 100 * (p_default * (1 - pbeta(t_grid, res$a_def2, res$b_def2)) +
                   (1 - p_default) * (1 - pbeta(t_grid, res$a_non2, res$b_non2)))
    y2 <- 100 * (1 - pbeta(t_grid, res$a_def2, res$b_def2))
    df_cap2 <- data.frame(x = x2, y = y2, Model = "Model 2")

    df_cap <- rbind(df_cap1, df_cap2)

    ggplot(df_cap, aes(x = x, y = y, color = Model)) +
      geom_line(size = 1) +
      geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey") +
      labs(title = "Population-level CAP Curves",
           x = "Cumulative % of Population",
           y = "Cumulative % of Defaults Captured") +
      theme_minimal() +
      scale_x_continuous(limits = c(0, 100)) +
      scale_y_continuous(limits = c(0, 100)) +
      scale_color_manual(values = c("Model 1" = "#009E73", "Model 2" = "#CC79A7"))
  })

  # --- AR Difference Histogram with Exact Asymptotic Normal Density ---
  output$histPlot <- renderPlot({
    res <- simResults()
    if(is.null(res)) return(NULL)

    sim_data <- res$sim_data
    AR_diff_vals <- sim_data$AR_diff
    bin_width <- 2 * IQR(AR_diff_vals) / (length(AR_diff_vals))^(1/3)
    df_diff <- data.frame(AR_diff = AR_diff_vals)
    quantiles <- quantile(AR_diff_vals, probs = c(0.025, 0.975))

    title_str <- glue("Sampling Distribution of AR Difference\n(Pop. ARs: P1 = {round(res$AR1_implied, 3)}, P2 = {round(res$AR2_implied, 3)}; Diff = {round(res$impliedAR_diff, 3)})")

    # Use the exact standard error from the integration-based computation for the asymptotic normal density.
    x_vals <- seq(-1, 1, length.out = 1000)
    normal_df <- data.frame(
      x = x_vals,
      density = dnorm(x_vals, mean = res$impliedAR_diff, sd = res$se_AR_diff_exact)
    )

    ggplot() +
      geom_step(data = df_diff,
                aes(x = AR_diff, y = after_stat(density), color = "Empirical Sampling Density"),
                stat = "bin", direction = "mid", binwidth = bin_width) +
      geom_line(data = normal_df,
                aes(x = x, y = density, color = "Asymptotic Normal Density"),
                size = 1, linetype = "solid") +
      geom_vline(xintercept = res$impliedAR_diff, color = "red", linetype = "dashed", size = 1) +
      geom_vline(xintercept = quantiles, color = "blue", linetype = "dotted", size = 1) +
      scale_color_manual(
        name = "Density",
        values = c("Empirical Sampling Density" = "black",
                   "Asymptotic Normal Density" = "purple")
      ) +
      labs(title = title_str,
           x = "Difference in Accuracy Ratio (AR1 - AR2)", y = "Density") +
      xlim(c(-1, 1)) +
      theme_minimal()
  })

  # --- Confidence Interval Plot (Interactive) ---
  output$ciPlot <- renderPlotly({
    res <- simResults()
    if(is.null(res)) return(NULL)

    sim_data <- res$sim_data
    impliedAR_diff <- res$impliedAR_diff
    sim_data <- sim_data %>%
      mutate(covered = (lower <= impliedAR_diff & upper >= impliedAR_diff),
             sim = 1:n())

    coverage_rate <- mean(sim_data$covered, na.rm = TRUE)
    power_rate <- mean(sim_data$lower > 0 | sim_data$upper < 0, na.rm = TRUE)

    title_str <- glue("CI for AR Difference\nEmpirical Coverage: {round(coverage_rate * 100, 1)}%\nEmpirical Power: {round(power_rate * 100, 1)}%")

    p <- ggplot(sim_data, aes(x = sim)) +
      geom_errorbar(aes(ymin = lower, ymax = upper, color = covered), width = 0.5) +
      geom_point(aes(y = AR_diff, color = covered), size = 1) +
      geom_hline(yintercept = impliedAR_diff, linetype = "solid", color = "blue") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      scale_color_manual(values = c("TRUE" = "forestgreen", "FALSE" = "red"),
                         labels = c("FALSE" = "Not Covered", "TRUE" = "Covered")) +
      ylim(c(-1, 1)) +
      labs(title = title_str,
           x = "Simulation Index",
           y = "Difference in Accuracy Ratio (AR1 - AR2)",
           color = "CI Coverage") +
      theme_minimal()

    ggplotly(p)
  })
}

# --- Launch the App ---
shinyApp(ui, server)
