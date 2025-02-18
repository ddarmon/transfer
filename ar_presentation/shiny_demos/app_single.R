# Written using `o3-mini-high` on 3 February 2025

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

# Set up parallel processing
plan(multisession)

# --- Function to compute the Accuracy Ratio via Mann–Whitney U statistic ---
calc_AR_MW <- function(scores, outcomes) {
  d <- sum(outcomes)         # number of defaults
  n <- length(outcomes)      # total sample size
  ranks <- rank(scores)      # compute ranks
  R_def <- sum(ranks[outcomes == 1])
  U <- R_def - d * (d + 1) / 2
  AR <- 2 * U / (d * (n - d)) - 1
  return(AR)
}

# --- Simulation Function (parallelized via future_lapply) ---
simulate_AR <- function(n, d, sims,
                        alpha_default, beta_default,
                        alpha_nondefault, beta_nondefault) {

  outcomes <- c(rep(1, d), rep(0, n - d))

  # For each simulation, compute AR and its 95% CI (via pROC)
  results <- future_lapply(1:sims, function(i) {
    # Generate scores for defaults and nondefaults
    scores_defaults <- rbeta(d, shape1 = alpha_default, shape2 = beta_default)
    scores_nondef   <- rbeta(n - d, shape1 = alpha_nondefault, shape2 = beta_nondefault)
    scores <- c(scores_defaults, scores_nondef)

    # Compute AR using our Mann–Whitney based function
    AR_sim <- calc_AR_MW(scores, outcomes)

    # Compute ROC (AUC) and its analytical CI via pROC (using the DeLong method)
    roc_obj <- roc(response = outcomes,
                   predictor = scores,
                   levels = c(0, 1),
                   direction = "<",
                   quiet = TRUE)
    ci_auc <- ci.auc(roc_obj, conf.level = 0.95, method = "delong")

    # Transform the CI for AUC into a CI for AR (using AR = 2*AUC - 1)
    AR_ci_lower <- 2 * ci_auc[1] - 1
    AR_ci_est   <- 2 * ci_auc[2] - 1
    AR_ci_upper <- 2 * ci_auc[3] - 1

    c(AR = AR_sim, lower = AR_ci_lower, upper = AR_ci_upper)
  }, future.seed = TRUE)

  res_mat <- do.call(rbind, results)
  return(as.data.frame(res_mat))
}

# --- Shiny UI using bslib for theming ---
ui <- fluidPage(
  # Use a Bootswatch theme (here: "minty")
  theme = bs_theme(bootswatch = "minty"),

  titlePanel("Monte Carlo Simulation of Accuracy Ratio with Population CAP Curve"),

  sidebarLayout(
    sidebarPanel(
      numericInput("n", "Total Sample Size (n):", value = 500, min = 1, step = 1),
      numericInput("d", "Number of Defaults (d):", value = 50, min = 1, step = 1),
      numericInput("sims", "Number of Simulations:", value = 1000, min = 1, step = 1),
      hr(),
      fluidRow(
        column(width = 6,
               h4("Default"),
               numericInput("alpha_default", "Alpha:", value = 7, min = 0, step = 10),
               numericInput("beta_default",  "Beta:",  value = 100, min = 0, step = 10)
        ),
        column(width = 6,
               h4("Nondefault"),
               numericInput("alpha_nondefault", "Alpha:", value = 4, min = 0, step = 10),
               numericInput("beta_nondefault",  "Beta:",  value = 100, min = 0, step = 10)
        )
      )
    ),

    mainPanel(
      # Top row: left shows the beta density plot, right shows the CAP curve.
      fluidRow(
        column(width = 6, plotOutput("densityPlot")),
        column(width = 6, plotOutput("capPlot"))
      ),
      br(),
      # Bottom row: AR histogram and CI plot side-by-side.
      fluidRow(
        column(width = 6, plotOutput("histPlot")),
        column(width = 6, plotlyOutput("ciPlot"))
      )
    )
  )
)

# --- Shiny Server ---
server <- function(input, output, session) {

  # Reactive simulation: any change in inputs re-runs the simulation.
  simResults <- reactive({
    n_val    <- input$n
    d_val    <- input$d
    sims_val <- input$sims
    a_def    <- input$alpha_default
    b_def    <- input$beta_default
    a_non    <- input$alpha_nondefault
    b_non    <- input$beta_nondefault

    # Input validations
    if (d_val > n_val) {
      showNotification("Number of defaults cannot exceed total sample size.", type = "error")
      return(NULL)
    }
    if (d_val <= 0) {
      showNotification("Number of defaults must be greater than zero.", type = "error")
      return(NULL)
    }

    # Run the Monte Carlo simulation for AR values and their confidence intervals.
    sim_data <- simulate_AR(n = n_val, d = d_val, sims = sims_val,
                            alpha_default = a_def, beta_default = b_def,
                            alpha_nondefault = a_non, beta_nondefault = b_non)

    # Exact Calculation of Population AUC and Its Variance (DeLong's Method)

    # Compute the exact (population) AUC:
    # AUC = ∫[0,1] f_def(x) * F_non(x) dx, where
    # f_def = dbeta(x, a_def, b_def) and F_non = pbeta(x, a_non, b_non)
    AUC_pop <- integrate(function(x) {
      dbeta(x, a_def, b_def) * pbeta(x, a_non, b_non)
    }, lower = 0, upper = 1)$value
    impliedAR <- 2 * AUC_pop - 1

    # Compute the exact variance components:
    # For defaults, V(x) = F_non(x), so
    sigma2_V <- integrate(function(x) {
      dbeta(x, a_def, b_def) * (pbeta(x, a_non, b_non) - AUC_pop)^2
    }, lower = 0, upper = 1)$value

    # For nondefaults, note that W(y) = 1 - F_def(y); equivalently,
    # we use F_def(y) and center it around AUC (since AUC = E[F_non(x)] = 1 - E[W])
    sigma2_W <- integrate(function(y) {
      dbeta(y, a_non, b_non) * (pbeta(y, a_def, b_def) - AUC_pop)^2
    }, lower = 0, upper = 1)$value

    n1 <- d_val             # defaults
    n0 <- n_val - d_val     # nondefaults

    var_auc_exact <- sigma2_V / n1 + sigma2_W / n0
    se_auc_exact <- sqrt(var_auc_exact)
    se_AR_exact <- 2 * se_auc_exact  # because AR = 2*AUC - 1

    # Approximate Calculation of Variance Using Hanley & McNeil's Formulas
    # under the assumption of binormality of the score distributions.
    #
    # Note: The following relationships (with Q1 and Q2) are approximations.
    Q1 <- AUC_pop / (2 - AUC_pop)
    Q2 <- 2 * AUC_pop^2 / (1 + AUC_pop)
    var_auc_approx <- (AUC_pop * (1 - AUC_pop) +
                         (n1 - 1) * (Q1 - AUC_pop^2) +
                         (n0 - 1) * (Q2 - AUC_pop^2)) / (n1 * n0)
    se_auc_approx <- sqrt(var_auc_approx)
    se_AR_approx <- 2 * se_auc_approx

    list(sim_data = sim_data,
         impliedAR = impliedAR,
         se_AR_exact = se_AR_exact,
         se_AR_approx = se_AR_approx,
         a_def = a_def,
         b_def = b_def,
         a_non = a_non,
         b_non = b_non)
  })

  # --- Population Density Plot (Static) ---
  output$densityPlot <- renderPlot({
    res <- simResults()
    if (is.null(res)) return(NULL)

    # Create a fine grid of x values between 0 and 1
    x_grid <- seq(0, 1, length.out = 2001)

    # Compute the exact Beta densities for defaults and nondefaults
    dens_default    <- dbeta(x_grid, res$a_def, res$b_def)
    dens_nondefault <- dbeta(x_grid, res$a_non, res$b_non)

    # Build a data frame for plotting
    df_dens <- data.frame(
      x = rep(x_grid, 2),
      density = c(dens_default, dens_nondefault),
      Group = factor(rep(c("Default", "Nondefault"), each = length(x_grid)))
    )

    ggplot(df_dens, aes(x = x, y = density, color = Group)) +
      geom_line(size = 1) +
      scale_color_manual(values = c("Default" = "red", "Nondefault" = "blue")) +
      labs(title = "Population Density of Scores",
           x = "Score", y = "Density") +
      theme_minimal()
  })

  # --- Population-level CAP Curve Plot ---
  output$capPlot <- renderPlot({
    res <- simResults()
    if (is.null(res)) return(NULL)

    # Base rate (fraction of defaults)
    p_default <- input$d / input$n
    # Create a grid of thresholds (from high to low scores)
    t_grid <- seq(1, 0, length.out = 1000)

    # For a threshold t, the cumulative fraction of the overall population with score ≥ t is:
    x <- 100 * (p_default * (1 - pbeta(t_grid, res$a_def, res$b_def)) +
                  (1 - p_default) * (1 - pbeta(t_grid, res$a_non, res$b_non)))
    # The fraction of defaults captured is:
    y <- 100 * (1 - pbeta(t_grid, res$a_def, res$b_def))
    df_cap <- data.frame(Population = x, Defaults = y)

    ggplot(df_cap, aes(x = Population, y = Defaults)) +
      geom_line(size = 1, color = "darkblue") +
      geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey") +
      labs(title = "Population-level CAP Curve",
           x = "Cumulative % of Population",
           y = "Cumulative % of Defaults Captured") +
      theme_minimal() +
      scale_x_continuous(limits = c(0, 100)) +
      scale_y_continuous(limits = c(0, 100))
  })

  # --- AR Histogram Plot (Static) ---
  output$histPlot <- renderPlot({
    res <- simResults()
    if (is.null(res)) return(NULL)

    sim_data <- res$sim_data
    AR_vals <- sim_data$AR
    # Compute bin width using: 2 * IQR / n^(1/3)
    bin_width <- 2 * IQR(AR_vals) / (length(AR_vals))^(1/3)
    df_AR <- data.frame(AR = AR_vals)

    # Compute the 2.5th and 97.5th percentiles
    quantiles <- quantile(AR_vals, probs = c(0.025, 0.975))

    # Create a data frame for the implied normal density
    x_vals <- seq(-1, 1, length.out = 1000)
    normal_df <- data.frame(
      x = x_vals,
      density = dnorm(x_vals, mean = res$impliedAR, sd = res$se_AR_exact)
    )

    # Build the plot, mapping a constant label for each density curve.
    ggplot() +
      # Empirical sampling density from histogram (using geom_step)
      geom_step(data = df_AR,
                aes(x = AR, y = after_stat(density),
                    color = "Empirical Sampling Density"),
                stat = "bin", direction = "mid", binwidth = bin_width) +
      # Asymptotic normal density line
      geom_line(data = normal_df,
                aes(x = x, y = density,
                    color = "Asymptotic Normal Density"),
                size = 1, linetype = "solid") +
      # Add vertical lines (these can be added without affecting the legend)
      geom_vline(xintercept = res$impliedAR, color = "red", linetype = "dashed", size = 1) +
      geom_vline(xintercept = quantiles, color = "blue", linetype = "dotted", size = 1) +
      # Define manual colors and legend title
      scale_color_manual(
        name = "Density",
        values = c(
          "Empirical Sampling Density" = "black",
          "Asymptotic Normal Density" = "purple"
        ),
        breaks = c("Empirical Sampling Density", "Asymptotic Normal Density")
      ) +
      labs(title = glue("Sampling Dist. of AR Statistic, Pop. AR = {round(res$impliedAR, 3)}"),
           x = "Accuracy Ratio (AR)", y = "Density") +
      xlim(c(-1, 1)) +
      theme_minimal()
  })

  # --- Confidence Interval Plot (Interactive) ---
  output$ciPlot <- renderPlotly({
    res <- simResults()
    if (is.null(res)) return(NULL)

    sim_data <- res$sim_data
    impliedAR <- res$impliedAR
    # For each simulation, mark whether the computed CI covers the implied AR.
    sim_data <- sim_data %>%
      mutate(covered = (lower <= impliedAR & upper >= impliedAR),
             sim = 1:n())

    # Compute the empirical coverage rate
    coverage_rate <- mean(sim_data$covered)

    # Include the empirical coverage probability in the plot title
    title_str <- glue("CI for AR - Emp. Coverage: {round(coverage_rate * 100, 1)}%")

    p <- ggplot(sim_data, aes(x = sim)) +
      geom_errorbar(aes(ymin = lower, ymax = upper, color = covered), width = 0.5) +
      geom_point(aes(y = AR, color = covered), size = 1) +
      geom_hline(yintercept = impliedAR, linetype = "solid", color = "blue") +
      scale_color_manual(values = c("TRUE" = "forestgreen", "FALSE" = "red"),
                         labels = c("FALSE" = "Not Covered", "TRUE" = "Covered")) +
      ylim(c(-1, 1)) +
      labs(title = title_str,
           x = "Simulation Index",
           y = "Accuracy Ratio (AR)",
           color = "CI Coverage") +
      theme_minimal()

    ggplotly(p)
  })
}

# --- Launch the App ---
shinyApp(ui, server)
