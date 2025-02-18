# Created with `o3-mini-high` on 15 February 2025.

###############################################################################
# Accuracy Ratio Difference Precision & Power Calculator
#
# This Shiny application allows users to analyze pilot study data for
# comparing two predictive models via their score outputs. It computes
# pilot study statistics, estimates the number of required events, implied
# confidence interval half-widths, and the statistical power for testing the
# difference in Accuracy Ratios (AR).
#
# Users can either upload pilot data (CSV with columns: score1, score2, outcome)
# or simulate data if none is provided.
#
###############################################################################

library(shiny)
library(bslib)
library(MASS)     # For simulating data
library(plotly)   # For interactive plots

# Increase the maximum allowed file upload size to 100 MB.
options(shiny.maxRequestSize = 100 * 1024^2)

###############################################################################
# Helper Functions
###############################################################################

# Adaptive root-finding function with additional error-checking.
# Attempts to bracket the root by expanding the upper bound until a sign change
# in f(x) is observed. If a valid bracket is found, uniroot() is used.
find_root_adaptive <- function(f, lower, initial_upper, tol = 1e-6, max_iter = 10) {
  # Check if the lower bound is already close enough to a root.
  if (!is.na(f(lower)) && abs(f(lower)) < tol) return(lower)

  x_lower <- lower
  x_upper <- initial_upper
  iter <- 0

  # Expand the upper bound until a sign change is found or max iterations reached.
  while (iter < max_iter && (is.na(f(x_lower)) || is.na(f(x_upper)) || f(x_lower) * f(x_upper) > 0)) {
    x_upper <- x_upper * 2
    iter <- iter + 1
  }

  # If a proper bracket wasn't found, return the lower bound.
  if (is.na(f(x_lower)) || is.na(f(x_upper)) || f(x_lower) * f(x_upper) > 0) {
    return(lower)
  }

  # Use uniroot() to find the root within the bracket.
  return(uniroot(f, lower = x_lower, upper = x_upper, tol = tol)$root)
}

# Compute DeLong components for two scoring models.
# This function computes the AUC for each model, derives the Accuracy Ratio (AR),
# and computes differences in the V and W components needed for subsequent variance
# and SE calculations.
compute_delong_components <- function(score1, score2, outcome) {
  # Identify indexes for events (defaults) and non-events.
  pos_idx <- which(outcome == 1)
  neg_idx <- which(outcome == 0)

  n1 <- length(pos_idx)
  n0 <- length(neg_idx)

  # Compute pairwise comparisons between events and non-events.
  I1 <- outer(score1[pos_idx], score1[neg_idx],
              FUN = function(x, y) as.numeric(x > y) + 0.5 * as.numeric(x == y))
  I2 <- outer(score2[pos_idx], score2[neg_idx],
              FUN = function(x, y) as.numeric(x > y) + 0.5 * as.numeric(x == y))

  # Calculate AUCs.
  AUC1 <- mean(I1)
  AUC2 <- mean(I2)

  # Compute Accuracy Ratios (AR).
  AR1 <- 2 * AUC1 - 1
  AR2 <- 2 * AUC2 - 1

  # Compute per-case deviations needed for variance estimation.
  V1 <- rowMeans(I1) - AUC1
  V2 <- rowMeans(I2) - AUC2
  W1 <- colMeans(I1) - AUC1
  W2 <- colMeans(I2) - AUC2

  # Difference in V and W components between models.
  V_diff <- V1 - V2
  W_diff <- W1 - W2

  list(
    AUC1 = AUC1, AUC2 = AUC2,
    AR1 = AR1, AR2 = AR2,
    delta_AR = AR1 - AR2,
    V_diff = V_diff, W_diff = W_diff,
    n1 = n1, n0 = n0
  )
}

# Calculate required number of events (defaults) to achieve a desired
# half-width for the AR difference confidence interval.
sample_size_calc <- function(comp, desired_halfwidth, conf_level = 0.95) {
  z <- qnorm(1 - (1 - conf_level) / 2)
  n1_pilot <- comp$n1
  n0_pilot <- comp$n0

  # Compute variance estimates for V_diff and W_diff.
  var_V_diff <- if (n1_pilot > 1) var(comp$V_diff) else 0
  var_W_diff <- if (n0_pilot > 1) var(comp$W_diff) else 0

  # Calculate the pilot SE for the AR difference.
  se_AR_diff_pilot <- sqrt(4 * (var_V_diff / n1_pilot + var_W_diff / n0_pilot))

  output_text <- c()
  output_text <- c(output_text, "Pilot study results:")
  output_text <- c(output_text, sprintf("  Model 1 AR: %.3f", comp$AR1))
  output_text <- c(output_text, sprintf("  Model 2 AR: %.3f", comp$AR2))
  output_text <- c(output_text, sprintf("  AR difference: %.3f", comp$AR1 - comp$AR2))
  output_text <- c(output_text, sprintf("  Pilot SE for AR difference: %.4f", se_AR_diff_pilot))
  output_text <- c(output_text, sprintf("  n_events (defaults) in pilot: %d", n1_pilot))
  output_text <- c(output_text, sprintf("  n_non-events in pilot: %d", n0_pilot))
  output_text <- c(output_text, "")

  # Compute denominator for the required events formula.
  denominator <- (desired_halfwidth / (2 * z))^2 - var_W_diff / n0_pilot
  if (denominator <= 0) {
    stop("Desired half-width is too small given the variability among non-events. Increase the half-width or the non-event sample size.")
  }

  n1_required <- ceiling(var_V_diff / denominator)
  total_required <- n0_pilot + n1_required

  output_text <- c(output_text, sprintf("For a desired half-width of %.3f on the AR difference CI:", desired_halfwidth))
  output_text <- c(output_text, sprintf("  Required number of events (defaults): %d", n1_required))
  output_text <- c(output_text, sprintf("  Total required sample size (with %d non-events fixed): %d", n0_pilot, total_required))

  list(
    text = paste(output_text, collapse = "\n"),
    n_events_required = n1_required,
    total_required = total_required,
    pilot_AR1 = comp$AR1,
    pilot_AR2 = comp$AR2,
    pilot_AR_diff = comp$AR1 - comp$AR2,
    se_AR_diff_pilot = se_AR_diff_pilot,
    var_V_diff = var_V_diff,
    var_W_diff = var_W_diff,
    n0 = n0_pilot,
    z = z
  )
}

# Compute the implied half-width for the AR difference CI based on a planned
# number of events.
calc_implied_halfwidth <- function(comp, planned_events, conf_level = 0.95) {
  z <- qnorm(1 - (1 - conf_level) / 2)
  n1_pilot <- comp$n1
  n0_pilot <- comp$n0

  var_V_diff <- if (n1_pilot > 1) var(comp$V_diff) else 0
  var_W_diff <- if (n0_pilot > 1) var(comp$W_diff) else 0

  implied_halfwidth <- 2 * z * sqrt(var_V_diff / planned_events + var_W_diff / n0_pilot)

  output_text <- c()
  output_text <- c(output_text, "Pilot study results:")
  output_text <- c(output_text, sprintf("  Model 1 AR: %.3f", comp$AR1))
  output_text <- c(output_text, sprintf("  Model 2 AR: %.3f", comp$AR2))
  output_text <- c(output_text, sprintf("  AR difference: %.3f", comp$AR1 - comp$AR2))
  se_AR_diff_pilot <- sqrt(4 * (var_V_diff / n1_pilot + var_W_diff / n0_pilot))
  output_text <- c(output_text, sprintf("  Pilot SE for AR difference: %.4f", se_AR_diff_pilot))
  output_text <- c(output_text, sprintf("  n_events (defaults) in pilot: %d", n1_pilot))
  output_text <- c(output_text, sprintf("  n_non-events in pilot: %d", n0_pilot))
  output_text <- c(output_text, "")
  output_text <- c(output_text, sprintf("For a planned %d events in the study:", planned_events))
  output_text <- c(output_text, sprintf("  Implied half-width for the AR difference CI: %.3f", implied_halfwidth))

  list(
    text = paste(output_text, collapse = "\n"),
    implied_halfwidth = implied_halfwidth,
    var_V_diff = var_V_diff,
    var_W_diff = var_W_diff,
    n0 = n0_pilot,
    z = z
  )
}

# Calculate the statistical power for testing H0: AR1 = AR2
# using the pilot study variances and a planned number of events.
calc_power <- function(comp, planned_events, effect_size = NA, conf_level = 0.95) {
  z <- qnorm(1 - (1 - conf_level) / 2)
  n0 <- comp$n0
  n1_pilot <- comp$n1

  var_V_diff <- if (n1_pilot > 1) var(comp$V_diff) else 0
  var_W_diff <- if (n0 > 1) var(comp$W_diff) else 0

  SE <- 2 * sqrt(var_V_diff / planned_events + var_W_diff / n0)
  pilot_effect <- comp$AR1 - comp$AR2
  if (is.na(effect_size)) {
    effect_size <- pilot_effect
  }

  # Noncentrality parameter and computed power.
  nc <- effect_size / SE
  power <- pnorm(nc - z) + 1 - pnorm(nc + z)

  output_text <- c()
  output_text <- c(output_text, "Pilot study results:")
  output_text <- c(output_text, sprintf("  Model 1 AR: %.3f", comp$AR1))
  output_text <- c(output_text, sprintf("  Model 2 AR: %.3f", comp$AR2))
  output_text <- c(output_text, sprintf("  AR difference (pilot estimate): %.3f", pilot_effect))
  output_text <- c(output_text, sprintf("  Using effect size: %.3f", effect_size))

  if (abs(effect_size - pilot_effect) > 1e-8) {
    output_text <- c(output_text,
                     "WARNING: You are overriding the effect size. This assumes that the variance remains fixed as the true mean shifts, which may not hold in practice.")
  }

  se_pilot <- sqrt(4 * (if(n1_pilot > 1) var(comp$V_diff)/n1_pilot else 0 +
                          if(n0 > 1) var(comp$W_diff)/n0 else 0))
  output_text <- c(output_text, sprintf("  Pilot SE for AR difference: %.4f", se_pilot))
  output_text <- c(output_text, sprintf("  n_events (defaults) in pilot: %d", comp$n1))
  output_text <- c(output_text, sprintf("  n_non-events in pilot: %d", n0))
  output_text <- c(output_text, "")
  output_text <- c(output_text, sprintf("For a planned %d events:", planned_events))
  output_text <- c(output_text, sprintf("  Computed SE for AR difference: %.4f", SE))
  output_text <- c(output_text, sprintf("  Noncentral parameter: %.3f", nc))
  output_text <- c(output_text, sprintf("  Significance level: %.3f", 1 - conf_level))
  output_text <- c(output_text, sprintf("  Power for testing H0: AR1 = AR2: %.2f%%", power * 100))

  list(
    text = paste(output_text, collapse = "\n"),
    power = power,
    SE = SE,
    nc = nc,
    z = z
  )
}

# Compute the empirical CAP (Cumulative Accuracy Profile) curve.
# This function orders the scores, computes the cumulative sum of defaults,
# and returns a data frame with the cumulative fraction of defaults.
compute_cap_curve <- function(score, outcome) {
  ord <- order(score, decreasing = TRUE)
  outcome_sorted <- outcome[ord]
  cum_defaults <- cumsum(outcome_sorted)
  total_defaults <- sum(outcome_sorted)
  if (total_defaults == 0) total_defaults <- 1  # Avoid division by zero.
  cap <- cum_defaults / total_defaults
  data.frame(population_fraction = seq_along(score) / length(score), cap = cap)
}

###############################################################################
# Define the UI
###############################################################################
ui <- fluidPage(
  theme = bs_theme(bootswatch = "flatly"),
  titlePanel("Accuracy Ratio Difference Precision & Power Calculator"),

  sidebarLayout(
    sidebarPanel(
      # Option to display pilot statistics in the main panel.
      checkboxInput("show_pilot_statistics", "Show Pilot Statistics", value = TRUE),

      # Tabset for different calculation modes.
      tabsetPanel(
        type = "pills",
        id = "calc_mode",
        tabPanel(
          "Required Events",
          value = "halfwidth",
          numericInput("desired_halfwidth", "Desired Half-Width for AR Difference CI",
                       value = 0.05, min = 0.001, step = 0.005)
        ),
        tabPanel(
          "Implied Half-Width",
          value = "n_events",
          numericInput("planned_events", "Planned Number of Events", value = 150, min = 1)
        ),
        tabPanel(
          "Power Calculation",
          value = "power",
          numericInput("planned_events_power", "Planned Number of Events for Power Calculation",
                       value = 150, min = 1),
          numericInput("user_effect_size", "Effect Size (AR Difference) [Override Pilot]",
                       value = NA, step = 0.01)
        )
      ),

      HTML("<hr>"),
      numericInput("conf_level", "Confidence Level", value = 0.95, min = 0.80, max = 0.99, step = 0.01),
      HTML("<hr>"),

      # File input for pilot data.
      fileInput("pilot_file", "Upload Pilot Data (CSV)",
                accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
      helpText(HTML("CSV should include columns:<br><ul>
        <li>score1</li>
        <li>score2</li>
        <li>outcome (1 = default, 0 = non-default)</li>
      </ul>")),

      HTML("<hr>"),
      # Option to simulate pilot data if no file is uploaded.
      checkboxInput("simulate", "Simulate pilot data if no file is uploaded", value = TRUE),
      conditionalPanel(
        condition = "input.simulate == true && input.pilot_file == null",
        numericInput("n_total", "Total Pilot Sample Size", value = 10000, min = 100),
        numericInput("n_defaults", "Number of Defaults (Events)", value = 100, min = 1),
        numericInput("score_corr", "Score Correlation", value = 0.9, min = 0, max = 1, step = 0.01),
        downloadButton("downloadData", "Download Simulated Data")
      ),
      HTML("<hr>"),
      # Option to show CAP curves.
      checkboxInput("show_cap", "Show Pilot CAP Curves", value = FALSE),
      conditionalPanel(
        condition = "input.show_cap == true",
        br(),
        plotlyOutput("capPlot", height = "400px")
      )
    ),

    mainPanel(
      conditionalPanel(
        condition = "input.show_pilot_statistics == true",
        verbatimTextOutput("results")
      ),
      br(),
      plotlyOutput("plot", height = "400px"),
      conditionalPanel(
        condition = "input.calc_mode == 'power'",
        br(),
        plotlyOutput("plot2", height = "400px")
      )
    )
  )
)

###############################################################################
# Define the Server Logic
###############################################################################
server <- function(input, output, session) {

  # Reactive function to load or simulate pilot data.
  # Progress indicators are used to signal long-running operations.
  pilotData <- reactive({
    withProgress(message = "Loading pilot data...", value = 0, {
      incProgress(0.2, detail = "Checking for uploaded file...")

      if (!is.null(input$pilot_file)) {
        # Attempt to read the uploaded CSV file.
        incProgress(0.3, detail = "Reading uploaded data...")
        tryCatch({
          df <- read.csv(input$pilot_file$datapath)
          req(all(c("score1", "score2", "outcome") %in% names(df)))
          incProgress(0.5, detail = "Data loaded successfully.")
          return(df)
        }, error = function(e) {
          showNotification("Error reading uploaded file.", type = "error")
          return(NULL)
        })
      } else if (input$simulate) {
        # Simulate pilot data if no file is provided.
        incProgress(0.2, detail = "Simulating pilot data...")
        n_total <- input$n_total
        n_defaults <- input$n_defaults
        n_nondefaults <- n_total - n_defaults
        outcome <- c(rep(1, n_defaults), rep(0, n_nondefaults))
        score_corr <- input$score_corr
        Sigma <- matrix(c(1, score_corr, score_corr, 1), nrow = 2)
        # Simulate scores for defaults and non-defaults.
        scores_defaults <- mvrnorm(n_defaults, mu = c(1, 1), Sigma = Sigma)
        scores_nondefaults <- mvrnorm(n_nondefaults, mu = c(0, 0), Sigma = Sigma)
        score1 <- c(scores_defaults[, 1], scores_nondefaults[, 1])
        score2 <- c(scores_defaults[, 2], scores_nondefaults[, 2])
        df <- data.frame(score1 = score1, score2 = score2, outcome = outcome)
        incProgress(0.5, detail = "Simulation complete.")
        return(df)
      } else {
        return(NULL)
      }
    })
  })

  # Reactive function to compute pilot statistics.
  # Wrapped in withProgress in case the computation becomes heavy.
  pilotComps <- reactive({
    df <- pilotData()
    validate(need(!is.null(df), "Pilot data is not available."))
    withProgress(message = "Computing pilot statistics...", value = 0.5, {
      comps <- compute_delong_components(df$score1, df$score2, df$outcome)
      incProgress(0.5)
      comps
    })
  })

  # Update the effect size input with the pilot AR difference,
  # unless the user has manually overridden it.
  observe({
    req(pilotComps())
    default_effect <- pilotComps()$AR1 - pilotComps()$AR2
    updateNumericInput(session, "user_effect_size", value = default_effect)
  })

  # A simple summary for plotting purposes.
  pilotSummary <- reactive({
    comps <- pilotComps()
    n1 <- comps$n1
    n0 <- comps$n0
    z <- qnorm(1 - (1 - input$conf_level) / 2)
    var_V <- if(n1 > 1) var(comps$V_diff) else 0
    var_W <- if(n0 > 1) var(comps$W_diff) else 0
    list(var_V = var_V, var_W = var_W, n0 = n0, z = z)
  })

  # Render the results (text output) based on the chosen calculation mode.
  output$results <- renderText({
    comps <- pilotComps()
    res <- NULL
    if (input$calc_mode == "halfwidth") {
      res <- tryCatch({
        sample_size_calc(comps, desired_halfwidth = input$desired_halfwidth, conf_level = input$conf_level)
      }, error = function(e) {
        list(text = paste("Error:", e$message))
      })
    } else if (input$calc_mode == "n_events") {
      res <- tryCatch({
        calc_implied_halfwidth(comps, planned_events = input$planned_events, conf_level = input$conf_level)
      }, error = function(e) {
        list(text = paste("Error:", e$message))
      })
    } else if (input$calc_mode == "power") {
      res <- tryCatch({
        calc_power(
          comp = comps,
          planned_events = input$planned_events_power,
          effect_size = input$user_effect_size,
          conf_level = input$conf_level
        )
      }, error = function(e) {
        list(text = paste("Error:", e$message))
      })
    }
    res$text
  })

  # Render the main interactive plot based on the chosen calculation mode.
  output$plot <- renderPlotly({
    df <- pilotData()
    validate(need(!is.null(df), "No pilot data available for plotting."))

    ps <- pilotSummary()
    z <- ps$z
    var_V <- ps$var_V
    var_W <- ps$var_W
    n0 <- ps$n0

    if (input$calc_mode == "halfwidth") {
      desired_hw <- input$desired_halfwidth
      hw_min <- max(2 * z * sqrt(var_W / n0) + 0.001, desired_hw * 0.5)
      hw_max <- desired_hw * 1.5
      hw_seq <- seq(hw_min, hw_max, length.out = 100)
      events_vals <- sapply(hw_seq, function(hw) {
        denom <- (hw / (2 * z))^2 - var_W / n0
        if (denom > 0) var_V / denom else NA
      })
      denom_ref <- (desired_hw / (2 * z))^2 - var_W / n0
      req_events_ref <- if (denom_ref > 0) var_V / denom_ref else NA

      p <- plot_ly() %>%
        add_lines(x = hw_seq, y = events_vals,
                  line = list(color = "steelblue"),
                  name = "Required Events") %>%
        layout(
          xaxis = list(title = "CI Half-Width for AR Difference",
                       fixedrange = TRUE,
                       range = c(hw_min, hw_max)),
          yaxis = list(title = "Required Number of Events (Defaults)",
                       fixedrange = TRUE)
        )

      shapes <- list(
        list(type = "line", x0 = desired_hw, x1 = desired_hw,
             y0 = min(events_vals, na.rm = TRUE), y1 = max(events_vals, na.rm = TRUE),
             line = list(color = "darkgreen", dash = "dash")),
        list(type = "line", x0 = hw_min, x1 = hw_max,
             y0 = req_events_ref, y1 = req_events_ref,
             line = list(color = "purple", dash = "dash"))
      )
      annotations <- list(
        list(
          x = desired_hw,
          y = max(events_vals, na.rm = TRUE),
          text = paste("Desired Half-Width =", desired_hw),
          showarrow = FALSE,
          xanchor = "center",
          yanchor = "top",
          font = list(color = "darkgreen"),
          textangle = 90,
          xshift = 10
        ),
        list(
          x = hw_max * 0.98,
          y = req_events_ref,
          text = paste("Required Events =", round(req_events_ref, 0)),
          showarrow = FALSE,
          xanchor = "right",
          yanchor = "bottom",
          font = list(color = "purple")
        )
      )
      p <- p %>% layout(shapes = shapes, annotations = annotations)

    } else if (input$calc_mode == "n_events") {
      planned_events <- input$planned_events
      x_min <- 1
      x_max <- max(planned_events * 1.5, 50)
      events_seq <- seq(x_min, x_max, length.out = 2001)
      hw_vals <- 2 * z * sqrt(var_V / events_seq + var_W / n0)

      p <- plot_ly() %>%
        add_lines(x = events_seq, y = hw_vals,
                  line = list(color = "steelblue"),
                  name = "Implied CI Half-Width") %>%
        layout(
          xaxis = list(title = "Number of Events (Defaults)",
                       fixedrange = TRUE,
                       range = c(x_min, x_max)),
          yaxis = list(title = "CI Half-Width for AR Difference",
                       fixedrange = TRUE)
        )

      shapes <- list(
        list(type = "line", x0 = planned_events, x1 = planned_events,
             y0 = min(hw_vals), y1 = max(hw_vals),
             line = list(color = "darkgreen", dash = "dash")),
        list(type = "line", x0 = x_min, x1 = x_max,
             y0 = 2 * z * sqrt(var_V / planned_events + var_W / n0),
             y1 = 2 * z * sqrt(var_V / planned_events + var_W / n0),
             line = list(color = "darkred", dash = "dash"))
      )
      annotations <- list(
        list(
          x = planned_events, y = max(hw_vals),
          text = paste("Planned Events =", planned_events),
          showarrow = FALSE,
          xanchor = "center", yanchor = "top",
          font = list(color = "darkgreen"),
          textangle = 90, xshift = 10
        ),
        list(
          x = x_max * 0.98, y = 2 * z * sqrt(var_V / planned_events + var_W / n0),
          text = paste("Implied Half-Width =", round(2 * z * sqrt(var_V / planned_events + var_W / n0), 3)),
          showarrow = FALSE,
          xanchor = "right", yanchor = "bottom",
          font = list(color = "darkred")
        )
      )
      p <- p %>% layout(shapes = shapes, annotations = annotations)

    } else if (input$calc_mode == "power") {
      planned <- input$planned_events_power
      comp <- pilotComps()
      delta_AR <- input$user_effect_size

      # Define a function to compute power for a given number of events.
      power_function <- function(n_evt) {
        SE <- 2 * sqrt(var_V / n_evt + var_W / n0)
        nc <- delta_AR / SE
        pnorm(nc - z) + 1 - pnorm(nc + z)
      }

      # Find the number of events that gives 80% power.
      n80 <- find_root_adaptive(function(n) power_function(n) - 0.8, lower = 1, initial_upper = planned * 2)

      if (!is.na(n80)) {
        x_max <- max(planned * 1.5, 50, n80 * 1.1)
      } else {
        x_max <- max(planned * 1.5, 50)
      }
      x_min <- 1
      events_seq <- seq(x_min, x_max, length.out = 2001)
      power_vals <- sapply(events_seq, power_function)

      SE_planned <- 2 * sqrt(var_V / planned + var_W / n0)
      nc_planned <- delta_AR / SE_planned
      power_at_planned <- pnorm(nc_planned - z) + 1 - pnorm(nc_planned + z)

      p <- plot_ly() %>%
        add_lines(x = events_seq, y = power_vals,
                  line = list(color = "steelblue"),
                  name = "Computed Power") %>%
        layout(
          xaxis = list(title = "Number of Events (Defaults)",
                       fixedrange = TRUE,
                       range = c(x_min, x_max)),
          yaxis = list(title = "Power for AR Difference Test",
                       fixedrange = TRUE,
                       tickformat = ".0%")
        )

      shapes <- list(
        list(type = "line", x0 = planned, x1 = planned,
             y0 = min(power_vals), y1 = max(power_vals),
             line = list(color = "darkgreen", dash = "dash")),
        list(type = "line", x0 = x_min, x1 = x_max,
             y0 = power_at_planned, y1 = power_at_planned,
             line = list(color = "darkred", dash = "dash"))
      )
      annotations <- list(
        list(
          x = planned, y = max(power_vals),
          text = paste("Planned Events =", planned),
          showarrow = FALSE,
          xanchor = "center", yanchor = "top",
          font = list(color = "darkgreen"),
          textangle = 90, xshift = 10
        ),
        list(
          x = x_max * 0.98, y = power_at_planned,
          text = paste("Computed Power =", round(power_at_planned * 100, 1), "%"),
          showarrow = FALSE,
          xanchor = "right", yanchor = "bottom",
          font = list(color = "darkred")
        )
      )
      p <- p %>% layout(shapes = shapes, annotations = annotations)
    }
    p
  })

  # Reactive expressions for CAP curves based on pilot data.
  capData <- reactive({
    df <- pilotData()
    req(df)
    list(
      score1 = compute_cap_curve(df$score1, df$outcome),
      score2 = compute_cap_curve(df$score2, df$outcome)
    )
  })

  # Render the CAP curve plot.
  output$capPlot <- renderPlotly({
    req(pilotData())
    df <- pilotData()
    cap <- capData()

    # Compute the optimal CAP curve.
    n_defaults <- sum(df$outcome)
    total <- nrow(df)
    optimal <- data.frame(
      population_fraction = c(0, n_defaults / total, 1),
      cap = c(0, 1, 1)
    )

    p <- plot_ly() %>%
      add_lines(data = cap$score1, x = ~population_fraction, y = ~cap,
                name = "score1 CAP", line = list(color = "blue")) %>%
      add_lines(data = cap$score2, x = ~population_fraction, y = ~cap,
                name = "score2 CAP", line = list(color = "red")) %>%
      add_lines(x = c(0, 1), y = c(0, 1),
                name = "Random Model", line = list(color = "black", dash = "dot")) %>%
      add_lines(data = optimal, x = ~population_fraction, y = ~cap,
                name = "Optimal CAP", line = list(color = "green", dash = "dash")) %>%
      layout(
        xaxis = list(title = "P(S > s)", fixedrange = TRUE),
        yaxis = list(title = "P(S > s | D = 1)", fixedrange = TRUE, range = c(0, 1))
      )
    p
  })

  # Additional power plot: Power vs. Effect Size.
  output$plot2 <- renderPlotly({
    if (input$calc_mode != "power") return(NULL)

    comp <- pilotComps()
    planned <- input$planned_events_power
    n0 <- comp$n0
    n1_pilot <- comp$n1
    var_V_diff <- if(n1_pilot > 1) var(comp$V_diff) else 0
    var_W_diff <- if(n0 > 1) var(comp$W_diff) else 0
    z <- qnorm(1 - (1 - input$conf_level) / 2)
    SE <- 2 * sqrt(var_V_diff / planned + var_W_diff / n0)

    # Define power as a function of effect size.
    power_effect <- function(eff) {
      nc <- eff / SE
      pnorm(nc - z) + 1 - pnorm(nc + z)
    }

    # Find the effect size that would give 80% power.
    effect80 <- find_root_adaptive(function(e) power_effect(e) - 0.8, lower = 0, initial_upper = 1)

    if (!is.na(effect80)) {
      range_max <- effect80 * 1.1
    } else {
      pilot_effect <- comp$AR1 - comp$AR2
      range_max <- max(abs(pilot_effect)*2, 0.1)
    }

    effect_seq <- seq(-range_max, range_max, length.out = 2001)
    power_vals <- sapply(effect_seq, power_effect)

    p2 <- plot_ly() %>%
      add_lines(x = effect_seq, y = power_vals,
                line = list(color = "darkorange"),
                name = "Power vs. Effect Size") %>%
      layout(
        xaxis = list(title = "Effect Size (AR Difference)", fixedrange = TRUE,
                     range = c(-range_max, range_max)),
        yaxis = list(title = "Power", fixedrange = TRUE, tickformat = ".0%")
      )

    current_eff <- input$user_effect_size
    current_power <- power_effect(current_eff)
    p2 <- p2 %>% layout(
      shapes = list(
        list(type = "line", x0 = current_eff, x1 = current_eff,
             y0 = min(power_vals), y1 = max(power_vals),
             line = list(color = "blue", dash = "dash"))
      ),
      annotations = list(
        list(x = current_eff, y = max(power_vals),
             text = paste("Current effect size =", round(current_eff, 3)),
             showarrow = FALSE, xanchor = "center", yanchor = "top", font = list(color = "blue"))
      )
    )
    p2
  })

  # Handler to download simulated pilot data.
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("simulated_pilot_data_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      df <- pilotData()
      if (!is.null(df)) {
        write.csv(df, file, row.names = FALSE)
      }
    }
  )
}

# Run the Shiny application.
shinyApp(ui = ui, server = server)
