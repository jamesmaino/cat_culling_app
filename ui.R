# Define UI for app that draws a histogram ----
ui <- fluidPage(

  # App title ----
  titlePanel("CAT-BE-DEAD"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      # Input: Sliders
      sliderInput(inputId = "harv.prop.init",
                  label = "Proportion of cats culled (initial):",
                  min = 0,
                  max = 1,
                  value = 0.5, step = 0.01),
      sliderInput(inputId = "harv.prop.maint",
                  label = "Proportion of cats culled (ongoing):",
                  min = 0,
                  max = 1,
                  value = 0.1, step = 0.01),
      sliderInput(inputId = "iter",
                  label = "Number of replicate simulations:",
                  min = 1,
                  max = 1000,
                  value = 100, step = 10),
      sliderInput(inputId = "tt",
                  label = "Time horizon (years):",
                  min = 0,
                  max = 100,
                  value = 10, step = 1)            

    ),

    # Main panel for displaying outputs ----
    mainPanel(

      # Output: Plot
      plotOutput(outputId = "sim_plot")

    )
  )
)