library(tidyverse)
source('cat_script.R')

server <- function(input, output) {
    output$sim_plot <- renderPlot({

        # get inputs
        harv.prop.init <- input$harv.prop.init # initial harvest rate
        harv.prop.maint <- input$harv.prop.maint # ongoing harvest
        iter <- input$iter # replicate sims
        tt <- input$tt  # time horizon (years)

        # run simulations
        n.sums.mat = cat_sim(harv.prop.init, harv.prop.maint, iter, tt)

        # plots results
        sims <- n.sums.mat %>%
            as_tibble() %>%
            mutate(rep = 1:iter) %>%
            pivot_longer(-rep) %>%
            mutate(year = as.numeric(str_extract(name, "\\d+")))
        sum_sims <- sims %>% 
            group_by(year) %>% 
            summarise(value = mean(value))
 
        ggplot() +
            geom_line(data=sims, aes(year, value, group = rep), alpha = 0.1) +
            geom_line(data=sum_sims, aes(year, value), size=1.5, linetype=2) + 
            ylab('Cat population size') + 
            theme_bw()

    })
}
