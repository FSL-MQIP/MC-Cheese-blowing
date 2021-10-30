library(shiny)
library(shinycssloaders)

# Define UI for application
shinyUI(fluidPage(

    # Application title
    #titlePanel("Prediction of cheese late blowing"),

    # Sidebar with checkbox for prevention strategies
    sidebarLayout(
        sidebarPanel(
            sliderInput("temp", "What is the ripening temperature?",11,15, value=13),
            numericInput("count", "What is the mean spore concentration (log10 MPN/L) in raw milk?", value = 1.85),
            h4("Select intervention strategies"),
            checkboxInput("mf","Microfiltration (reduces 98% spores)", value=FALSE),
            checkboxInput("bf","Bactofugation (reduces 98% spores)", value=FALSE),
            checkboxInput("nitrate","Addition of nitrate at 2.5g/ 100L milk", value = FALSE),
            checkboxInput("lysozyme","Addition of lysozyme at 2.5g/ 100L milk", value = FALSE),
            checkboxInput("lab","Use of bacteriocinogenic LAB at 0.3% level", value = FALSE),
            #h4("Sensitivity analysis"),
            #sliderInput("mu_change","Proportional change in optimal growth rate", -0.4,0.4,value=0),
            #sliderInput("cf_change", "Absolute change in concentration factor", -2,2,value=0),
            #sliderInput("tl_change","Absolute change (log10) in threshold level", -0.5,0.5,value=0),
            #sliderInput("pH_change","Absolute change in cheese pH", -0.2,0.2,value = 0),
            submitButton("Submit"),
        ),

        # Show a histgram of predicted concentration and a textbox with predicted proportion of late blown cheese
        mainPanel(
            withSpinner(plotOutput("plot", width = 480, height = 240)),
            h5("Predicted proportion of cheese with LBD"), 
            verbatimTextOutput("prop"),
            plotOutput("base_plot", width = 480, height = 240),
            h5("Predicted proportion of cheese with LBD (base model)"), 
            verbatimTextOutput("base_prop"),
        )
    )
))
