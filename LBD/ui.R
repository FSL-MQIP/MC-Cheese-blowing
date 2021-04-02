library(shiny)

# Define UI for application
shinyUI(fluidPage(

    # Application title
    titlePanel("Prediction of cheese late blowing"),

    # Sidebar with checkbox for prevention strategies
    sidebarLayout(
        sidebarPanel(
            sliderInput("days", "How long is the ripening time?", 30, 120, value=60),
            sliderInput("temp", "What is the ripening temperature?",11,15, value=13),
            numericInput("count", "What is the mean spore concentration (log10 MPN/L) in raw milk?", value = 1.85),
            h4("Select intervention strategies"),
            checkboxInput("mf","Microfiltration/bactofugation", value=FALSE),
            checkboxInput("nitrate","Addition of nitrate at 2.5g/ 100L milk", value = FALSE),
            checkboxInput("lysozyme","Addition of lysozyme at 2.5g/ 100L milk", value = FALSE),
            checkboxInput("lab","Use of bacteriocinogenic LAB at 0.3% level", value = FALSE),
            submitButton("Submit")
        ),

        # Show a histgram of predicted concentration and a textbox with predicted proportion of late blown cheese
        mainPanel(
            plotOutput("hist"),
            h5("Predicted proportion of cheese with LBD"),
            textOutput("prop")
        )
    )
))
