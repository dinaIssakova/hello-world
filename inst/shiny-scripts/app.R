
# library(shiny)
# ui <- fluidPage(
#   #tags$head(tags$script(src = "message-handler.js")),
#   # Application title
#   titlePanel("rgenesconverged"),
#
#   # Sidebar with a slider input for number of bins
#   sidebarLayout(
#     sidebarPanel(
#
#     fileInput("tree", h3("File input: tree") ),
#     fileInput("phydat", h3("File input: phydat object") ),
#
#     textInput("species", h3("Enter name of reference species"),
#                     value = "Human"),
#
#     numericInput("pos", h3("Enter position at which to align"),
#               value = "1", min=1),
#
#     sliderInput("t", h3("Convergence Score Threshold"),
#                      min = 0, max = 20, value = 5),
#
#     actionButton("submit", "Submit"),
#       ),
#   mainPanel(
#    # withSpinner(plotlyOutput("mainPlot")),
#     plotOutput(outputId = "mainPlot")
#
#
#   )
#
#   )
#
#
# )
#
#
# server <- function(input, output) {
#   randomVals <- eventReactive(input$submit, {
#     runif(input$tree && input$phydat)
#   })
#   observeEvent(input$submit, {
#     print("Clicked!")
#   })
#   output$mainPlot <- renderPlot({
#
#     #rgenesconverged::rgenesconvergedPlot(input$tree, input$phydat, input$species, input$pos, input$t)
#     convNodes <- getConvergent(tree, phydat, spe, pos, t)
#     convNodes <- c(convNodes, spe)
#
#     groupInfo <- split(tree$tip.label, tree$tip.label%in%convNodes)
#     tree <- ggtree::groupOTU(tree, groupInfo)
#
#     # Guangchuang Yu, David Smith, Huachen Zhu, Yi Guan, Tommy Tsan-Yuk Lam.
#     # ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates
#     # and other associated data. Methods in Ecology and Evolution 2017, 8(1):28-36, doi:10.1111/2041-210X.12628
#     ggtree::ggtree(tree, aes(color=group)) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()
#
#
#   })
# }
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# #END
# shinyApp(ui = ui, server = server)
