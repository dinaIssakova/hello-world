
 library(shiny)
 library(rgenesconverged)
 ui <- fluidPage(

   # Application title
   titlePanel("rgenesconverged"),

   # Sidebar with a slider input for number of bins
   sidebarLayout(
    sidebarPanel(

    fileInput("phydat", h3("Upload FASTA alignment.") ),

    textInput("species", h3("Enter name of reference species"),
                    value = "Human"),

    numericInput("pos", h3("Enter position at which to align"),
              value = "1", min=1),

    sliderInput("th", h3("Convergence Score Threshold"),
                     min = 0, max = 20, value = 5),
    radioButtons("type", "Convergence Options", choices = c("Absolute", "Score")),
    checkboxInput("pcheck", "Calculate p if convergent (beta)", FALSE),

    actionButton("submit", "Submit"),
      ),
    mainPanel(

    plotOutput(outputId = "mainPlot")


  )

  )


)


server <- function(input, output, session) {
  plot_reactive <- eventReactive(input$submit, {
    isolate(input$phydat)
    isolate(input$species)
    isolate(input$pos)
    isolate(input$th)
  })

  output$mainPlot <- renderPlot({
    # Plot only when submit is selected.
    plot_reactive()
    #data(BLOSUM62)

    #Import input variables.
    phydat <- phangorn::read.phyDat(isolate(input$phydat$datapath), format="fasta")
    typeIn <- isolate(input$type)
    pos <- isolate(input$pos)
    spe = isolate(input$species)
    threshold <- isolate(input$th)

    # Want to change input for sake of variable naming while providing a more descriptive vers. to user
    if (typeIn == "Absolute"){
      type="abs"
    } else {
      type="score"
    }

    #Build phylogenetic tree
    dm <- phangorn::dist.ml(phydat)
    treeNJ <- ape::nj(dm)

    if (typeIn == "Absolute"){
      type="abs"
    } else {
      type="score"
    }

    # We 'unfold' the rgenesconverged plot and getConvergent here to be able to insert a progress bar.
    withProgress(message = 'Making plot', value = 0, {

    # List of species from leaves of tree
    species <- treeNJ$tip.label
    convSpe <- c(spe)
    for (s in species){
      # Increment progress bar for every species evaluated for convergence.
      print(s)
      incProgress(1/length(species), detail = paste("Evaluating species", s))

      # Verify if that species is convergent
      cond = areCondSatisfied(treeNJ, phydat, s, spe, pos, type, threshold=threshold, BLOSUM62)
      if (!rapportools::is.boolean(cond)){
        cond = FALSE
      }


      if (s != spe && cond){
        convSpe = c(convSpe, s)
      }
      if (s != spe && input$pcheck){
        p = probOfSiteConfig(treeNJ, phydat, s, spe, pos)
        print(sprintf("Species %s is potentially convergent with p %e.", s, p))
      }

    }


    # Form groups based on species in the convergent subgroup, and outside it.
    groupInfo <- split(treeNJ$tip.label, treeNJ$tip.label%in%convSpe)

    tree <- ggtree::groupOTU(treeNJ, groupInfo)
    library(ggtree)
    # Guangchuang Yu, David Smith, Huachen Zhu, Yi Guan, Tommy Tsan-Yuk Lam.
    # ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates
    # and other associated data. Methods in Ecology and Evolution 2017, 8(1):28-36, doi:10.1111/2041-210X.12628
    ggtree::ggtree(tree, ggtree::aes(color=group)) + ggtree::geom_text2(aes(subset=!isTip, label=node), hjust=-.3) +
      ggtree::geom_tiplab()
    })
    })
}

#END
shinyApp(ui = ui, server = server)
