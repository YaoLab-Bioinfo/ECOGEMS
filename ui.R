
shinyUI(
  fluidPage(
    titlePanel(title=div(
      img(src="headerN.png"),
      span("ECOGEMS:", style = "font-size:36px;color:white;"), 
      span("genotypes of 2058 rice accessions at 8,584,244 SNP sites", style = "font-size:28px;color:white;"),
      style = "background-color:#0073B7;margin-left: -15px;margin-right: -15px;margin-top: -20px;margin-bottom: -10px;"
    ), windowTitle = "Welcome to ECOGEMS!"
    ),
    
    includeCSS("www/footer.css"),
    
    shinydisconnect::disconnectMessage(
      text = "Your session timed out, reload the application!",
      refresh = "Reload now",
      background = "#f89f43",
      colour = "white",
      overlayColour = "grey",
      overlayOpacity = 0.75,
      top = 250,
      refreshColour = "brown"
    ),
    
    navbarPage(
      title = "", 
      #    theme = shinytheme("darkly"), 
      windowTitle = "efficient compression of genotype matrix",
      
      ## About
      tabPanel(title = HTML("<strong style='font-size:20px'>Home</strong>"), icon = icon("home"), includeMarkdown("About.md")),
      
      # Genome browser
      tabPanel(
        title = HTML("<strong style='font-size:20px'>Browse</strong>"), icon = icon("folder-open"),
        
        tags$head(tags$script(HTML('Shiny.addCustomMessageHandler("jsCode",function(message) {eval(message.value);});')),
                  tags$style("
                 input[type='file'] {width:5em;}
                 .toggleButton {width:100%;}
                 .clearButton {float:right; font-size:12px;}
                 .fa-angle-down:before, .fa-angle-up:before {float:right;}
                 .popover{text-align:left;width:500px;background-color:#000000;}
                 .popover-title{color:#FFFFFF;font-size:16px;background-color:#000000;border-color:#000000;}
                 .jhr{display: inline; vertical-align: top; padding-left: 10px;}

                 #sidebarPanel_1 {width:25em;}
                 #mainPanel_1 {left:28em; position:absolute; min-width:27em;}
                 .popover{max-width: 60%;}
              "),
                  tags$style(HTML(".shiny-output-error-validation {color: red;}")),
                  tags$style(HTML(
                    ".checkbox {margin: 0}
                 .checkbox p {margin: 0;}
                 .shiny-input-container {margin-bottom: 0;}
                 .navbar-default .navbar-brand {color: black; font-size:150%;}
                 .navbar-default .navbar-nav > li > a {color:black; font-size:120%;}
                 .shiny-input-container:not(.shiny-input-container-inline) {width: 100%;}
               ")),
                  tags$script(HTML('Shiny.addCustomMessageHandler("jsCode",function(message) {eval(message.value);});')),
                  
                  tags$style(
                    HTML(
                      "
            #inputs-table {
            border-collapse: collapse;
            }
            
            #inputs-table td {
            padding: 3px;
            vertical-align: bottom;
            }

            .multicol .shiny-options-group{
                            -webkit-column-count: 2; /* Chrome, Safari, Opera */
            -moz-column-count: 2;    /* Firefox */
            column-count: 2;
            -moz-column-fill: balanced;
            -column-fill: balanced;
            }
            .checkbox{
            margin-top: 0px !important;
            -webkit-margin-after: 1px !important; 
            }

            "
                    ) #/ HTML
                  ) #/ style
        ), #/ head
        
        fluidRow(column(12,
                        class = "col-md-5",
                        style = "margin: 1px 1px 1px 1px",
                        tags$table(id = "inputs-table",
                                   style = "width: 100%",
                                   tags$tr(
                                     tags$td(style = "width: 40%",
                                             textInput("regB", label = tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Genomic region:</font>'),
                                                                          bsButton("q6", label="", icon=icon("question"), style="info", size="small")),
                                                       value = ""),
                                             
                                             bsPopover("q6", "A genomic region can be determined by chromosome positions or gene locus. For example, chr07:29611303-29669223 or LOC_Os11g35500.",
                                                       trigger = "focus")
                                     ), #/ column 1
                                     tags$td(style = "width: 40%; text-align: right",
                                             div(class = "form-group shiny-input-container",
                                                 shinysky::actionButton("submit1", strong("Submit!",
                                                                                          bsButton("q7", label="", icon=icon("question"), style="info", size="small")
                                                 ), width = "90%", styleclass = "success"),
                                                 conditionalPanel(condition="input.submit1 != '0'", shinysky::busyIndicator(HTML("<div style='color:red;font-size:30px'>Calculation In progress...</div>"), wait = 0)),
                                                 bsPopover("q7", "Whenever the genomic region is updated, please click Submit!",
                                                           trigger = "focus")
                                             )
                                     ), #/ column 2
                                     tags$td(style = "width: 80%; text-align: right",
                                             div(class = "form-group shiny-input-container",
                                                 shinysky::actionButton("clearGB", strong("Reset"), styleclass = "warning"),
                                                 shinysky::actionButton("GBExam", strong("Load example"), styleclass = "info")
                                             )
                                     )
                                   ) #/ tr
                        ) #/ table
        )),
        
        #      br(),
        
        fluidRow(
          column(3,
                 sliderInput("GBUP", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Upstream:</font>'),
                                        bsButton("qg2", label="", icon=icon("question"), style="info", size="small")
                 ), min = 0, max = 50000, value = 0, ticks = FALSE),
                 bsPopover("qg2", "SNPs in the upstream of the specified genomic region will be used.",
                           trigger = "focus"),
                 sliderInput("GBDOWN", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Downstream:</font>'),
                                          bsButton("qg4", label="", icon=icon("question"), style="info", size="small")
                 ), min = 0, max = 50000, value = 0, ticks = FALSE),
                 bsPopover("qg4", "SNPs in the downstream of the specified genomic region will be used.",
                           trigger = "focus")
          ),
          
          column(4,
                 p(tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Select rice accessions:</font>'),
                      bsButton("qg3", label="", icon=icon("question"), style="info", size="small"))),
                 bsPopover("qg3", "Only the chosen rice accessions will be used. Select from the box on the left to the box on the right. By default, all accessions are chosen.",
                           trigger = "focus"),
                 
                 chooserInput("mychooserB", "Available frobs", "Selected frobs", c(),
                              all.acc.cho, size = 10, multiple = TRUE)
          ),
          
          column(5,
                 tags$div(align = 'left',
                          class = 'multicol', style = "width: 100%",
                          checkboxGroupInput("GB_mut_group", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Mutation types:</font>'),
                                                                bsButton("qg1", label="", icon=icon("question"), style="info", size="small")),
                                             choices = c("Stop_lost","Stop_gained","Start_lost","Start_gained",
                                                         "Intron", "Upstream", "Downstream", "Intergenic",
                                                         "five_prime_UTR","three_prime_UTR",
                                                         "Non_synonymous_start","Non_synonymous_coding",
                                                         "Splice_site_acceptor","Splice_site_donor",
                                                         "Synonymous_stop","Synonymous_coding"
                                             ),
                                             selected = c("Stop_lost","Stop_gained","Start_lost","Start_gained",
                                                          "Intron", "Upstream", "Downstream", "Intergenic",
                                                          "five_prime_UTR","three_prime_UTR",
                                                          "Non_synonymous_start","Non_synonymous_coding",
                                                          "Splice_site_acceptor","Splice_site_donor",
                                                          "Synonymous_stop","Synonymous_coding"
                                             )),
                          bsPopover("qg1", "Only SNPs with selected mutation effects will be used.",
                                    trigger = "focus")
                 )
          )
        ),
        
        downloadButton("downloadsnp.txt", "Download genotype data"),
        downloadButton("downloadsnpInfo.txt", "Download SNPs information"),
        downloadButton("downloadGB.pdf", "Download pdf-file"),
        
        plotly::plotlyOutput("gbrowser", height = '100%', width = '100%'),
        
        br()
        
      ),
      
      # LDheatmap
      tabPanel(
        title = HTML("<strong style='font-size:20px'>LDheatmap</strong>"), icon = icon("project-diagram"),
        
        sidebarPanel(
          textInput("regL", label = tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Genomic region:</font>'),
                                       bsButton("q5", label="", icon=icon("question"), style="info", size="small")),
                    value = ""),
          
          bsPopover("q5", "A genomic region can be determined by chromosome positions or gene locus. For example, chr07:29611303-29669223 or LOC_Os11g35500.",
                    trigger = "focus"),
          
          shinysky::actionButton("submit2", strong("Submit!",
                                                   bsButton("q8", label="", icon=icon("question"), style="info", size="small")
          ), styleclass = "success"),
          shinysky::actionButton("clearLD", strong("Reset"), styleclass = "warning"),
          shinysky::actionButton("LDExam", strong("Load example"), styleclass = "info"),
          conditionalPanel(condition="input.submit2 != '0'", shinysky::busyIndicator(HTML("<div style='color:red;font-size:30px'>Calculation In progress...</div>"), wait = 0)),
          bsPopover("q8", "Whenever the genomic region is updated, please click Submit!",
                    trigger = "focus"),
          
          br(),
          tags$div(HTML('<i class="fa fa-circle" aria-hidden="true"></i> <font size="4" color="black">Plot options:</font>')),
          
          radioButtons("flip", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Flip the figure</font>')), 
                       list("FALSE" = 0, "TRUE" = 1)),
          
          conditionalPanel(
            condition = "input.flip==1",
            checkboxInput("LDshowGene", "Show gene model", FALSE),
            conditionalPanel(
              condition = "input.LDshowGene",
              numericInput("ldY", "Y:", value = 72),
              numericInput("ldW", "W:", value = 72)
            )
          ),
          
          conditionalPanel(
            condition = "input.flip==0",
            radioButtons("showText", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Print LD measurements</font>')), 
                         list("FALSE" =  0, "TRUE" = 1)),
            textInput("ldpos", "Label SNPs:", value = "5, 8")
          ),
          
          radioButtons("ldcol",
                       "Color", list("grey.colors(20)" = 1, "heat.colors(20)" = 2)
          ),
          
          numericInput("ldUp", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Upstream (kb):</font>'),
                                  bsButton("ql4", label="", icon=icon("question"), style="info", size="small")
          ), value = 0),
          bsPopover("ql4", "SNPs in the upstream of the specified genomic region will be used.",
                    trigger = "focus"),
          numericInput("ldDown", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Downstream (kb):</font>'),
                                    bsButton("ql5", label="", icon=icon("question"), style="info", size="small")
          ), value = 0),
          bsPopover("ql5", "SNPs in the downstream of the specified genomic region will be used.",
                    trigger = "focus"),
          
          tags$div(align = 'left',
                   class = 'multicol', style = "width: 100%",
                   checkboxGroupInput("ld_mut_group", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Mutation types:</font>'),
                                                         bsButton("ql1", label="", icon=icon("question"), style="info", size="small")),
                                      choices = c("Stop_lost","Stop_gained","Start_lost","Start_gained",
                                                  "Intron", "Upstream", "Downstream", "Intergenic",
                                                  "five_prime_UTR","three_prime_UTR",
                                                  "Non_synonymous_start","Non_synonymous_coding",
                                                  "Splice_site_acceptor","Splice_site_donor",
                                                  "Synonymous_stop","Synonymous_coding"
                                      ),
                                      selected = c("Stop_lost","Stop_gained","Start_lost","Start_gained",
                                                   "Intron", "Upstream", "Downstream", "Intergenic",
                                                   "five_prime_UTR","three_prime_UTR",
                                                   "Non_synonymous_start","Non_synonymous_coding",
                                                   "Splice_site_acceptor","Splice_site_donor",
                                                   "Synonymous_stop","Synonymous_coding"
                                      )),
                   bsPopover("ql1", "Only SNPs with selected mutation effects will be used.",
                             trigger = "focus")
          ),
          
          p(tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Select rice accessions:</font>'),
               bsButton("ql3", label="", icon=icon("question"), style="info", size="small"))),
          bsPopover("ql3", "Only the chosen rice accessions will be used. Select from the box on the left to the box on the right. By default, all accessions are chosen.",
                    trigger = "focus"),
          
          chooserInput("mychooserLD", "Available frobs", "Selected frobs",
                       c(), all.acc.cho, size = 10, multiple = TRUE
          ),
          
          radioButtons("uploadLD", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">SNP sites to be retained:</font>'),
                                      bsButton("ql2", label="", icon=icon("question"), style="info", size="small")), 
                       c("ALL" = "1", "Upload SNP sites file" = "2"), "1"),
          bsPopover("ql2", "A text file with SNP IDs (one ID per row) could be uploaded to screen the SNPs used in the analysis. Or else, all the SNPs in the specifid genomic region will be used.",
                    trigger = "focus"),
          conditionalPanel(condition="input.uploadLD == '2'",
                           fileInput("LD.snpsite", NULL, multiple = FALSE)),
          
          checkboxInput("ldSize", "Adjust plot size", FALSE),
          conditionalPanel(
            condition = "input.ldSize",
            numericInput("ldHeight", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Plot height:</font>')), value = 550),
            numericInput("ldWidth", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Plot width:</font>')), value = 750)
          )
        ),
        
        mainPanel(
          downloadButton("downloadLD.pdf", "Download pdf-file"),
          downloadButton("downloadLD.svg", "Download svg-file"),
          downloadButton("downloadLD.txt", "Download TXT-file"),
          plotOutput("ldheatmap", height = '100%', width = '100%')
          
        )
      ),
      
      # Haplotype network
      tabPanel(
        title = HTML("<strong style='font-size:20px'>Haplotype</strong>"),icon = icon("cogs"),
        
        sidebarPanel(
          textInput("regH", label = tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Genomic region:</font>'),
                                       bsButton("q4", label="", icon=icon("question"), style="info", size="small")),
                    value = ""),
          
          bsPopover("q4", "A genomic region can be determined by chromosome positions or gene locus. For example, chr07:29611303-29669223 or LOC_Os11g35500.",
                    trigger = "focus"),
          
          shinysky::actionButton("submit3", strong("Submit!",
                                                   bsButton("q9", label="", icon=icon("question"), style="info", size="small")
          ), styleclass = "success"),
          shinysky::actionButton("clearHap", strong("Reset"), styleclass = "warning"),
          shinysky::actionButton("HapExam", strong("Load example"), styleclass = "info"),
          conditionalPanel(condition="input.submit3 != '0'", shinysky::busyIndicator(HTML("<div style='color:red;font-size:30px'>Calculation In progress...</div>"), wait = 0)),
          bsPopover("q9", "Whenever the genomic region or any option is updated, please click Submit!",
                    trigger = "focus"),
          
          tags$div(HTML('<i class="fa fa-circle" aria-hidden="true"></i> <font size="4" color="black">Plot options:</font>')),
          
          selectInput("hapPop", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Populations</font>')), 
                      list("Ind,Jap,TeJ,TrJ,Aus,Or-I,Or-II,Or-III,Other" = 1,
                                                    "Ind,Jap,TeJ,TrJ,Aus,Wild,Other" = 2,
                                                    "Ind,Jap,Aus,Wild,Other" = 3
          )),
          numericInput("hapScale", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Scale ratio:</font>'),
                                      bsButton("qh6", label="", icon=icon("question"), style="info", size="small"))
                       , value = 100),
          bsPopover("qh6", "The ratio of the scale of the links representing the number of steps on the scale of the circles representing the haplotypes. It may be needed to give a value greater than one to avoid overlapping circles.",
                    trigger = "focus"),
          
          selectInput("hapMut", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Show mutations</font>'),
                                   bsButton("qh7", label="", icon=icon("question"), style="info", size="small")), list(
                                     "0" = 0, "1" = 1, "2" = 2, "3" = 3
                                   )),
          bsPopover("qh7", "If 0, nothing is drawn on the links; if 1, the mutations are shown with small segments on the links; if 2, they are shown with small dots; if 3, the number of mutations are printed on the links.",
                    trigger = "focus"),
          
          radioButtons("hapLab", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Show labels</font>'),
                                    bsButton("qh8", label="", icon=icon("question"), style="info", size="small")
          ), list("FALSE" = 1, "TRUE" = 2)),
          bsPopover("qh8", "A logical specifying whether to identify the haplotypes with their labels.",
                    trigger = "focus"),
          checkboxInput("legendPos", "Modify legend position", FALSE),
          conditionalPanel(
            condition = "input.legendPos",
            textInput("hapLenX", "Legend X position:", value = "bottomleft"),
            textInput("hapLenY", "Legend Y position:", value =  "NULL")
          ),
          checkboxInput("hapFreq", "Filter haplotype frequency", FALSE),
          conditionalPanel(
            condition = "input.hapFreq",
            numericInput("hapMin", "Min freq:", value = 50),
            numericInput("hapMax", "Max freq:", value = 2508)
          ),
          
          checkboxInput("hapLink", "Modify links", FALSE),
          conditionalPanel(
            condition = "input.hapLink",
            numericInput("hapLinkWd", "Link width:", value = 1),
            textInput("hapLinkCol", "Link color:", value = "black")
          ),
          
          numericInput("hapUp", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Upstream (kb):</font>'),
                                   bsButton("qh4", label="", icon=icon("question"), style="info", size="small")
          ), value = 0),
          bsPopover("qh4", "SNPs in the upstream of the specified genomic region will be used.",
                    trigger = "focus"),
          numericInput("hapDown", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Downstream (kb):</font>'),
                                     bsButton("qh5", label="", icon=icon("question"), style="info", size="small")
          ), value = 0),
          bsPopover("qh5", "SNPs in the downstream of the specified genomic region will be used.",
                    trigger = "focus"),
          
          tags$div(align = 'left',
                   class = 'multicol', style = "width: 100%",
                   checkboxGroupInput("hap_mut_group", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Mutation types:</font>'),
                                                          bsButton("qh1", label="", icon=icon("question"), style="info", size="small")),
                                      choices = c("Stop_lost","Stop_gained","Start_lost","Start_gained",
                                                  "Intron", "Upstream", "Downstream", "Intergenic",
                                                  "five_prime_UTR","three_prime_UTR",
                                                  "Non_synonymous_start","Non_synonymous_coding",
                                                  "Splice_site_acceptor","Splice_site_donor",
                                                  "Synonymous_stop","Synonymous_coding"
                                      ),
                                      selected = c("Stop_lost","Stop_gained","Start_lost","Start_gained",
                                                   "Intron", "Upstream", "Downstream", "Intergenic",
                                                   "five_prime_UTR","three_prime_UTR",
                                                   "Non_synonymous_start","Non_synonymous_coding",
                                                   "Splice_site_acceptor","Splice_site_donor",
                                                   "Synonymous_stop","Synonymous_coding"
                                      )),
                   bsPopover("qh1", "Only SNPs with selected mutation effects will be used.",
                             trigger = "focus")
          ),
          
          radioButtons("uploadHAP", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">SNP sites to be retained:</font>'),
                                       bsButton("qh2", label="", icon=icon("question"), style="info", size="small")), 
                       c("ALL" = "1", "Upload SNP sites file" = "2"), "1"),
          bsPopover("qh2", "A text file with SNP IDs (one ID per row) could be uploaded to screen the SNPs used in the analysis. Or else, all the SNPs in the specifid genomic region will be used.",
                    trigger = "focus"),
          conditionalPanel(condition="input.uploadHAP == '2'",
                           fileInput("HAP.snpsite", NULL, multiple = FALSE))
        ),
        
        mainPanel(
          downloadButton("downloadHap.pdf", "Download pdf-file"),
          downloadButton("downloadHap.svg", "Download svg-file"),
          downloadButton("downloadHap.nex", "Download NEXUS-file"),
          plotOutput("haplotype", height = '100%', width = '100%'),
          
          br(),
          
          downloadButton("downloadHapSta.pdf", "Download pdf-file"),
          downloadButton("downloadHapSta.svg", "Download svg-file"),
          plotly::plotlyOutput("hapGeo")
          
        )
      ),
      
      # Nucleotide diversity
      tabPanel(
        title = HTML("<strong style='font-size:20px'>Diversity</strong>"), icon = icon("chart-area"),
        
        sidebarPanel(
          textInput("regD", label = tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Genomic region:</font>'),
                                       bsButton("q3", label="", icon=icon("question"), style="info", size="small")),
                    value = ""),
          
          bsPopover("q3", "A genomic region can be determined by chromosome positions or gene locus. For example, chr07:29611303-29669223 or LOC_Os11g35500.",
                    trigger = "focus"),
          
          shinysky::actionButton("submit4", strong("Submit!",
                                                   bsButton("q10", label="", icon=icon("question"), style="info", size="small")
          ), styleclass = "success"),
          shinysky::actionButton("clearDiv", strong("Reset"), styleclass = "warning"),
          shinysky::actionButton("DivExam", strong("Load example"), styleclass = "info"),
          conditionalPanel(condition="input.submit4 != '0'", shinysky::busyIndicator(HTML("<div style='color:red;font-size:30px'>Calculation In progress...</div>"), wait = 0)),
          bsPopover("q10", "Whenever the genomic region or any option is updated, please click Submit!!",
                    trigger = "focus"),
          
          tags$div(HTML('<i class="fa fa-circle" aria-hidden="true"></i> <font size="4" color="black">Plot options:</font>')),
          
          numericInput("snpnumD", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Number of SNPs in each window:</font>'),
                                     bsButton("qd6", label="", icon=icon("question"), style="info", size="small")
          ), value = 10, min = 5, max = 20),
          bsPopover("qd6", "A specified genomic region would be split into non-overlapping window so that each window contains specified number of SNPs. The nucleotide diversity of all rice accessions belong to the specified ecotypes in each window would be calculated.",
                    trigger = "focus"),
          
          tags$div(align = 'left',
                   class = 'multicol', style = "width: 100%",
                   checkboxGroupInput("div_acc_group", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Ecotypes to calculate diversity:</font>')),
                                      choices = c("Wild", "Cultivar", "Aus", "Indica", "IndicaI", "IndicaII",
                                                  "Japonica", "TrJ", "TeJ", "Or-I", "Or-II", "Or-III"),
                                      selected = c("Wild", "Cultivar")) 
          ),
          
          selectInput("nuc_numerator", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Numerator ecotype:</font>'),
                                          bsButton("qd7", label="", icon=icon("question"), style="info", size="small")
          ), choices = 
            c("Wild", "Cultivar", "Aus", "Indica", "IndicaI", "IndicaII",
              "Japonica", "TrJ", "TeJ", "Or-I", "Or-II", "Or-III")),
          bsPopover("qd7", "The nucleotide diversity of rice accessions belong to the Numerator ecotype would be divided by the nucleotide diversity of rice accessions belong to the Denominator ecotype for comparison.",
                    trigger = "focus"),
          selectInput("nuc_denominator", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Denominator ecotype:</font>')), choices = 
                        c("Cultivar", "Wild", "Aus", "Indica", "IndicaI", "IndicaII",
                          "Japonica", "TrJ", "TeJ", "Or-I", "Or-II", "Or-III")),
          
          tags$div(align = 'left',
                   class = 'multicol', style = "width: 100%",
                   checkboxGroupInput("div_mut_group", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Mutation types:</font>'),
                                                          bsButton("qd1", label="", icon=icon("question"), style="info", size="small")),
                                      choices = c("Stop_lost","Stop_gained","Start_lost","Start_gained",
                                                  "Intron", "Upstream", "Downstream", "Intergenic",
                                                  "five_prime_UTR","three_prime_UTR",
                                                  "Non_synonymous_start","Non_synonymous_coding",
                                                  "Splice_site_acceptor","Splice_site_donor",
                                                  "Synonymous_stop","Synonymous_coding"
                                      ),
                                      selected = c("Stop_lost","Stop_gained","Start_lost","Start_gained",
                                                   "Intron", "Upstream", "Downstream", "Intergenic",
                                                   "five_prime_UTR","three_prime_UTR",
                                                   "Non_synonymous_start","Non_synonymous_coding",
                                                   "Splice_site_acceptor","Splice_site_donor",
                                                   "Synonymous_stop","Synonymous_coding"
                                      )),
                   bsPopover("qd1", "Only SNPs with selected mutation effects will be used.",
                             trigger = "focus")
          ),
          
          
          numericInput("divUp", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Upstream (kb):</font>'),
                                   bsButton("qd4", label="", icon=icon("question"), style="info", size="small")
          ), value = 20),
          bsPopover("qd4", "SNPs in the upstream of the specified genomic region will be used.",
                    trigger = "focus"),
          numericInput("divDown", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Downstream (kb):</font>'),
                                     bsButton("qd5", label="", icon=icon("question"), style="info", size="small")
          ), value = 20),
          bsPopover("qd5", "SNPs in the downstream of the specified genomic region will be used.",
                    trigger = "focus"),
          
          radioButtons("uploadDIV", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">SNP sites to be retained:</font>'),
                                       bsButton("qd2", label="", icon=icon("question"), style="info", size="small")), 
                       c("ALL" = "1", "Upload SNP sites file" = "2"), "1"),
          bsPopover("qd2", "A text file with SNP IDs (one ID per row) could be uploaded to screen the SNPs used in the analysis. Or else, all the SNPs in the specifid genomic region will be used.",
                    trigger = "focus"),
          conditionalPanel(condition="input.uploadDIV == '2'",
                           fileInput("DIV.snpsite", NULL, multiple = FALSE)),
          
          checkboxInput("divSize", "Adjust plot size", FALSE),
          conditionalPanel(
            condition = "input.divSize",
            numericInput("divHeight", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Plot height:</font>')), value = 550),
            numericInput("divWidth", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Plot width:</font>')), value = 750)
          )
        ),
        
        mainPanel(
          fluidRow(
            column(3, uiOutput("downloadDiv01")),
            column(3, uiOutput("downloadDiv02")),
            column(3, uiOutput("downloadDiv03"))
          ),
          
          plotOutput("diversity", height = '100%', width = '100%')
        )
        
      ),
      
      # Phylogenetic tree
      tabPanel(
        title = HTML("<strong style='font-size:20px'>Phylogenetic</strong>"),icon = icon("bezier-curve"),
        
        sidebarPanel(
          textInput("regP", label = tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Genomic region:</font>'),
                                       bsButton("q2", label="", icon=icon("question"), style="info", size="small")),
                    value = ""),
          
          bsPopover("q2", "A genomic region can be determined by chromosome positions or gene locus. For example, chr07:29611303-29669223 or LOC_Os11g35500.",
                    trigger = "focus"),
          
          shinysky::actionButton("submit5", strong("Submit!",
                                                   bsButton("q11", label="", icon=icon("question"), style="info", size="small")
          ), styleclass = "success"),
          shinysky::actionButton("clearPhy", strong("Reset"), styleclass = "warning"),
          shinysky::actionButton("PhyExam", strong("Load example"), styleclass = "info"),
          conditionalPanel(condition="input.submit5 != '0'", shinysky::busyIndicator(HTML("<p style='color:red;font-size:30px;'>Calculation In progress...</p>"), wait = 0)),
          bsPopover("q11", "Whenever the genomic region or any option is updated, please click Submit!",
                    trigger = "focus"),
          
          tags$div(HTML('<i class="fa fa-circle" aria-hidden="true"></i> <font size="4" color="black">Plot options:</font>')),
          numericInput("phyUp", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Upstream (kb):</font>'),
                                   bsButton("qp4", label="", icon=icon("question"), style="info", size="small")
          ), value = 20),
          bsPopover("qp4", "SNPs in the upstream of the specified genomic region will be used.",
                    trigger = "focus"),
          numericInput("phyDown", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Downstream (kb):</font>'),
                                     bsButton("qp5", label="", icon=icon("question"), style="info", size="small")
          ), value = 20),
          bsPopover("qp5", "SNPs in the downstream of the specified genomic region will be used.",
                    trigger = "focus"),
          
          tags$div(align = 'left',
                   class = 'multicol', style = "width: 100%",
                   checkboxGroupInput("phy_mut_group", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Mutation types:</font>'),
                                                          bsButton("qp1", label="", icon=icon("question"), style="info", size="small")),
                                      choices = c("Stop_lost","Stop_gained","Start_lost","Start_gained",
                                                  "Intron", "Upstream", "Downstream", "Intergenic",
                                                  "five_prime_UTR","three_prime_UTR",
                                                  "Non_synonymous_start","Non_synonymous_coding",
                                                  "Splice_site_acceptor","Splice_site_donor",
                                                  "Synonymous_stop","Synonymous_coding"
                                      ),
                                      selected = c("Stop_lost","Stop_gained","Start_lost","Start_gained",
                                                   "Intron", "Upstream", "Downstream", "Intergenic",
                                                   "five_prime_UTR","three_prime_UTR",
                                                   "Non_synonymous_start","Non_synonymous_coding",
                                                   "Splice_site_acceptor","Splice_site_donor",
                                                   "Synonymous_stop","Synonymous_coding"
                                      )),
                   bsPopover("qp1", "Only SNPs with selected mutation effects will be used.",
                             trigger = "focus")
          ),
          
          p(tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Select rice accessions:</font>'),
               bsButton("qp3", label="", icon=icon("question"), style="info", size="small"))),
          bsPopover("qp3", "Only the chosen rice accessions will be used. Select from the box on the left to the box on the right. By default, all accessions are chosen.",
                    trigger = "focus"),
          
          chooserInput("mychooserPhy", "Available frobs", "Selected frobs",
                       c(), all.acc.cho, size = 10, multiple = TRUE
          ),
          
          radioButtons("uploadPHY", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">SNP sites to be retained:</font>'),
                                       bsButton("qp2", label="", icon=icon("question"), style="info", size="small")), 
                       c("ALL" = "1", "Upload SNP sites file" = "2"), "1"),
          bsPopover("qp2", "A text file with SNP IDs (one ID per row) could be uploaded to screen the SNPs used in the analysis. Or else, all the SNPs in the specifid genomic region will be used.",
                    trigger = "focus"),
          conditionalPanel(condition="input.uploadPHY == '2'",
                           fileInput("PHY.snpsite", NULL, multiple = FALSE)),
          
          checkboxInput("phySize", "Adjust plot size", FALSE),
          conditionalPanel(
            condition = "input.phySize",
            numericInput("phyHeight", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Plot height:</font>')), value = 700),
            numericInput("phyWidth", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Plot width:</font>')), value = 750)
          )
        ),
        
        mainPanel(
          fluidRow(
            column(5, uiOutput("downloadPhy01")),
            column(5, uiOutput("downloadPhy02"))
          ),
          
          plotOutput("phylo", height = '100%', width = '100%')
          
        )
        
      ),
      
      # Allele frequency
      tabPanel(
        title = HTML("<strong style='font-size:20px'>AlleleFreq</strong>"), icon = icon("chart-pie"),
        
        sidebarPanel(
          textAreaInput("af_snp_site", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">SNP sites to calculate allele frequency:</font>'),
                                          bsButton("qaf3", label="", icon=icon("question"), style="info", size="small")
          ), 
          width="100%", resize="vertical", height="150px", 
          placeholder = "One SNP site in one row", 
          value = ""
          ),
          bsPopover("qaf3", "Each SNP site should be a 10-digits integer and the first two digits represent the chromosome ID while the rest eight digits represent the genomic position of each SNP site. Each SNP site should take only one row!",
                    trigger = "focus"),
          
          tags$div(align = 'left',
                   class = 'multicol', style = "width: 100%",
                   checkboxGroupInput("af_acc_group", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Ecotypes to calculate allele frequency:</font>')),
                                      choices = c("Wild", "Cultivar", "Aus", "Indica", "IndicaI", "IndicaII",
                                                  "Japonica", "TrJ", "TeJ", "Or-I", "Or-II", "Or-III"),
                                      selected = c("Aus", "Indica", "TeJ", "TrJ", "Wild")) 
          ),
          textInput("afCol", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Allele colors:</font>'),
                                bsButton("qaf2", label="", icon=icon("question"), style="info", size="small")
          ), value = "cornflowerblue, forestgreen"),
          bsPopover("qaf2", "Colors for the major and minor allele in the pie chart respectively!",
                    trigger = "focus"),
          
          numericInput("afHeight", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Plot height:</font>')), value = 550),
          numericInput("afWidth", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Plot width:</font>')), value = 700),
          
          shinysky::actionButton("submitaf1", strong("Submit!",
                                                     bsButton("qaf1", label="", icon=icon("question"), style="info", size="small")
          ), styleclass = "success"),
          shinysky::actionButton("clearAf", strong("Reset"), styleclass = "warning"),
          shinysky::actionButton("AfExam", strong("Load example"), styleclass = "info"),
          conditionalPanel(condition="input.submitaf1 != '0'", shinysky::busyIndicator(HTML("<p style='color:red;font-size:30px;'>Calculation In progress...</p>"), wait = 0)),
          bsPopover("qaf1", "Whenever the SNP sites or any option is updated, please click Submit!",
                    trigger = "focus")
        ),
        
        mainPanel(
          fluidRow(
            column(4, uiOutput("downloadAfq01")),
            column(4, uiOutput("downloadAfq02")),
            column(4, uiOutput("downloadAfq03"))
          ),
          
          br(),
          
          plotOutput("alleleFreq", height = "550px", width = "700px")
        )
      ),
      
      # Accession
      tabPanel(
        title = HTML("<strong style='font-size:20px'>Accession</strong>"), icon = icon("list"),
        
        sidebarPanel(
          downloadButton("acc.info.txt", "Download information of all accessions"),
          
          p(tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Select rice accessions:</font>'),
               bsButton("qa1", label="", icon=icon("question"), style="info", size="small"))),
          bsPopover("qa1", "Only the chosen rice accessions will be used. Select from the box on the left to the box on the right. By default, all accessions are chosen.",
                    trigger = "focus"),
          
          chooserInput("mychooserA", "Available frobs", "Selected frobs",
                       c(), all.acc.cho, size = 10, multiple = TRUE
          ),
          
          br(),
          downloadButton("sel.acc.info.txt", "Download information of selected accessions")
        ),
        
        mainPanel(
          downloadButton("downloadAccDis.pdf", "Download pdf-file"),
          downloadButton("downloadAccDis.svg", "Download svg-file"),
          shinycssloaders::withSpinner(plotly::plotlyOutput("accDis")),
          
          h4("Information of selected rice accessions"),
          shinycssloaders::withSpinner(dataTableOutput("mytable1"))
        )
      ),
      
      # Bulk download of data
      tabPanel(
        title = HTML("<strong style='font-size:20px'>Download</strong>"), icon = icon("download"),
        
        sidebarPanel(
          textInput("regBB", label = tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Genomic region:</font>'),
                                        bsButton("q1", label="", icon=icon("question"), style="info", size="small")),
                    value = ""),
          
          bsPopover("q1", "A genomic region can be determined by chromosome positions or gene locus. For example, chr07:29611303-29669223 or LOC_Os11g35500.",
                    trigger = "focus"),
          
          shinysky::actionButton("submit6", strong("Submit!",
                                                   bsButton("q12", label="", icon=icon("question"), style="info", size="small")
          ), width = "60%", styleclass = "success"),
          shinysky::actionButton("clearBB", strong("Reset"), styleclass = "warning"),
          shinysky::actionButton("BBExam", strong("Load example"), styleclass = "info"),
          conditionalPanel(condition="input.submit6 != '0'", shinysky::busyIndicator(HTML("<p style='color:red;font-size:30px;'>Calculation In progress...</p>"), wait = 0)),
          bsPopover("q12", "Whenever the genomic region or any option is updated, please click Submit!",
                    trigger = "focus"),
          
          br(),
          
          p(tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Select rice accessions:</font>'),
               bsButton("qdl2", label="", icon=icon("question"), style="info", size="small"))),
          bsPopover("qdl2", "Only the chosen rice accessions will be used. Select from the box on the left to the box on the right. By default, all accessions are chosen.",
                    trigger = "focus"),
          
          chooserInput("mychooserD", "Available frobs", "Selected frobs",
                       c(), all.acc.cho, size = 10, multiple = TRUE
          ),
          
          tags$div(align = 'left',
                   class = 'multicol', style = "width: 100%",
                   checkboxGroupInput("down_mut_group", tags$div(HTML('<i class="fa fa-play" aria-hidden="true"></i> <font size="4" color="red">Mutation types:</font>'),
                                                           bsButton("qdl1", label="", icon=icon("question"), style="info", size="small")),
                                      choices = c("Stop_lost","Stop_gained","Start_lost","Start_gained",
                                                  "Intron", "Upstream", "Downstream", "Intergenic",
                                                  "five_prime_UTR","three_prime_UTR",
                                                  "Non_synonymous_start","Non_synonymous_coding",
                                                  "Splice_site_acceptor","Splice_site_donor",
                                                  "Synonymous_stop","Synonymous_coding"
                                      ),
                                      selected = c("Stop_lost","Stop_gained","Start_lost","Start_gained",
                                                   "Intron", "Upstream", "Downstream", "Intergenic",
                                                   "five_prime_UTR","three_prime_UTR",
                                                   "Non_synonymous_start","Non_synonymous_coding",
                                                   "Splice_site_acceptor","Splice_site_donor",
                                                   "Synonymous_stop","Synonymous_coding"
                                      )),
                   bsPopover("qdl1", "Only SNPs with selected mutation effects will be used.",
                             trigger = "focus")
          )
        ),
        
        mainPanel(
          downloadButton("bulkdownloadsnpInfo.txt", "Download SNPs information"),
          downloadButton("bulkdownloadsnp.txt", "Download genotype data"),
          downloadButton("bulkdownloadgene.txt", "Download gene annotation"),
          dataTableOutput("mytable2")
          
        )
      ),
      
      ## Help
      tabPanel(
               title = HTML("<strong style='font-size:20px'>Help</strong>"), icon = icon("book"),
               includeMarkdown("README.md")
      ),
      
      footer = footerTagList
    )
  )
  
)

