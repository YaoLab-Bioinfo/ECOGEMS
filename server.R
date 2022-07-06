
# options(shiny.maxRequestSize = 200*1024^2)

shinyServer(function(input, output, session) {
  
  # GBrowser
  observe({
    if (input$submit1>0) {
      isolate({
        myPos <- anaReg(input$regB)
        
        if (validReg(myPos)) {
          if (!is.null(myPos)) {
            snp.info <- snpInfo(chr=myPos$chr, start=myPos$start - input$GBUP, end=myPos$end + input$GBDOWN, 
                                accession = input$mychooserB$selected, mutType = input$GB_mut_group)
          } else {
            snp.info <- NULL
          }
          
          if (is.null(snp.info) || nrow(snp.info[[1]][[1]]) < 1) {
            shinyWidgets::sendSweetAlert(
              session = session,
              title = "Error input!", type = "error",
              text = "No SNPs are detected in the specified genomic region or the specified genomic region is too large!"
            )
          } else {
            GBplot <<- NULL
            output$gbrowser <- plotly::renderPlotly({
              GBplot <<- GBrowser(chr=myPos$chr, start=myPos$start - input$GBUP, 
                                  end=myPos$end + input$GBDOWN,
                                  accession = input$mychooserB$selected,
                                  mutType = input$GB_mut_group)
              GBplot[[2]]
            })
            
            ## Download PDF file of GBrowser
            output$downloadGB.pdf <- downloadHandler(
              filename <- function() { paste('GBrowser.pdf') },
              content <- function(file) {
                pdf(file, width = 900/72, height = 300/72)
                grid::grid.draw(GBplot[[1]])
                dev.off()
              }, contentType = 'application/pdf')
            
            # Download genotypes of seleceted SNPs
            output$downloadsnp.txt <- downloadHandler(
              filename = function() { "snp.geno.txt" },
              content = function(file) {
                write.table(snp.info[[1]][[1]], file, sep="\t", quote=F)
              })
            
            # Download information of SNPs
            output$downloadsnpInfo.txt <- downloadHandler(
              filename = function() { "snp.info.txt" },
              content = function(file) {
                write.table(snp.info[[2]], file, sep="\t", quote=F, row.names=F)
              })
          }
        } else {
          shinyWidgets::sendSweetAlert(
            session = session,
            title = "Error input!", type = "error",
            text = "Please input genomic region or gene model in appropriate format!"
          )
        }
      })
    } else {
      NULL
    }
  })
  
  observe({
    if (input$clearGB>0) {
      isolate({
        updateTextInput(session, "regB", value="")
      })
    } else {NULL}
  })
  
  observe({
    if (input$GBExam >0) {
      isolate({
        updateTextInput(session, "regB", value="chr07:29611303-29669223")
      })
    } else {NULL}
  })
    
  # LDheatmap
  observe({
    if (input$submit2>0) {
      isolate({
        ld.height <<- input$ldHeight
        ld.width <<- input$ldWidth
        myPos <- anaReg(input$regL)
        
        if (validReg(myPos)) {
          if (!is.null(myPos)) {
            snp.reg <- fetchSnp(chr=myPos$chr, start=myPos$start - input$ldUp * 1000, 
                                end=myPos$end + input$ldDown * 1000, accession = input$mychooserLD$selected,
                                mutType = input$ld_mut_group)[[1]]
          } else {
            snp.reg <- NULL
          }
          
          if (is.null(snp.reg) || nrow(snp.reg) < 5) {
            shinyWidgets::sendSweetAlert(
              session = session,
              title = "Error input!", type = "error",
              text = "Too few SNPs are detected in the specified genomic region or the specified genomic region is too large!"
            )
          } else {
            snp.pos <- as.numeric(unlist(strsplit(input$ldpos, split=",")))
            
            if (input$uploadLD == 1) {
              ld.snp.site <- NULL
            } else {
              if (!is.null(input$LD.snpsite)) {
                ld.snp.site <- readLines(input$LD.snpsite$datapath)
              } else {
                ld.snp.site <- NULL
              }
            }
            
            output$ldheatmap <- renderPlot({
              if (input$flip == "0") {
                ld.heatmap(chr=myPos$chr, start=myPos$start - input$ldUp * 1000, end=myPos$end + input$ldDown * 1000, text=c(FALSE, TRUE)[as.numeric(input$showText)+1],
                           snp.pos=snp.pos, gene=FALSE, flip=c(FALSE, TRUE)[as.numeric(input$flip)+1],
                           col=list(grey.colors(20), heat.colors(20))[[as.numeric(input$ldcol)]],
                           mutType = input$ld_mut_group, accession = input$mychooserLD$selected, 
                           snpSites = ld.snp.site)
              } else if (input$flip == "1") {
                if (input$LDshowGene) {
                  ld.heatmap(chr=myPos$chr, start=myPos$start - input$ldUp * 1000, end=myPos$end + input$ldDown * 1000, text=FALSE,
                             snp.pos=snp.pos, ld.y=input$ldY/100, ld.w=input$ldW/100, gene=TRUE, 
                             col=list(grey.colors(20), heat.colors(20))[[as.numeric(input$ldcol)]],
                             mutType = input$ld_mut_group, accession = input$mychooserLD$selected, 
                             snpSites = ld.snp.site)
                } else {
                  ld.heatmap(chr=myPos$chr, start=myPos$start - input$ldUp * 1000, end=myPos$end + input$ldDown * 1000, text=FALSE,
                             gene=FALSE, flip=c(FALSE, TRUE)[as.numeric(input$flip)+1],
                             col=list(grey.colors(20), heat.colors(20))[[as.numeric(input$ldcol)]],
                             mutType = input$ld_mut_group, accession = input$mychooserLD$selected, 
                             snpSites = ld.snp.site)
                }
              }
              
            }, height = ld.height, width = ld.width)
          }
        } else {
          shinyWidgets::sendSweetAlert(
            session = session,
            title = "Error input!", type = "error",
            text = "Please input genomic region or gene model in appropriate format!"
          )
        }
      })
    } else {
      NULL
    }
  })
  
  observe({
    if (input$clearLD>0) {
      isolate({
        updateTextInput(session, "regL", value="")
      })
    } else {NULL}
  })
  
  observe({
    if (input$LDExam >0) {
      isolate({
        updateTextInput(session, "regL", value="LOC_Os11g35500")
      })
    } else {NULL}
  })

	## Download PDF file of LDheatmap
	output$downloadLD.pdf <- downloadHandler(
	  filename <- function() { paste('LDheatmap.pdf') },
	  content <- function(file) {
	      myPos <- anaReg(input$regL)
	    
	    snp.reg <- fetchSnp(chr=myPos$chr, start=myPos$start - input$ldUp * 1000, 
	                        end=myPos$end + input$ldDown * 1000, accession = input$mychooserLD$selected,
	                        mutType = input$ld_mut_group)[[1]]
	    if (nrow(snp.reg) < 5) {
	      js_string <- 'alert("Too few SNPs in specified genomic region!");'
	      session$sendCustomMessage(type='jsCode', list(value = js_string))
	    } else {
	      pdf(file, width = input$ldWidth/72, height = input$ldHeight/72, onefile = FALSE)
	      snp.pos <- as.numeric(unlist(strsplit(input$ldpos, split=",")))
	      
	      if (input$uploadLD == 1) {
	        ld.snp.site <- NULL
	      } else {
	        if (!is.null(input$LD.snpsite)) {
	          ld.snp.site <- readLines(input$LD.snpsite$datapath)
	        } else {
	          ld.snp.site <- NULL
	        }
	      }
	      
	      if (input$flip == "0") {
	        ld.heatmap(chr=myPos$chr, start=myPos$start - input$ldUp * 1000, end=myPos$end + input$ldDown * 1000, text=c(FALSE, TRUE)[as.numeric(input$showText)+1],
	                   snp.pos=snp.pos, gene=FALSE, flip=c(FALSE, TRUE)[as.numeric(input$flip)+1],
	                   col=list(grey.colors(20), heat.colors(20))[[as.numeric(input$ldcol)]],
	                   mutType = input$ld_mut_group, accession = input$mychooserLD$selected, 
	                   snpSites = ld.snp.site)
	      } else if (input$flip == "1") {
	        if (input$LDshowGene) {
	          ld.heatmap(chr=myPos$chr, start=myPos$start - input$ldUp * 1000, end=myPos$end + input$ldDown * 1000, text=FALSE,
	                     snp.pos=snp.pos, ld.y=input$ldY/100, ld.w=input$ldW/100, gene=TRUE, 
	                     col=list(grey.colors(20), heat.colors(20))[[as.numeric(input$ldcol)]],
	                     mutType = input$ld_mut_group, accession = input$mychooserLD$selected, 
	                     snpSites = ld.snp.site)
	        } else {
	          ld.heatmap(chr=myPos$chr, start=myPos$start - input$ldUp * 1000, end=myPos$end + input$ldDown * 1000, text=FALSE,
	                     gene=FALSE, flip=c(FALSE, TRUE)[as.numeric(input$flip)+1],
	                     col=list(grey.colors(20), heat.colors(20))[[as.numeric(input$ldcol)]],
	                     mutType = input$ld_mut_group, accession = input$mychooserLD$selected, 
	                     snpSites = ld.snp.site)
	        }
	      }
	      dev.off()
	    }
	    
	  }, contentType = 'application/pdf')
	
	## Download SVG file of LDheatmap
	output$downloadLD.svg <- downloadHandler(
	  filename <- function() { paste('LDheatmap.svg') },
	  content <- function(file) {
	      myPos <- anaReg(input$regL)
	    
	    snp.reg <- fetchSnp(chr=myPos$chr, start=myPos$start - input$ldUp * 1000, 
	                        end=myPos$end + input$ldDown * 1000, accession = input$mychooserLD$selected,
	                        mutType = input$ld_mut_group)[[1]]
	    if (nrow(snp.reg) < 5) {
	      js_string <- 'alert("Too few SNPs in specified genomic region!");'
	      session$sendCustomMessage(type='jsCode', list(value = js_string))
	    } else {
	      svg(file, width = input$ldWidth/72, height = input$ldHeight/72)
	      snp.pos <- as.numeric(unlist(strsplit(input$ldpos, split=",")))
	      
	      if (input$uploadLD == 1) {
	        ld.snp.site <- NULL
	      } else {
	        if (!is.null(input$LD.snpsite)) {
	          ld.snp.site <- readLines(input$LD.snpsite$datapath)
	        } else {
	          ld.snp.site <- NULL
	        }
	      }
	      
	      if (input$flip == "0") {
	        ld.heatmap(chr=myPos$chr, start=myPos$start - input$ldUp * 1000, end=myPos$end + input$ldDown * 1000, text=c(FALSE, TRUE)[as.numeric(input$showText)+1],
	                   snp.pos=snp.pos, gene=FALSE, flip=c(FALSE, TRUE)[as.numeric(input$flip)+1],
	                   col=list(grey.colors(20), heat.colors(20))[[as.numeric(input$ldcol)]],
	                   mutType = input$ld_mut_group, accession = input$mychooserLD$selected, 
	                   snpSites = ld.snp.site)
	      } else if (input$flip == "1") {
	        if (input$LDshowGene) {
	          ld.heatmap(chr=myPos$chr, start=myPos$start - input$ldUp * 1000, end=myPos$end + input$ldDown * 1000, text=FALSE,
	                     snp.pos=snp.pos, ld.y=input$ldY/100, ld.w=input$ldW/100, gene=TRUE, 
	                     col=list(grey.colors(20), heat.colors(20))[[as.numeric(input$ldcol)]],
	                     mutType = input$ld_mut_group, accession = input$mychooserLD$selected, 
	                     snpSites = ld.snp.site)
	        } else {
	          ld.heatmap(chr=myPos$chr, start=myPos$start - input$ldUp * 1000, end=myPos$end + input$ldDown * 1000, text=FALSE,
	                     gene=FALSE, flip=c(FALSE, TRUE)[as.numeric(input$flip)+1],
	                     col=list(grey.colors(20), heat.colors(20))[[as.numeric(input$ldcol)]],
	                     mutType = input$ld_mut_group, accession = input$mychooserLD$selected, 
	                     snpSites = ld.snp.site)
	        }
	      }
	      dev.off()
	    }
	    
	  }, contentType = 'image/svg')
	
	# Haplotype
	observe({
	  if (input$submit3>0) {
	    isolate({
	      if (input$hapLenY == "NULL") {
	        hap.leg.x <- input$hapLenX
	        hap.leg.y <- NULL
	      } else {
	        hap.leg.x <- as.numeric(input$hapLenX)
	        hap.leg.y <- as.numeric(input$hapLenY)
	      }
	      
	      hap.lab <- c(FALSE, TRUE)[as.numeric(input$hapLab)]
	      hap.mut <- as.numeric(input$hapMut)
	      hap.pop <- as.numeric(input$hapPop)
	      hap.scale <- input$hapScale
	      hap.min <- input$hapMin
	      hap.max <- input$hapMax
	      link.wd <- input$hapLinkWd
	      link.col <- gsub("0x","#", input$hapLinkCol)
	      
	      hap.up <- input$hapUp * 1000
	      hap.down <- input$hapDown * 1000
	      
	      hap.mut.grp <- input$hap_mut_group
	      
	      myPos <- anaReg(input$regH)
	      
	      if (validReg(myPos)) {
	        if (!is.null(myPos)) {
	          snp.reg <- fetchSnp(chr=myPos$chr, start=myPos$start - hap.up, end=myPos$end + hap.down, mutType=hap.mut.grp)[[1]]
	        } else {
	          snp.reg <- NULL
	        }
	        
	        if (is.null(snp.reg) || nrow(snp.reg) < 5) {
	          shinyWidgets::sendSweetAlert(
	            session = session,
	            title = "Error input!", type = "error",
	            text = "Too few SNPs are detected in the specified genomic region or the specified genomic region is too large!"
	          )
	        } else {
	          if (input$uploadHAP == 1) {
	            hap.snp.site <- NULL
	          } else {
	            if (!is.null(input$HAP.snpsite)) {
	              hap.snp.site <- readLines(input$HAP.snpsite$datapath)
	            } else {
	              hap.snp.site <- NULL
	            }
	          }
	          
	          output$haplotype <- renderPlot({
	            hapNet(data = snp.reg, legend.x=hap.leg.x, legend.y=hap.leg.y, pop.list=hap.pop,
	                   labels=hap.lab, show.mutation=hap.mut, scale.ratio=hap.scale,
	                   min.freq=hap.min, max.freq=hap.max, lwd=link.wd, col.link=link.col, snpSites = hap.snp.site)
	          }, height = 550, width = 750)
	        }
	      } else {
	        shinyWidgets::sendSweetAlert(
	          session = session,
	          title = "Error input!", type = "error",
	          text = "Please input genomic region or gene model in appropriate format!"
	        )
	      }
	    })
	  } else {
	    NULL
	  }
	})
	
	observe({
	  if (input$clearHap>0) {
	    isolate({
	      updateTextInput(session, "regH", value="")
	    })
	  } else {NULL}
	})
	
	observe({
	  if (input$HapExam >0) {
	    isolate({
	      updateTextInput(session, "regH", value="LOC_Os10g32600")
	    })
	  } else {NULL}
	})
	
	## Download PDF file of haplotype
	output$downloadHap.pdf <- downloadHandler(
	  filename <- function() { paste('haplotype.pdf') },
	  content <- function(file) {
	    withProgress(message='Calculation in progress...',value = 0, detail = 'This may take a while...', {
	      pdf(file, width = 750/72, height = 550/72)
	    if (input$hapLenY == "NULL") {
	      hap.leg.x <- input$hapLenX
	      hap.leg.y <- NULL
	    } else {
	      hap.leg.x <- as.numeric(input$hapLenX)
	      hap.leg.y <- as.numeric(input$hapLenY)
	    }
	    
	    myPos <- anaReg(input$regH)
	    snp.reg <- fetchSnp(chr=myPos$chr, start=myPos$start - input$hapUp * 1000, end=myPos$end + input$hapDown * 1000,
	                        mutType=input$hap_mut_group)[[1]]
	    
	    if (input$uploadHAP == 1) {
	      hap.snp.site <- NULL
	    } else {
	      if (!is.null(input$HAP.snpsite)) {
	        hap.snp.site <- readLines(input$HAP.snpsite$datapath)
	      } else {
	        hap.snp.site <- NULL
	      }
	    }
	    
	    hapNet(data = snp.reg, legend.x=hap.leg.x, legend.y=hap.leg.y, 
	           labels=c(FALSE, TRUE)[as.numeric(input$hapLab)],
	           show.mutation=as.numeric(input$hapMut), scale.ratio=input$hapScale,
	           min.freq=input$hapMin, max.freq=input$hapMax, lwd=input$hapLinkWd,
	           col.link=gsub("0x","#", input$hapLinkCol), pop.list=as.numeric(input$hapPop), snpSites = hap.snp.site)
	    dev.off()
	      
	    })
	    
	  }, contentType = 'application/pdf')
	
	## Download SVG file of haplotype
	output$downloadHap.svg <- downloadHandler(
	  filename <- function() { paste('haplotype.svg') },
	  content <- function(file) {
	    withProgress(message='Calculation in progress...',value = 0, detail = 'This may take a while...', {
	      svg(file, width = 750/72, height = 550/72)
	    if (input$hapLenY == "NULL") {
	      hap.leg.x <- input$hapLenX
	      hap.leg.y <- NULL
	    } else {
	      hap.leg.x <- as.numeric(input$hapLenX)
	      hap.leg.y <- as.numeric(input$hapLenY)
	    }
	    myPos <- anaReg(input$regH)
	    snp.reg <- fetchSnp(chr=myPos$chr, start=myPos$start - input$hapUp * 1000, end=myPos$end + input$hapDown * 1000,
	                        mutType=input$hap_mut_group)[[1]]
	    
	    if (input$uploadHAP == 1) {
	      hap.snp.site <- NULL
	    } else {
	      if (!is.null(input$HAP.snpsite)) {
	        hap.snp.site <- readLines(input$HAP.snpsite$datapath)
	      } else {
	        hap.snp.site <- NULL
	      }
	    }
	    
	    hapNet(data = snp.reg, legend.x=hap.leg.x, legend.y=hap.leg.y,  
	           labels=c(FALSE, TRUE)[as.numeric(input$hapLab)],
	           show.mutation=as.numeric(input$hapMut), scale.ratio=input$hapScale,
	           min.freq=input$hapMin, max.freq=input$hapMax, lwd=input$hapLinkWd,
	           col.link=gsub("0x","#", input$hapLinkCol), pop.list=as.numeric(input$hapPop), snpSites = hap.snp.site)
	    dev.off()
	      
	    })
	    
	  }, contentType = 'image/svg')
	
	# Download haplotypes in NEXUS format
	output$downloadHap.nex <- downloadHandler(
	  filename = function() { "hap.res.nex" },
	  content = function(file) {
	    withProgress(message='Calculation in progress...',value = 0, detail = 'This may take a while...', {
	      myPos <- anaReg(input$regH)
	    snp.reg <- fetchSnp(chr=myPos$chr, start=myPos$start - input$hapUp * 1000, end=myPos$end + input$hapDown * 1000,
	                        mutType=input$hap_mut_group)[[1]]
	    
	    if (input$uploadHAP == 1) {
	      hap.snp.site <- NULL
	    } else {
	      if (!is.null(input$HAP.snpsite)) {
	        hap.snp.site <- readLines(input$HAP.snpsite$datapath)
	      } else {
	        hap.snp.site <- NULL
	      }
	    }
	    
	    writeLines(writeHap(hapConten(data = snp.reg, min.freq=input$hapMin, 
	                                  max.freq=input$hapMax, 
	                                  pop.list=as.numeric(input$hapPop), snpSites = hap.snp.site)), file)
	      
	    })
	    
	  })
	
	# haplotype geo distribution
	observe({
	  if (input$submit3>0) {
	    isolate({
	      myPos <- anaReg(input$regH)
	      hap.pop <- as.numeric(input$hapPop)
	      hap.min <- input$hapMin
	      hap.max <- input$hapMax
	      
	      if (validReg(myPos)) {
	        snp.reg <- fetchSnp(chr=myPos$chr, start=myPos$start - input$hapUp * 1000, end=myPos$end + input$hapDown * 1000,
	                            mutType=input$hap_mut_group)[[1]]
	        
	        if (nrow(snp.reg) < 5) {
	          shinyWidgets::sendSweetAlert(
	            session = session,
	            title = "Error input!", type = "error",
	            text = "Too few SNPs in specified genomic region!"
	          )
	        } else {
	          output$hapGeo <- plotly::renderPlotly({
	            
	            
	            if (input$uploadHAP == 1) {
	              hap.snp.site <- NULL
	            } else {
	              if (!is.null(input$HAP.snpsite)) {
	                hap.snp.site <- readLines(input$HAP.snpsite$datapath)
	              } else {
	                hap.snp.site <- NULL
	              }
	            }
	            
	            hapGeo(haplotype = hapConten(data = snp.reg, min.freq=hap.min, 
	                                         max.freq=hap.max, 
	                                         pop.list=hap.pop, snpSites = hap.snp.site))
	          })
	        }
	      } else {
	        shinyWidgets::sendSweetAlert(
	          session = session,
	          title = "Error input!", type = "error",
	          text = "No SNPs are detected in the specified genomic region or the specified genomic region is too large!"
	        )
	      }
	    })
	  } else {
	      NULL
	  }
	})
	
	## Download PDF file of haplotype geo distribution
	output$downloadHapSta.pdf <- downloadHandler(
	  filename <- function() { paste('hapGeoDis.pdf') },
	  content <- function(file) {
	    withProgress(message='Calculation in progress...',value = 0, detail = 'This may take a while...', {
	      pdf(file, width = 750/72, height = 550/72)
	    
	    myPos <- anaReg(input$regH)
	    snp.reg <- fetchSnp(chr=myPos$chr, start=myPos$start - input$hapUp * 1000, end=myPos$end + input$hapDown * 1000,
	                        mutType=input$hap_mut_group)[[1]]
	    
	    if (input$uploadHAP == 1) {
	      hap.snp.site <- NULL
	    } else {
	      if (!is.null(input$HAP.snpsite)) {
	        hap.snp.site <- readLines(input$HAP.snpsite$datapath)
	      } else {
	        hap.snp.site <- NULL
	      }
	    }
	    
	    grid::grid.draw(hapGeoStatic(haplotype = hapConten(data = snp.reg, min.freq=input$hapMin, 
	                                           max.freq=input$hapMax, 
	                                           pop.list=as.numeric(input$hapPop), snpSites = hap.snp.site)))
	    
	    dev.off()
	      
	    })
	    
	  }, contentType = 'application/pdf')
	
	## Download SVG file of haplotype geo distribution
	output$downloadHapSta.svg <- downloadHandler(
	  filename <- function() { paste('hapGeoDis.svg') },
	  content <- function(file) {
	    withProgress(message='Calculation in progress...',value = 0, detail = 'This may take a while...', {
	      svg(file, width = 750/72, height = 550/72)
	    
	    myPos <- anaReg(input$regH)
	    snp.reg <- fetchSnp(chr=myPos$chr, start=myPos$start - input$hapUp * 1000, end=myPos$end + input$hapDown * 1000,
	                        mutType=input$hap_mut_group)[[1]]
	    
	    if (input$uploadHAP == 1) {
	      hap.snp.site <- NULL
	    } else {
	      if (!is.null(input$HAP.snpsite)) {
	        hap.snp.site <- readLines(input$HAP.snpsite$datapath)
	      } else {
	        hap.snp.site <- NULL
	      }
	    }
	    
	    grid::grid.draw(hapGeoStatic(haplotype = hapConten(data = snp.reg, min.freq=input$hapMin, 
	                                                 max.freq=input$hapMax, 
	                                                 pop.list=as.numeric(input$hapPop), snpSites = hap.snp.site)))
	    
	    dev.off()
	      
	    })
	    
	  }, contentType = 'image/svg')
	
	# Diversity
	observe({
	  if (input$submit4>0) {
	    isolate({
	      div.height <<- input$divHeight
	      div.width <<- input$divWidth
	      
	      myPos <- anaReg(input$regD)
	      
	      if (validReg(myPos)) {
	        div.up <- input$divUp * 1000
	        div.down <- input$divDown * 1000
	        div.group <- input$div_acc_group
	        div.step <- input$snpnumD
	        div.numerator <- input$nuc_numerator
	        div.denominator <- input$nuc_denominator
	        div.mut.group <- input$div_mut_group
	        
	        if (input$uploadDIV == 1) {
	          div.snp.site <- NULL
	        } else {
	          if (!is.null(input$DIV.snpsite)) {
	            div.snp.site <- readLines(input$DIV.snpsite$datapath)
	          } else {
	            div.snp.site <- NULL
	          }
	        }
	        
	        if (!is.null(myPos)) {
	          snp.reg <- fetchSnp(chr=myPos$chr, start=myPos$start - div.up, end=myPos$end + div.down,
	                              mutType=input$div_mut_group)[[1]]
	        } else {
	          snp.reg <- NULL
	        }
	        
	        if (is.null(snp.reg) || nrow(snp.reg) < 10) {
	          js_string <- 'alert("No SNPs are detected in the specified genomic region or the specified genomic region is too large!");'
	          session$sendCustomMessage(type='jsCode', list(value = js_string))
	        } else {
	          nuc.div.plot <<- NULL
	          output$diversity <- renderPlot({
	            nuc.div.plot <<- nucDiv(chr=myPos$chr, nuc.start=myPos$start - div.up, nuc.end=myPos$end + div.down, 
	                                    groups = div.group, step = div.step,
	                                    numerator = div.numerator, denominator = div.denominator, 
	                                    mutType = div.mut.group, snpSites = div.snp.site)
	            grid::grid.draw(gridExtra::grid.arrange(nuc.div.plot[[1]], nuc.div.plot[[2]], ncol=1, heights=c(2.3, 1)))
	          }, height = div.height, width = div.width)
	          
	          ## Download PDF file of Diversity
	          output$downloadDiv01 <- renderUI({
	            req(input$submit4, nuc.div.plot)
	            downloadButton("downloadDiv.pdf", "Download pdf-file")
	          })
	          
	          output$downloadDiv.pdf <- downloadHandler(
	            filename <- function() { paste('diversity.pdf') },
	            content <- function(file) {
	              pdf(file, width = input$divWidth/72, height = input$divHeight/72)
	              grid::grid.draw(gridExtra::grid.arrange(nuc.div.plot[[1]], nuc.div.plot[[2]], ncol=1, heights=c(2.3, 1)))
	              
	              dev.off()
	            }, contentType = 'application/pdf')
	          
	          ## Download SVG file of Diversity
	          output$downloadDiv02 <- renderUI({
	            req(input$submit4, nuc.div.plot)
	            downloadButton("downloadDiv.svg", "Download svg-file")
	          })
	          
	          output$downloadDiv.svg <- downloadHandler(
	            filename <- function() { paste('diversity.svg') },
	            content <- function(file) {
	              svg(file, width = input$divWidth/72, height = input$divHeight/72)
	              grid::grid.draw(gridExtra::grid.arrange(nuc.div.plot[[1]], nuc.div.plot[[2]], ncol=1, heights=c(2.3, 1)))
	              
	              dev.off()
	            }, contentType = 'image/svg')
	          
	        }
	      } else {
	        shinyWidgets::sendSweetAlert(
	          session = session,
	          title = "Error input!", type = "error",
	          text = "Please input genomic region or gene model in appropriate format!"
	        )
	      }
	    })
	  } else {
	    NULL
	  }
	})
	
	observe({
	  if (input$clearDiv>0) {
	    isolate({
	      updateTextInput(session, "regD", value="")
	    })
	  } else {NULL}
	})
	
	observe({
	  if (input$DivExam >0) {
	    isolate({
	      updateTextInput(session, "regD", value="chr07:2839000-2840000")
	    })
	  } else {NULL}
	})
	
	## Download TXT file of diversity
	output$downloadDiv03 <- renderUI({
	  req(input$submit4)
	  downloadButton("downloadDiv.txt", "Download TXT-file")
	})
	
	output$downloadDiv.txt <- downloadHandler(
	  filename <- function() { paste('diversity.txt') },
	  content <- function(file) {
	    write.table(diVTxt, file, sep="\t", quote=F, row.names = F)
	  }, contentType = 'text/plain')
	
	
	# phylogenetics
	observe({
	  if (input$submit5>0) {
	    isolate({
	      phy.height <<- input$phyHeight
	      phy.width <<- input$phyWidth
	      phy.up <- input$phyUp * 1000
	      phy.down <- input$phyDown * 1000
	      
	      myPos <- anaReg(input$regP)
	      
	      if (validReg(myPos)) {
	        phy.acc <- input$mychooserPhy$selected
	        phy.mut.group <- input$phy_mut_group
	        
	        if (input$uploadPHY == 1) {
	          phy.snp.site <- NULL
	        } else {
	          if (!is.null(input$PHY.snpsite)) {
	            phy.snp.site <- readLines(input$PHY.snpsite$datapath)
	          } else {
	            phy.snp.site <- NULL
	          }
	        }
	        
	        if (!is.null(myPos)) {
	          snp.reg <- fetchSnp(chr=myPos$chr, start=myPos$start - phy.up, end=myPos$end + phy.down,
	                              accession=phy.acc, mutType=phy.mut.group)[[1]]
	        } else {
	          snp.reg <- NULL
	        }
	        
	        if (is.null(snp.reg) || nrow(snp.reg) < 10) {
	          js_string <- 'alert("No SNPs are detected in the specified genomic region or the specified genomic region is too large!");'
	          session$sendCustomMessage(type='jsCode', list(value = js_string))
	        } else {
	          output$phylo <- renderPlot({
	            phylo(chr=myPos$chr, start=myPos$start - phy.up, end=myPos$end + phy.down,
	                  accession=phy.acc, mutType=phy.mut.group, snpSites = phy.snp.site)
	          }, height = phy.height, width = phy.width)
	        }
	      } else {
	        shinyWidgets::sendSweetAlert(
	          session = session,
	          title = "Error input!", type = "error",
	          text = "Please input genomic region or gene model in appropriate format!"
	        )
	      }
	    })
	  } else {
	    NULL
	  }
	})
	
	observe({
	  if (input$clearPhy>0) {
	    isolate({
	      updateTextInput(session, "regP", value="")
	    })
	  } else {NULL}
	})
	
	observe({
	  if (input$PhyExam >0) {
	    isolate({
	      updateTextInput(session, "regP", value="chr07:2839000-2840000")
	    })
	  } else {NULL}
	})
	
	## Download PDF file of phylogenetics
	output$downloadPhy01 <- renderUI({
	  req(input$submit5)
	  downloadButton("downloadPhylo.pdf", "Download pdf-file")
	})
	
	output$downloadPhylo.pdf <- downloadHandler(
	  filename <- function() { paste('phylogenetics.pdf') },
	  content <- function(file) {
	    pdf(file, width = input$phyWidth/72, height = input$phyHeight/72)
	    print(figurecp)
	    dev.off()
	  }, contentType = 'application/pdf')
	
	## Download NWK file of phylogenetics
	output$downloadPhy02 <- renderUI({
	  req(input$submit5)
	  downloadButton("downloadPhylo.nwk", "Download Newick-file")
	})
	
	output$downloadPhylo.nwk <- downloadHandler(
	  filename <- function() { paste('phylogenetics.nwk') },
	  content <- function(file) {
	    ape::write.tree(treNwk, file)
	  }, contentType = 'text/plain')
	
	# allele frequency
	observe({
	  if (input$submitaf1>0) {
	    isolate({
	      in.snpid <- unlist(strsplit(input$af_snp_site, split="\\n"))
	      in.snpid <- gsub("^\\s+", "", in.snpid)
	      in.snpid <- gsub("\\s+$", "", in.snpid)
	      in.snpid <- in.snpid[in.snpid!=""]
	      
	      af.group <- input$af_acc_group
	      in.af.col <- unlist(strsplit(input$afCol, split=","))
	      in.af.col <- gsub("^\\s+", "", in.af.col)
	      in.af.col <- gsub("\\s+$", "", in.af.col)
	      
	      af.height <<- input$afHeight
	      af.width <<- input$afWidth
	      
	      output$alleleFreq <- renderPlot({
	        alleleFreq(
	          snpSite = in.snpid,
	          accGroup = af.group,
	          pieCols = in.af.col
	        )
	      }, height = af.height, width = af.width)
	        
	    })
	  } else {
	    NULL
	  }
	})
	
	observe({
	  if (input$clearAf>0) {
	    isolate({
	      updateTextAreaInput(session, "af_snp_site", value="")
	    })
	  } else {NULL}
	})
	
	observe({
	  if (input$AfExam >0) {
	    isolate({
	      updateTextAreaInput(session, "af_snp_site", value="0602942293\n0138383182\n0329584501\n0316733111")
	    })
	  } else {NULL}
	})
	
	## Download PDF file of allele frequency
	output$downloadAfq01 <- renderUI({
	  req(input$submitaf1)
	  downloadButton("downloadAlleleFreq.pdf", "Download pdf-file")
	})
	
	output$downloadAlleleFreq.pdf <- downloadHandler(
	  filename <- function() { paste('alleleFreq.pdf') },
	  content <- function(file) {
	    pdf(file, width = input$afWidth/72, height = input$afHeight/72)
	    
	    in.snpid <- unlist(strsplit(input$af_snp_site, split="\\n"))
	    in.snpid <- gsub("^\\s+", "", in.snpid)
	    in.snpid <- gsub("\\s+$", "", in.snpid)
	    in.snpid <- in.snpid[in.snpid!=""]
	    
	    af.group <- input$af_acc_group
	    in.af.col <- unlist(strsplit(input$afCol, split=","))
	    in.af.col <- gsub("^\\s+", "", in.af.col)
	    in.af.col <- gsub("\\s+$", "", in.af.col)
	    alleleFreq(
	      snpSite = in.snpid,
	      accGroup = af.group,
	      pieCols = in.af.col
	    )
	    
	    dev.off()
	  }, contentType = 'application/pdf')
	
	## Download SVG file of allele frequency
	output$downloadAfq02 <- renderUI({
	  req(input$submitaf1)
	  downloadButton("downloadAlleleFreq.svg", "Download svg-file")
	})
	
	output$downloadAlleleFreq.svg <- downloadHandler(
	  filename <- function() { paste('alleleFreq.svg') },
	  content <- function(file) {
	    svg(file, width = input$afWidth/72, height = input$afHeight/72)
	    
	    in.snpid <- unlist(strsplit(input$af_snp_site, split="\\n"))
	    in.snpid <- gsub("^\\s+", "", in.snpid)
	    in.snpid <- gsub("\\s+$", "", in.snpid)
	    in.snpid <- in.snpid[in.snpid!=""]
	    
	    af.group <- input$af_acc_group
	    in.af.col <- unlist(strsplit(input$afCol, split=","))
	    in.af.col <- gsub("^\\s+", "", in.af.col)
	    in.af.col <- gsub("\\s+$", "", in.af.col)
	    alleleFreq(
	      snpSite = in.snpid,
	      accGroup = af.group,
	      pieCols = in.af.col
	    )
	    
	    dev.off()
	  }, contentType = 'image/svg')
	
	## Download TXT file of allele frequency
	output$downloadAfq03 <- renderUI({
	  req(input$submitaf1)
	  downloadButton("downloadAlleleFreq.txt", "Download TXT-file")
	})
	
	output$downloadAlleleFreq.txt <- downloadHandler(
	  filename <- function() { paste('alleleFreq.txt') },
	  content <- function(file) {
	    in.snpid <- unlist(strsplit(input$af_snp_site, split="\\n"))
	    in.snpid <- gsub("^\\s+", "", in.snpid)
	    in.snpid <- gsub("\\s+$", "", in.snpid)
	    in.snpid <- in.snpid[in.snpid!=""]
	    
	    af.group <- input$af_acc_group
	    in.af.col <- unlist(strsplit(input$afCol, split=","))
	    in.af.col <- gsub("^\\s+", "", in.af.col)
	    in.af.col <- gsub("\\s+$", "", in.af.col)
	    AF.txt <- alleleFreq(
	      snpSite = in.snpid,
	      accGroup = af.group,
	      pieCols = in.af.col
	    )
	    
	    AF.txt.mat <-do.call(cbind, AF.txt)
	    colnames(AF.txt.mat) <- paste0(rep(in.snpid, each=2), ":", colnames(AF.txt.mat))
	    write.table(AF.txt.mat, file, sep = "\t", quote=FALSE, row.names = T, col.names = T)
	  }, contentType = 'text/plain')
  
	# accession information
	output$acc.info.txt <- downloadHandler(
	  filename = function() { "acc.info.txt" },
	  content = function(file) {
	    write.table(acc.info, file, sep = "\t", quote=FALSE, row.names = FALSE)
	}, contentType = 'text/plain')
	
	
	output$sel.acc.info.txt <- downloadHandler(
	  filename = function() { "sel.acc.info.txt" },
	  content = function(file) {
	    accession <- input$mychooserA$selected
	    accession <- gsub(",.+", "", accession)
	    accession <- sapply(accession, function(x){
	      if (x %in% c("Aus", "Indica", "IndicaI", "IndicaII", "Japonica", "TeJ", "TrJ", "Or-I", "Or-II", "Or-III")) {
	        x.dat <- readLines(paste0("./data/", x, ".acc.txt"))
	        return(x.dat)
	      } else {
	        return(x)
	      }
	    })
	    accession <- unique(unlist(accession))
	    
	    write.table(acc.info[acc.info$ID %in% accession, ], 
	                file, sep = "\t", quote=FALSE, row.names = FALSE)
	  }, contentType = 'text/plain')
	
	output$accDis <- plotly::renderPlotly({
	  g <- list(
	    scope = 'world',
	    projection = list(type = 'Equirectangular'),
	    showland = TRUE,
	    showocean=TRUE,
	    showcountries = TRUE,
	    showsubunits = TRUE,
	    landcolor = "white",
	    oceancolor = plotly::toRGB("gray90"),
	    subunitwidth = 1,
	    countrywidth = 0.5,
	    subunitcolor = "blue",
	    countrycolor = "gray85"
	  )
	  
	  acc.info <- acc.info[!is.na(acc.info$Latitude), ]
	  accession <- input$mychooserA$selected
	  accession <- gsub(",.+", "", accession)
	  accession <- sapply(accession, function(x){
	    if (x %in% c("Aus", "Indica", "IndicaI", "IndicaII", "Japonica", "TeJ", "TrJ", "Or-I", "Or-II", "Or-III")) {
	      x.dat <- readLines(paste0("./data/", x, ".acc.txt"))
	      return(x.dat)
	    } else {
	      return(x)
	    }
	  })
	  accession <- unique(unlist(accession))
	  
	  acc.info <- acc.info[acc.info$ID %in% accession, ]
	  acc.info$Ecotype.n <- acc.info$Ecotype
	  acc.info$Ecotype.n[acc.info$Ecotype.n %in% c("Ind_Int", "IndI", "indica", "IndII")] <- "Ind"
	  acc.info$Ecotype.n[acc.info$Ecotype.n %in% c("Jap_Int", "TrJ", "TeJ")] <- "Jap"
	  acc.info$Ecotype.n[acc.info$Ecotype.n %in% c("Int", "VI")] <- "Other"
	  
	  
	  acc.info$Name[is.na(acc.info$Name)] <- ""
	  acc.info <- acc.info[!is.na(acc.info$Ecotype), ]
	  acc.info$text <- paste(acc.info[, 1], acc.info[,2], acc.info[,3], acc.info[,4], acc.info[,7], sep=", ")
	  acc.info.n <- acc.info %>% group_by(Latitude, Longitude) %>% summarise(maptext = paste(text, collapse = "<br>"))
	  acc.info.nn <- merge(acc.info.n, acc.info, by=c("Latitude", "Longitude"))
	  
	  plotly::plot_geo(acc.info.nn, lat = ~Latitude, lon = ~Longitude) %>%
	    plotly::add_markers(
	      marker=list(size=4, color="red"),
	      hovertext = ~maptext
	    ) %>% plotly::layout(title = 'Geographic distribution of selected rice accessions', geo = g)
	 
	})
	
	## Download PDF file of accession distribution
	output$downloadAccDis.pdf <- downloadHandler(
	  filename <- function() { paste('accDis.pdf') },
	  content <- function(file) {
	    withProgress(message='Calculation in progress...',value = 0, detail = 'This may take a while...', {
	      pdf(file, width = 750/72, height = 550/72)
	    
	    acc.info <- acc.info[!is.na(acc.info$Latitude), ]
	    accession <- input$mychooserA$selected
	    accession <- gsub(",.+", "", accession)
	    accession <- sapply(accession, function(x){
	      if (x %in% c("Aus", "Indica", "IndicaI", "IndicaII", "Japonica", "TeJ", "TrJ", "Or-I", "Or-II", "Or-III")) {
	        x.dat <- readLines(paste0("./data/", x, ".acc.txt"))
	        return(x.dat)
	      } else {
	        return(x)
	      }
	    })
	    accession <- unique(unlist(accession))
	    
	    acc.info <- acc.info[acc.info$ID %in% accession, ]
	    acc.info$Ecotype.n <- acc.info$Ecotype
	    acc.info$Ecotype.n[acc.info$Ecotype.n %in% c("Ind_Int", "IndI", "indica", "IndII")] <- "Ind"
	    acc.info$Ecotype.n[acc.info$Ecotype.n %in% c("Jap_Int", "TrJ", "TeJ")] <- "Jap"
	    acc.info$Ecotype.n[acc.info$Ecotype.n %in% c("Int", "VI")] <- "Other"
	    
	    load("./data/worldmap.RData")
	    
	    acc.info$Name[is.na(acc.info$Name)] <- ""
	    acc.info <- acc.info[!is.na(acc.info$Ecotype), ]
	    
	    mp <- mp + ggplot2::geom_point(ggplot2::aes(x=Longitude, y=Latitude, color=Ecotype.n), size=0.5, data=acc.info) + 
	      ggplot2::scale_x_continuous("", breaks=NULL) + ggplot2::scale_y_continuous("", breaks=NULL) 
	    mp <- mp + ggplot2::guides(color=ggplot2::guide_legend(title=NULL))
	    grid::grid.draw(mp)
	    
	    dev.off()
	      
	    })
	    
	  }, contentType = 'application/pdf')
	
	## Download SVG file of accession distribution
	output$downloadAccDis.svg <- downloadHandler(
	  filename <- function() { paste('accDis.svg') },
	  content <- function(file) {
	    withProgress(message='Calculation in progress...',value = 0, detail = 'This may take a while...', {
	      svg(file, width = 750/72, height = 550/72)
	    
	    acc.info <- acc.info[!is.na(acc.info$Latitude), ]
	    accession <- input$mychooserA$selected
	    accession <- gsub(",.+", "", accession)
	    accession <- sapply(accession, function(x){
	      if (x %in% c("Aus", "Indica", "IndicaI", "IndicaII", "Japonica", "TeJ", "TrJ", "Or-I", "Or-II", "Or-III")) {
	        x.dat <- readLines(paste0("./data/", x, ".acc.txt"))
	        return(x.dat)
	      } else {
	        return(x)
	      }
	    })
	    accession <- unique(unlist(accession))
	    
	    acc.info <- acc.info[acc.info$ID %in% accession, ]
	    acc.info$Ecotype.n <- acc.info$Ecotype
	    acc.info$Ecotype.n[acc.info$Ecotype.n %in% c("Ind_Int", "IndI", "indica", "IndII")] <- "Ind"
	    acc.info$Ecotype.n[acc.info$Ecotype.n %in% c("Jap_Int", "TrJ", "TeJ")] <- "Jap"
	    acc.info$Ecotype.n[acc.info$Ecotype.n %in% c("Int", "VI")] <- "Other"
	    
	    load("./data/worldmap.RData")
	    
	    acc.info$Name[is.na(acc.info$Name)] <- ""
	    acc.info <- acc.info[!is.na(acc.info$Ecotype), ]
	    
	    mp <- mp + ggplot2::geom_point(ggplot2::aes(x=Longitude, y=Latitude, color=Ecotype.n), size=0.1, data=acc.info) + 
	      ggplot2::scale_x_continuous("", breaks=NULL) + ggplot2::scale_y_continuous("", breaks=NULL) 
	    mp <- mp + ggplot2::guides(color=ggplot2::guide_legend(title=NULL))
	    grid::grid.draw(mp)
	    
	    dev.off()
	      
	    })
	    
	  }, contentType = 'application/svg')
	
	output$mytable1 = renderDataTable({
	  accession <- input$mychooserA$selected
	  accession <- gsub(",.+", "", accession)
	  accession <- sapply(accession, function(x){
	    if (x %in% c("Aus", "Indica", "IndicaI", "IndicaII", "Japonica", "TeJ", "TrJ", "Or-I", "Or-II", "Or-III")) {
	      x.dat <- readLines(paste0("./data/", x, ".acc.txt"))
	      return(x.dat)
	    } else {
	      return(x)
	    }
	  })
	  accession <- unique(unlist(accession))
	  
	  acc.info[acc.info$ID %in% accession, ]
	}, options = list(lengthMenu = c(5, 8, 10), pageLength = 5, searching = TRUE, autoWidth = FALSE), escape = FALSE
	)
	
	# Bulk download genotypes of seleceted SNPs
	observe({
	  if (input$submit6>0) {
	    isolate({
	      myPos <- anaReg(input$regBB)
	      
	      if (validReg(myPos)) {
	        snp.info.down <<- NULL
	        
	        output$mytable2 = renderDataTable({
	          snp.info.down <<- snpInfo(chr=myPos$chr, start=myPos$start, end=myPos$end, 
	                                    accession = input$mychooserD$selected, mutType = input$down_mut_group)
	          snp.info.down[[2]]
	        }, options = list(lengthMenu = c(5, 8, 10), pageLength = 5, searching = TRUE, autoWidth = TRUE), escape = FALSE
	        )
	        
	        output$bulkdownloadsnp.txt <- downloadHandler(
	          filename = function() { "down.snp.geno.txt" },
	          content = function(file) {
	            write.table(snp.info.down[[1]][[1]], file, sep="\t", quote=F)
	          })
	        
	        # Bulk download information of SNPs
	        output$bulkdownloadsnpInfo.txt <- downloadHandler(
	          filename = function() { "down.snp.info.txt" },
	          content = function(file) {
	            write.table(snp.info.down[[2]], file, sep="\t", quote=F, row.names=F)
	          })
	        
	        # Bulk download gene annotation
	        output$bulkdownloadgene.txt <- downloadHandler(
	          filename = function() { "down.gene.info.txt" },
	          content = function(file) {
	            gene.info <- gff[gff$chr==myPos$chr & gff$start>=myPos$start & gff$end<=myPos$end, ]
	            write.table(gene.info, file, sep="\t", quote=F, row.names=F)
	          })
	      } else {
	        shinyWidgets::sendSweetAlert(
	          session = session,
	          title = "Error input!", type = "error",
	          text = "Please input genomic region or gene model in appropriate format!"
	        )
	      }
	    })
	  } else {
	      NULL
	  }
	})
	
	observe({
	  if (input$clearBB>0) {
	    isolate({
	      updateTextInput(session, "regBB", value="")
	    })
	  } else {NULL}
	})
	
	observe({
	  if (input$BBExam >0) {
	    isolate({
	      updateTextInput(session, "regBB", value="chr07:29611303-29669223")
	    })
	  } else {NULL}
	})

})

