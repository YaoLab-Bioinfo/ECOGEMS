
# options(shiny.maxRequestSize = 200*1024^2)

shinyServer(function(input, output, session) {
  
  # GBrowser
  observe({
    if (input$submit1>0) {
      isolate({
        myPos <- anaReg(input$regB)
        
        if (validReg(myPos)) {
        } else {
          js_string <- 'alert("Please input genomic region or gene model in appropriate format!");'
          session$sendCustomMessage(type='jsCode', list(value = js_string))
          #myPos <- anaReg("chr07:29611303-29669223")
          myPos <- NULL
        }
        
        snp.info <- snpInfo(chr=myPos$chr, start=myPos$start - input$GBUP, end=myPos$end + input$GBDOWN, 
                            accession = input$mychooserB$selected, mutType = input$GB_mut_group)
        if (nrow(snp.info[[1]][[1]]) < 1) {
          js_string <- 'alert("No SNPs in specified genomic region!");'
          session$sendCustomMessage(type='jsCode', list(value = js_string))
        } else {
          GBplot <<- NULL
          output$gbrowser <- renderPlotly({
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
              grid.draw(GBplot[[1]])
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
        
      })
    } else {
      if (input$regB == "chr07:29611303-29669223") {
        myPos <- anaReg(input$regB)
        
        GBplot <<- NULL
        output$gbrowser <- renderPlotly({
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
            grid.draw(GBplot[[1]])
            dev.off()
          }, contentType = 'application/pdf')
        
        snp.info <- snpInfo(chr=myPos$chr, start=myPos$start - input$GBUP, end=myPos$end + input$GBDOWN, 
                            accession = input$mychooserB$selected, mutType = input$GB_mut_group)
        
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
        
      } else {
        NULL
      }
    }
  })
  

  
  
  
  # LDheatmap
  observe({
    if (input$submit2>0) {
      isolate({
        ld.height <<- input$ldHeight
        ld.width <<- input$ldWidth
        myPos <- anaReg(input$regL)
        
        if (validReg(myPos)) {
        } else {
          js_string <- 'alert("Please input genomic region or gene model in appropriate format!");'
          session$sendCustomMessage(type='jsCode', list(value = js_string))
          myPos <- NULL
        }
        
        snp.reg <- fetchSnp(chr=myPos$chr, start=myPos$start - input$ldUp * 1000, 
                            end=myPos$end + input$ldDown * 1000, accession = input$mychooserLD$selected,
                            mutType = input$ld_mut_group)[[1]]
        if (nrow(snp.reg) < 5) {
          js_string <- 'alert("Too few SNPs in specified genomic region!");'
          session$sendCustomMessage(type='jsCode', list(value = js_string))
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
        
      })
    } else {
      ld.height <<- input$ldHeight
      ld.width <<- input$ldWidth
      
      if (input$uploadLD == 1) {
        ld.snp.site <- NULL
      } else {
        if (!is.null(input$LD.snpsite)) {
          ld.snp.site <- readLines(input$LD.snpsite$datapath)
        } else {
          ld.snp.site <- NULL
        }
      }
      
      if (input$regL == "LOC_Os11g35500") {
        myPos <- anaReg(input$regL)
        snp.pos <- as.numeric(unlist(strsplit(input$ldpos, split=",")))
        
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
      } else {
        NULL
      }
    }
  })

	## Download PDF file of LDheatmap
	output$downloadLD.pdf <- downloadHandler(
	  filename <- function() { paste('LDheatmap.pdf') },
	  content <- function(file) {
	    withProgress(message='Calculation in progress...',value = 0, detail = 'This may take a while...', {
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
	      
	    })
	    
	  }, contentType = 'application/pdf')
	
	## Download SVG file of LDheatmap
	output$downloadLD.svg <- downloadHandler(
	  filename <- function() { paste('LDheatmap.svg') },
	  content <- function(file) {
	    withProgress(message='Calculation in progress...',value = 0, detail = 'This may take a while...', {
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
	      
	    })
	    
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
	      } else {
	        js_string <- 'alert("Please input genomic region or gene model in appropriate format!");'
	        session$sendCustomMessage(type='jsCode', list(value = js_string))
	        myPos <- NULL
	      }
	      
	      snp.reg <- fetchSnp(chr=myPos$chr, start=myPos$start - hap.up, end=myPos$end + hap.down, mutType=hap.mut.grp)[[1]]
	      
	      if (nrow(snp.reg) < 5) {
	        js_string <- 'alert("Too few SNPs in specified genomic region!");'
	        session$sendCustomMessage(type='jsCode', list(value = js_string))
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
	     
	    })
	  } else {
	    isolate({
	      if (input$regH == "LOC_Os10g32600") {
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
	        
	        output$haplotype <- renderPlot({
	          hapNet(data = snp.reg, labels=FALSE, show.mutation=0, scale.ratio=100,
	                 min.freq=50, max.freq=2508, legend.x="bottomleft", snpSites = hap.snp.site)
	        }, height = 550, width = 750)
	      } else {
	        NULL
	      }
	    })
	  }
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
	      
	      snp.reg <- fetchSnp(chr=myPos$chr, start=myPos$start - input$hapUp * 1000, end=myPos$end + input$hapDown * 1000,
	                          mutType=input$hap_mut_group)[[1]]
	      
	      if (nrow(snp.reg) < 5) {
	        js_string <- 'alert("Too few SNPs in specified genomic region!");'
	        session$sendCustomMessage(type='jsCode', list(value = js_string))
	      } else {
	        output$hapGeo <- renderPlotly({
	          
	          
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
	      
	    })
	  } else {
	    isolate({
	      if (input$regH == "LOC_Os10g32600") { 
	        myPos <- anaReg(input$regH)
	        if (input$uploadHAP == 1) {
	          hap.snp.site <- NULL
	        } else {
	          if (!is.null(input$HAP.snpsite)) {
	            hap.snp.site <- readLines(input$HAP.snpsite$datapath)
	          } else {
	            hap.snp.site <- NULL
	          }
	        }
	        
	        output$hapGeo <- renderPlotly({
	          snp.reg <- fetchSnp(chr=myPos$chr, start=myPos$start - input$hapUp * 1000, end=myPos$end + input$hapDown * 1000,
	                              mutType=input$hap_mut_group)[[1]]
	          
	          hapGeo(haplotype = hapConten(data = snp.reg, min.freq=50, max.freq=2508, snpSites = hap.snp.site))
	        })
	      } else {
	        NULL
	      }
	    })
	  }
	})
	
	## Download PDF file of haplotype geo distribution
	output$downloadHapSta.pdf <- downloadHandler(
	  filename <- function() { paste('hapGeoDis.pdf') },
	  content <- function(file) {
	    withProgress(message='Calculation in progress...',value = 0, detail = 'This may take a while...', {
	      pdf(file, width = input$divWidth/72, height = input$divHeight/72)
	    
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
	    
	    grid.draw(hapGeoStatic(haplotype = hapConten(data = snp.reg, min.freq=input$hapMin, 
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
	      svg(file, width = input$divWidth/72, height = input$divHeight/72)
	    
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
	    
	    grid.draw(hapGeoStatic(haplotype = hapConten(data = snp.reg, min.freq=input$hapMin, 
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
	      div.height2 <<- input$divHeight2
	      div.width2 <<- input$divWidth2
	      
	      myPos <- anaReg(input$regD)
	      
	      if (validReg(myPos)) {
	      } else {
	        js_string <- 'alert("Please input genomic region or gene model in appropriate format!");'
	        session$sendCustomMessage(type='jsCode', list(value = js_string))
	        myPos <- NULL
	      }
	      
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
	      
	      snp.reg <- fetchSnp(chr=myPos$chr, start=myPos$start - div.up, end=myPos$end + div.down,
	                          mutType=input$hap_mut_group)[[1]]
	      
	      if (nrow(snp.reg) < 10) {
	        js_string <- 'alert("Too few SNPs in specified genomic region!");'
	        session$sendCustomMessage(type='jsCode', list(value = js_string))
	      } else {
	        nuc.div.plot <<- NULL
	        output$diversity <- renderPlot({
	          nuc.div.plot <<- nucDiv(chr=myPos$chr, nuc.start=myPos$start - div.up, nuc.end=myPos$end + div.down, 
	                   groups = div.group, step = div.step,
	                   numerator = div.numerator, denominator = div.denominator, 
	                   mutType = div.mut.group, snpSites = div.snp.site)
	          grid.draw(grid.arrange(nuc.div.plot[[1]], nuc.div.plot[[2]], ncol=1, heights=c(2.3, 1)))
	        }, height = div.height, width = div.width)
	        
	        ## Download PDF file of Diversity
	        output$downloadDiv.pdf <- downloadHandler(
	          filename <- function() { paste('diversity.pdf') },
	          content <- function(file) {
	            pdf(file, width = input$divWidth/72, height = input$divHeight/72)
	            grid.draw(grid.arrange(nuc.div.plot[[1]], nuc.div.plot[[2]], ncol=1, heights=c(2.3, 1)))
	            
	            dev.off()
	          }, contentType = 'application/pdf')
	        
	        ## Download SVG file of Diversity
	        output$downloadDiv.svg <- downloadHandler(
	          filename <- function() { paste('diversity.svg') },
	          content <- function(file) {
	            svg(file, width = input$divWidth/72, height = input$divHeight/72)
	            grid.draw(grid.arrange(nuc.div.plot[[1]], nuc.div.plot[[2]], ncol=1, heights=c(2.3, 1)))
	            
	            dev.off()
	          }, contentType = 'image/svg')
	        
	      }

	    })
	  } else {
	    NULL
	  }
	})
	
	## Download TXT file of diversity
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
	      
#	      withProgress(message='Calculation in progress...',value = 0, detail = 'This may take a while...', {
	        myPos <- anaReg(input$regP)
	      
	      if (validReg(myPos)) {
	      } else {
	        js_string <- 'alert("Please input genomic region or gene model in appropriate format!");'
	        session$sendCustomMessage(type='jsCode', list(value = js_string))
	        myPos <- NULL
	      }
	      
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
	      
	      snp.reg <- fetchSnp(chr=myPos$chr, start=myPos$start - phy.up, end=myPos$end + phy.down,
	                          accession=phy.acc, mutType=phy.mut.group)[[1]]
	      
	      if (nrow(snp.reg) < 10) {
	        js_string <- 'alert("Too few SNPs in specified genomic region!");'
	        session$sendCustomMessage(type='jsCode', list(value = js_string))
	      } else {
	        output$phylo <- renderPlot({
	            phylo(chr=myPos$chr, start=myPos$start - phy.up, end=myPos$end + phy.down,
	                  accession=phy.acc, mutType=phy.mut.group, snpSites = phy.snp.site)
	        }, height = phy.height, width = phy.width)
	      }
	        
#	      })
	      
	    })
	  } else {
	    NULL
	  }
	})
	
	## Download PDF file of phylogenetics
	output$downloadPhylo.pdf <- downloadHandler(
	  filename <- function() { paste('phylogenetics.pdf') },
	  content <- function(file) {
	    pdf(file, width = input$phyWidth/72, height = input$phyHeight/72)
	    print(figurecp)
	    dev.off()
	  }, contentType = 'application/pdf')
	
	## Download NWK file of phylogenetics
	output$downloadPhylo.nwk <- downloadHandler(
	  filename <- function() { paste('phylogenetics.nwk') },
	  content <- function(file) {
	    write.tree(treNwk, file)
	  }, contentType = 'text/plain')
	
	# allele frequency
	observe({
	  if (input$submitaf1>0) {
	    isolate({
#	      withProgress(message='Calculation in progress...',value = 0, detail = 'This may take a while...', {
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
	        
#	      })
	      
	    })
	  } else {
	    NULL
	  }
	})
	
	## Download PDF file of allele frequency
	output$downloadAlleleFreq.pdf <- downloadHandler(
	  filename <- function() { paste('alleleFreq.pdf') },
	  content <- function(file) {
	    pdf(file, width = input$afWidth/72, height = input$afHeight/72)
	    
#	    withProgress(message='Calculation in progress...',value = 0, detail = 'This may take a while...', {
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
	      
#	    })
	    
	    dev.off()
	  }, contentType = 'application/pdf')
	
	## Download SVG file of allele frequency
	output$downloadAlleleFreq.svg <- downloadHandler(
	  filename <- function() { paste('alleleFreq.svg') },
	  content <- function(file) {
	    svg(file, width = input$afWidth/72, height = input$afHeight/72)
	    
#	    withProgress(message='Calculation in progress...',value = 0, detail = 'This may take a while...', {
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
	      
#	    })
	    
	    dev.off()
	  }, contentType = 'image/svg')
  
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
	
	output$accDis <- renderPlotly({
	  g <- list(
	    scope = 'world',
	    projection = list(type = 'Equirectangular'),
	    showland = TRUE,
	    showocean=TRUE,
	    showcountries = TRUE,
	    showsubunits = TRUE,
	    landcolor = "white",
	    oceancolor = toRGB("gray90"),
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
	  
	  plot_geo(acc.info.nn, lat = ~Latitude, lon = ~Longitude) %>%
	    add_markers(
	      marker=list(size=4, color="red"),
	      hovertext = ~maptext
	    ) %>% layout(title = 'Geographic distribution of selected rice accessions', geo = g)
	 
	})
	
	## Download PDF file of accession distribution
	output$downloadAccDis.pdf <- downloadHandler(
	  filename <- function() { paste('accDis.pdf') },
	  content <- function(file) {
	    withProgress(message='Calculation in progress...',value = 0, detail = 'This may take a while...', {
	      pdf(file, width = input$divWidth/72, height = input$divHeight/72)
	    
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
	    
	    mp <- mp + geom_point(aes(x=Longitude, y=Latitude, color=Ecotype.n), size=0.5, data=acc.info) + 
	      scale_x_continuous("", breaks=NULL) + scale_y_continuous("", breaks=NULL) 
	    mp <- mp + guides(color=guide_legend(title=NULL))
	    grid.draw(mp)
	    
	    dev.off()
	      
	    })
	    
	  }, contentType = 'application/pdf')
	
	## Download SVG file of accession distribution
	output$downloadAccDis.svg <- downloadHandler(
	  filename <- function() { paste('accDis.svg') },
	  content <- function(file) {
	    withProgress(message='Calculation in progress...',value = 0, detail = 'This may take a while...', {
	      svg(file, width = input$divWidth/72, height = input$divHeight/72)
	    
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
	    
	    mp <- mp + geom_point(aes(x=Longitude, y=Latitude, color=Ecotype.n), size=0.1, data=acc.info) + 
	      scale_x_continuous("", breaks=NULL) + scale_y_continuous("", breaks=NULL) 
	    mp <- mp + guides(color=guide_legend(title=NULL))
	    grid.draw(mp)
	    
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
	      } else {
	        js_string <- 'alert("Please input genomic region or gene model in appropriate format!");'
	        session$sendCustomMessage(type='jsCode', list(value = js_string))
	        myPos <- NULL
	      }
	      
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
	      
	    })
	  } else {
	    if (input$regBB == "chr07:29611303-29669223") {
	      isolate({
	        myPos <- anaReg(input$regBB)
	        
	        snp.info.down <<- NULL
	        
	        output$mytable2 = renderDataTable({
	          snp.info.down <<- snpInfo(chr=myPos$chr, start=myPos$start, end=myPos$end, 
	                              accession = input$mychooserD$selected, mutType = input$down_mut_group)
	          snp.info.down[[2]]
	        }, options = list(lengthMenu = c(5, 8, 10), pageLength = 5, searching = TRUE, autoWidth = TRUE), escape = FALSE)
	        
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
	        
	      })
	    } else {
	      NULL
	    }
	  }
	})

})

