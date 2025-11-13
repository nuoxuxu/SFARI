######EDITING code from IsoformSwitchAnalyzeR TO COLOUR FOR AGO UTRs
evalSig <- function (pValue, alphas) {
  sapply(pValue, function(x) {
    if (is.na(x)) {
      return("NA")
    }
    else if (x == -1) {
      sigLevel <- "*"
    }
    else if (x < min(alphas)) {
      sigLevel <- "***"
    }
    else if (x < max(alphas)) {
      sigLevel <- "*"
    }
    else {
      sigLevel <- "ns"
    }
    return(sigLevel)
  })
}

startCapitalLetter <- function (aVec)  {
  paste(toupper(substr(aVec, 1, 1)), substr(aVec, 2, nchar(aVec)), 
        sep = "")
}





AGO1_plotting <- function() {

switchAnalyzeRlist = aSwitchList_part2
  gene = 'AGO1'
  isoform_id = NULL
  
  ### Advanced arguments
  rescaleTranscripts = TRUE
  rescaleRoot = 3
  plotXaxis = !rescaleTranscripts
  reverseMinus = TRUE
  ifMultipleIdenticalAnnotation = 'summarize'
  annotationImportance = c('signal_peptide','protein_domain','idr')
  plotTopology = F
  IFcutoff = 0.05
  abbreviateLocations = TRUE
  rectHegith = 0.2
  codingWidthFactor = 2
  nrArrows = 20
  arrowSize = 0.2
  optimizeForCombinedPlot = FALSE
  condition1 = "t00"
  condition2 = "t30"
  dIFcutoff = 0.1
  alphas = c(0.05, 0.001)
  localTheme = theme_bw()

  ### Check input
  if (TRUE) {
    if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
      stop(
        'The object supplied to \'switchAnalyzeRlist\' must be a \'switchAnalyzeRlist\''
      )
    }
    
    # check isoform and gene name input
    idInfoCheck <- sum(c(is.null(gene), is.null(isoform_id)))
    if (idInfoCheck != 1) {
      if (idInfoCheck == 2) {
        stop('One of \'gene\' or \'isoform_id\' must be given as input')
      }
      if (idInfoCheck == 0) {
        stop('Only one of \'gene\' or \'isoform_id\' can be supplied')
      }
    }
    
    if (!is.logical(rescaleTranscripts)) {
      stop('The \'transformCoordinats\' argument must be either be TRUE or FALSE')
    }
    if (!ifMultipleIdenticalAnnotation %in% c('summarize', 'number','ignore')) {
      stop(
        'the argument \'ifMultipleIdenticalAnnotation\' must be either \'summarize\', \'number\' or \'ignore\' - please see documentation for details'
      )
    }
    okImportance <- c('signal_peptide','protein_domain','idr')
    if( ! all( okImportance %in% intersect(okImportance, annotationImportance)) ) {
      stop('The \'annotationImportance\' must specify the order of all features annotated and only that')
    }
    
    if (optimizeForCombinedPlot) {
      if (is.null(condition1) | is.null(condition2)) {
        stop(
          'When optimizeForCombinedPlot is TRUE, both condition1 and condition2 must also be suppplied'
        )
      }
    }
    
    if( switchAnalyzeRlist$sourceId == 'preDefinedSwitches') {
      if(is.null(condition1)) {
        condition1 <- switchAnalyzeRlist$isoformFeatures$condition_1[1]
      }
      if(is.null(condition2)) {
        condition2 <- switchAnalyzeRlist$isoformFeatures$condition_2[1]
      }
    }
    
    if(plotTopology) {
      if( ! 'topologyAnalysis' %in% names(switchAnalyzeRlist) ) {
        message('Omitting toplogy visualization as it has not been added. You can add this analysis through analyzeDeepTMHMM(). To avoid this message set \"plotTopology=FALSE\"')
        plotTopology <- FALSE
      }
      
    }
    
    isConditional <- ! is.null(condition1)
    hasQuant <- ! all(is.na(switchAnalyzeRlist$isoformFeatures$IF1))
    
    ### Don't plot topology if it is not annotated
    if(plotTopology) {
      plotTopology <- 'topologyAnalysis' %in% names(switchAnalyzeRlist)
    }
    
  }
  
  ### Check for what annotation are stored in the switchAnalyzeRlist
  if (TRUE) {
    if (!is.null(switchAnalyzeRlist$orfAnalysis)) {
      inclORF <- TRUE
      
      if( ! is.null(switchAnalyzeRlist$orfAnalysis$orf_origin) ) {
        if ( any( switchAnalyzeRlist$orfAnalysis$orf_origin == 'not_annotated_yet' )) {
          stop('Some ORFs have not been annotated yet. Please return to the analyzeNovelIsoformORF() step and start again.')
        }
      }
    } else {
      inclORF <- FALSE
      warning(
        'ORFs have not to have been annoated. If ORF should be visualized it can be annoated with the \'annotatePTC()\' function'
      )
    }
    
    if (!is.null(switchAnalyzeRlist$domainAnalysis)) {
      inclDomainAnalysis <- TRUE
      if (!inclORF) {
        warning('Cannot plot annoated protein domians when no ORF are annoated')
        inclDomainAnalysis <- FALSE
      }
    } else {
      inclDomainAnalysis <- FALSE
    }
    
    if (!is.null(switchAnalyzeRlist$idrAnalysis)) {
      inclIdrAnalysis <- TRUE
      if (!inclORF) {
        warning('Cannot plot annoated IDR when no ORF are annoated')
        inclIdrAnalysis <- FALSE
      }
    } else {
      inclIdrAnalysis <- FALSE
    }
    
    if (!is.null(switchAnalyzeRlist$signalPeptideAnalysis)) {
      inclSignalPAnalysis <- TRUE
      if (!inclORF) {
        warning('Cannot plot annoated signal peptides when no ORF are annoated')
        inclSignalPAnalysis <- FALSE
      }
    } else {
      inclSignalPAnalysis <- FALSE
    }
    
    if ('codingPotential' %in%           
        colnames(switchAnalyzeRlist$isoformFeatures)
    ) {
      inclCodingPotential <- TRUE
    } else {
      inclCodingPotential <- FALSE
    }
    
    if ('PTC' %in% colnames(switchAnalyzeRlist$isoformFeatures)) {
      inclPTC <- TRUE
    } else {
      inclPTC <- FALSE
    }
  }
  
  
  
  ### interpret gene and isoform_id input
  if (TRUE) {
    ### Handle gene and isoform_id input
    if (!is.null(gene)) {
      # Decode gene supplied
      if (
        tolower(gene) %in%
        tolower(switchAnalyzeRlist$isoformFeatures$gene_id)
      ) {
        gene_id <- gene
        isoform_id <- unique(
          switchAnalyzeRlist$isoformFeatures$isoform_id[which(
            switchAnalyzeRlist$isoformFeatures$gene_id %in% gene_id
          )]
        )
      } else if (
        tolower(gene) %in%
        tolower(switchAnalyzeRlist$isoformFeatures$gene_name)
      ) {
        gene_id <- unique(
          switchAnalyzeRlist$isoformFeatures$gene_id[which(
            tolower(
              switchAnalyzeRlist$isoformFeatures$gene_name
            ) %in% tolower(gene)
          )])
        
        if (length(gene_id) > 1) {
          stop(
            paste(
              'The gene supplied covers multiple gene_ids (usually due to gene duplications). Currently multigene plotting is not supported. Please use either of the following gene_ids: \'',
              paste(gene_id, collapse = '\', \''),
              '\', and try again',
              sep = ''
            )
          )
        }
        isoform_id <-
          unique(
            switchAnalyzeRlist$isoformFeatures$isoform_id[which(
              switchAnalyzeRlist$isoformFeatures$gene_id %in%
                gene_id
            )]
          )
      } else {
        similarGenes <- c(
          unique(
            switchAnalyzeRlist$isoformFeatures$gene_id[which(
              agrepl(
                gene,
                switchAnalyzeRlist$isoformFeatures$gene_id
              )
            )]
          ),
          unique(
            switchAnalyzeRlist$isoformFeatures$gene_name[which(
              agrepl(
                gene,
                switchAnalyzeRlist$isoformFeatures$gene_name
              )
            )]
          )
        )
        if (length(similarGenes)) {
          stop(
            paste(
              'The gene supplied is not pressent in the switchAnalyzeRlist. did you mean any of: \'',
              paste(similarGenes, collapse = '\', \''),
              '\'',
              sep = ''
            )
          )
        } else {
          stop(
            'The gene supplied is not pressent in the switchAnalyzeRlist, please re-check the name and try again.'
          )
        }
      }
      if (!length(isoform_id)) {
        stop(
          'No isoforms annotated to the supplied gene was found. re-check the name and try again.'
        )
      }
    } else {
      if (any(
        isoform_id %in% switchAnalyzeRlist$isoformFeatures$isoform_id
      )) {
        if (!all(
          isoform_id %in%
          switchAnalyzeRlist$isoformFeatures$isoform_id
        )) {
          notFound <-
            setdiff(isoform_id,
                    switchAnalyzeRlist$isoformFeatures$isoform_id)
          warning(
            paste(
              '\nThe following isoform was not found: \'',
              paste(notFound, collapse = '\', \''),
              '\'. Only the other isoforms will be used\n',
              sep = ''
            )
          )
        }
        
        gene_id <- unique(
          switchAnalyzeRlist$isoformFeatures$gene_id[which(
            switchAnalyzeRlist$isoformFeatures$isoform_id %in%
              isoform_id
          )
          ])
        if (length(gene_id) > 1) {
          stop(
            paste(
              'The isoforms supplied covers multiple gene_ids. Currently multigene plotting is not supported. Please use either of the following gene_ids: \'',
              paste(gene_id, collapse = '\', \''),
              '\', and try again',
              sep = ''
            )
          )
        }
        
      } else {
        stop(
          'Non of the supplied isoforms were found in the switchAnalyzeRlist, please re-check the name and try again'
        )
      }
    }
    
  }
  
  ### Extract isoform and annotation data
  if (TRUE) {
    ### Extract iso annoation
    if(TRUE) {
      columnsToExtract <-
        c(
          'isoform_id',
          'gene_id',
          'gene_name',
          'codingPotential',
          'PTC',
          'class_code',
          'sub_cell_location',
          'solubility_status',
          'isoform_switch_q_value',
          'IF_overall',
          'IF1',
          'IF2',
          'dIF'
        )
      columnsToExtract <-
        na.omit(match(
          columnsToExtract,
          colnames(switchAnalyzeRlist$isoformFeatures)
        ))
      
      ### Extract the rows corresponding to features of interes
      if( isConditional ) {
        rowsToExtract <-
          which(
            switchAnalyzeRlist$isoformFeatures$isoform_id %in% isoform_id &
              switchAnalyzeRlist$isoformFeatures$condition_1 == condition1 &
              switchAnalyzeRlist$isoformFeatures$condition_2 == condition2
          )
      } else {
        rowsToExtract <-
          which(switchAnalyzeRlist$isoformFeatures$isoform_id %in% isoform_id)
      }
      
      isoInfo <-
        unique(switchAnalyzeRlist$isoformFeatures[
          rowsToExtract,
          columnsToExtract
        ])
      
      ### Subset to used isoforms
      if(hasQuant) {
        isoInfo$minIF <- apply(
          isoInfo[,na.omit(match(c('IF1','IF2'), colnames(isoInfo)) ),drop=FALSE],
          1,
          function(x) {
            max(x, na.rm = TRUE)
          }
        )
        isoInfo <- isoInfo[which(
          isoInfo$minIF >= IFcutoff
        ),]
        if(nrow(isoInfo) == 0) {
          stop('No isoforms left after filtering via the "IFcutoff" argument.')
        }
        
        isoform_id <- isoInfo$isoform_id
      }
      if( switchAnalyzeRlist$sourceId == 'preDefinedSwitches' ) {
        isoInfo <- isoInfo[which(
          !is.na(isoInfo$dIF)
        ),]
        
        isoform_id <- intersect(
          isoform_id, isoInfo$isoform_id
        )
      }
      
      ### Remove if all is annotated as NA
      if(!is.null(isoInfo$sub_cell_location)) {
        if(all(is.na(isoInfo$sub_cell_location))) {
          isoInfo$sub_cell_location <- NULL
        }
      }
      if(!is.null(isoInfo$solubility_status)) {
        if(all(is.na(isoInfo$solubility_status))) {
          isoInfo$solubility_status <- NULL
        }
      }
      
      ### Extract ORF info
      columnsToExtract <-
        c('isoform_id', 'orfStartGenomic', 'orfEndGenomic','wasTrimmed', 'trimmedStartGenomic')
      columnsToExtract <-
        na.omit(match(
          columnsToExtract,
          colnames(switchAnalyzeRlist$orfAnalysis)
        ))
      rowsToExtract <-
        which(switchAnalyzeRlist$orfAnalysis$isoform_id %in% isoform_id)
      
      orfInfo <-
        switchAnalyzeRlist$orfAnalysis[rowsToExtract, columnsToExtract]
      
      if(nrow(orfInfo) == 0) {
        warning(
          paste(
            'There might somthing wrong with the switchAnalyzeRlist',
            '- there are no ORF annoation matching the isoforms of interest:',
            paste(isoform_id, collapse = ', '),
            '. These isoforoms will be plotted as non-coding.',
            sep=' '
          )
        )
        
        ### Add NAs manually
        isoInfo$orfStartGenomic <- NA
        isoInfo$orfEndGenomic <- NA
        isoInfo$wasTrimmed <- FALSE
        isoInfo$trimmedStartGenomic <- NA
        
      } else {
        inclTrimmedIsoforms <- FALSE
        if( 'wasTrimmed' %in% colnames(orfInfo) ) {
          if(any( orfInfo$wasTrimmed, na.rm = TRUE) ) {
            inclTrimmedIsoforms <- TRUE
          }
        }
        
        ### Combine
        isoInfo <- merge(isoInfo, orfInfo, by = 'isoform_id')
      }
      
      if (length(unique(isoInfo$gene_id)) > 1) {
        stop(
          'The isoforms supplied originates from more than one gene - a feature currently not supported. Please revise accordingly'
        )
      }
      geneName <- paste(unique(isoInfo$gene_name), collapse = ',')
      
    }
    
    ### Exon data
    if(TRUE) {
      exonInfo <-
        switchAnalyzeRlist$exons[which(
          switchAnalyzeRlist$exons$isoform_id %in% isoform_id
        ), ]
      if (!length(exonInfo)) {
        stop(
          'It seems like there is no exon information advailable for the gene/transcripts supplied. Please update to the latest version of IsoformSwitchAnalyzeR and try again (re-import the data). If the problem persist contact the developer'
        )
      }
      exonInfoSplit <- split(exonInfo, exonInfo$isoform_id)
      
      chrName <- as.character(seqnames(exonInfo)[1])
    }
    
    ### ORF data
    if(TRUE) {
      if (inclORF) {
        orfStart <-
          lapply(split(isoInfo$orfStartGenomic,  isoInfo$isoform_id),
                 unique)
        orfEnd   <-
          lapply(split(isoInfo$orfEndGenomic, isoInfo$isoform_id),
                 unique)
      } else {
        orfStart <- split(rep(NA, length(isoform_id)), isoform_id)
        orfEnd   <-
          split(rep(NA, length(isoform_id)), isoform_id)
      }
    }
    
    ### Transcript type
    if(TRUE) {
      if (inclCodingPotential) {
        isCoding <-
          lapply(split(isoInfo$codingPotential, isoInfo$isoform_id),
                 unique)
      }
      if (inclPTC) {
        isPTC    <- lapply(split(isoInfo$PTC, isoInfo$isoform_id), unique)
      }
    }
    
    ### Extract basis annoation data
    if(TRUE) {
      ### Make list with annotation
      annotationList <- list()
      
      ### Domain data
      if (inclDomainAnalysis) {
        if (any(
          isoInfo$isoform_id %in%
          switchAnalyzeRlist$domainAnalysis$isoform_id
        )) {
          DomainAnalysis <-
            switchAnalyzeRlist$domainAnalysis[which(
              switchAnalyzeRlist$domainAnalysis$isoform_id %in%
                isoInfo$isoform_id
            ), ]
          DomainAnalysis$isoform_id <-
            factor(DomainAnalysis$isoform_id,
                   levels = unique(isoInfo$isoform_id))
          DomainAnalysis$id <- 1:nrow(DomainAnalysis)
          
          ### Annotate structural variants
          if( !is.null(DomainAnalysis$domain_isotype_simple)) {
            DomainAnalysis$domain_sv <- ifelse(
              DomainAnalysis$domain_isotype_simple == 'Reference',
              yes = '',
              no = ' (Non-ref Isotype)'
            )
            DomainAnalysis$hmm_name <- paste0(
              DomainAnalysis$hmm_name,
              DomainAnalysis$domain_sv
            )
          }
          
          annotationList$protein_domain <- DomainAnalysis[,c('isoform_id','pfamStartGenomic','pfamEndGenomic','hmm_name','id')]
        }
        
      }
      
      if (inclIdrAnalysis) {
        if (any(
          isoInfo$isoform_id %in%
          switchAnalyzeRlist$idrAnalysis$isoform_id
        )) {
          idrAnalysis <-
            switchAnalyzeRlist$idrAnalysis[which(
              switchAnalyzeRlist$idrAnalysis$isoform_id %in%
                isoInfo$isoform_id
            ), ]
          idrAnalysis$isoform_id <-
            factor(idrAnalysis$isoform_id,
                   levels = unique(isoInfo$isoform_id))
          
          #idrAnalysis$idrName <- 'IDR'
          idrAnalysis$id <- 1:nrow(idrAnalysis)
          
          
          annotationList$idr <- idrAnalysis[,c('isoform_id','idrStartGenomic','idrEndGenomic','id')]
          
        }
      }
      
      if (inclSignalPAnalysis) {
        if (any(
          isoInfo$isoform_id %in%
          switchAnalyzeRlist$signalPeptideAnalysis$isoform_id
        )) {
          signalPanalysis <-
            switchAnalyzeRlist$signalPeptideAnalysis[which(
              switchAnalyzeRlist$signalPeptideAnalysis$isoform_id %in%
                isoInfo$isoform_id
            ), ]
          signalPanalysis$genomicClevageAfter <-
            unlist(signalPanalysis$genomicClevageAfter)
          
          ### Signal P region
          signalPdata <- isoInfo[,c('isoform_id','orfStartGenomic')]
          colnames(signalPdata)[2] <- 'spStartGenomic'
          signalPdata$spEndGenomic <- signalPanalysis$genomicClevageAfter[match(
            signalPdata$isoform_id, signalPanalysis$isoform_id
          )]
          signalPdata$id <- signalPdata$isoform_id
          
          signalPanalysis$isoform_id <-
            factor(signalPanalysis$isoform_id,
                   levels = unique(isoInfo$isoform_id))
          
          annotationList$signal_peptide <- signalPdata
        }
      }
      
      if (plotTopology) {
        if (any(
          isoInfo$isoform_id %in%
          switchAnalyzeRlist$topologyAnalysis$isoform_id
        )) {
          tolologyAnalysis <-
            switchAnalyzeRlist$topologyAnalysis[which(
              switchAnalyzeRlist$topologyAnalysis$isoform_id %in%
                isoInfo$isoform_id
            ), ]
          tolologyAnalysis$isoform_id <-
            factor(tolologyAnalysis$isoform_id,
                   levels = unique(isoInfo$isoform_id))
          tolologyAnalysis$id <- 1:nrow(tolologyAnalysis)
          
          annotationList$topology <- tolologyAnalysis[,c('isoform_id','regionStartGenomic','regionEndGenomic','region_type','id')]
        }
        
      }
      
      
    }
    
    ### Domain sites
    if(TRUE) {
      if (inclDomainAnalysis) {
        if (any(
          isoInfo$isoform_id %in%
          switchAnalyzeRlist$domainAnalysis$isoform_id
        )) {
          domainStart <-
            split(
              DomainAnalysis$pfamStartGenomic,
              DomainAnalysis$isoform_id,
              drop = FALSE
            )
          domainEnd   <-
            split(DomainAnalysis$pfamEndGenomic,
                  DomainAnalysis$isoform_id,
                  drop = FALSE)
          domainName  <-
            split(DomainAnalysis$hmm_name,
                  DomainAnalysis$isoform_id,
                  drop = FALSE)
          
          ### Handle if there are any none-unique domain names:
          if (ifMultipleIdenticalAnnotation == 'number') {
            domainName <- lapply(domainName, function(aVec) {
              duplicatedIndex <- which(duplicated(aVec))
              if (length(duplicatedIndex)) {
                aVec[duplicatedIndex] <-
                  paste(aVec[duplicatedIndex],
                        1:length(duplicatedIndex),
                        sep = '.')
              }
              return(aVec)
            })
          } else if (ifMultipleIdenticalAnnotation == 'summarize') {
            domainName <- lapply(domainName, function(aVec) {
              nameTable <- as.data.frame(
                table(aVec),
                stringsAsFactors = FALSE
              )
              
              # only modify those with mutiple instances
              nameTable$newName <- nameTable$aVec
              modifyIndex <- which(nameTable$Freq > 1)
              nameTable$newName[modifyIndex] <-
                paste(nameTable$aVec[modifyIndex],
                      ' (x',
                      nameTable$Freq[modifyIndex],
                      ')',
                      sep = '')
              
              newVec <-
                nameTable$newName[match(aVec, nameTable$aVec)]
            })
          } else if (ifMultipleIdenticalAnnotation != 'ignore') {
            stop(
              'An error occured with ifMultipleIdenticalAnnotation - please contact developers with reproducible example'
            )
          }
        } else {
          domainStart <- NULL # By setting them to null they are ignore from hereon out
          domainEnd   <- NULL
        }
      } else {
        domainStart <- NULL # By setting them to null they are ignore from hereon out
        domainEnd   <- NULL
      }
    }
    
    ### IDR Sites
    if(TRUE) {
      if (inclIdrAnalysis) {
        if(
          any(
            isoInfo$isoform_id %in%
            switchAnalyzeRlist$idrAnalysis$isoform_id
          )
        ) {
          idrStart <-
            split(
              idrAnalysis$idrStartGenomic,
              idrAnalysis$isoform_id,
              drop = FALSE
            )
          idrEnd   <-
            split(idrAnalysis$idrEndGenomic,
                  idrAnalysis$isoform_id,
                  drop = FALSE)
          idrName  <-
            split(idrAnalysis$idr_type,
                  idrAnalysis$isoform_id,
                  drop = FALSE)
          
          ### Handle if there are any none-unique domain names
          if(TRUE) {
            if (ifMultipleIdenticalAnnotation == 'number') {
              idrName <- lapply(idrName, function(aVec) {
                duplicatedIndex <- which(duplicated(aVec))
                if (length(duplicatedIndex)) {
                  aVec[duplicatedIndex] <-
                    paste(aVec[duplicatedIndex],
                          1:length(duplicatedIndex),
                          sep = '.')
                }
                return(aVec)
              })
            } else if (ifMultipleIdenticalAnnotation == 'summarize') {
              idrName <- lapply(idrName, function(aVec) {
                nameTable <- as.data.frame(
                  table(aVec),
                  stringsAsFactors = FALSE
                )
                
                # only modify those with mutiple instances
                nameTable$newName <- nameTable$aVec
                modifyIndex <- which(nameTable$Freq > 1)
                nameTable$newName[modifyIndex] <-
                  paste(nameTable$aVec[modifyIndex],
                        ' (x',
                        nameTable$Freq[modifyIndex],
                        ')',
                        sep = '')
                
                newVec <-
                  nameTable$newName[match(aVec, nameTable$aVec)]
              })
            } else if (ifMultipleIdenticalAnnotation != 'ignore') {
              stop(
                'An error occured with ifMultipleIdenticalAnnotation - please contact developers with reproducible example'
              )
            }
          }
        } else {
          idrStart <- NULL # By setting them to null they are ignore from hereon out
          idrEnd   <- NULL
        }
      } else {
        idrStart <- NULL # By setting them to null they are ignore from hereon out
        idrEnd   <- NULL
      }
    }
    
    ### Signal peptide data
    if(TRUE) {
      if (inclSignalPAnalysis) {
        if (any(
          isoInfo$isoform_id %in%
          switchAnalyzeRlist$signalPeptideAnalysis$isoform_id
        )) {
          cleaveageAfter <-
            split(
              signalPanalysis$genomicClevageAfter,
              signalPanalysis$isoform_id,
              drop = FALSE
            )
        } else {
          cleaveageAfter <-
            NULL # By setting them to null they are ignore from hereon out
        }
      } else {
        cleaveageAfter <-
          NULL # By setting them to null they are ignore from hereon out
      }
    }
    
    ### Toplogy
    if(TRUE) {
      if (plotTopology) {
        if (any(
          isoInfo$isoform_id %in%
          switchAnalyzeRlist$topologyAnalysis$isoform_id
        )) {
          toplologyStart <-
            split(
              tolologyAnalysis$regionStartGenomic,
              tolologyAnalysis$isoform_id,
              drop = FALSE
            )
          toplologyEnd   <-
            split(tolologyAnalysis$regionEndGenomic,
                  tolologyAnalysis$isoform_id,
                  drop = FALSE)
          toplologyType  <-
            split(tolologyAnalysis$region_type,
                  tolologyAnalysis$isoform_id,
                  drop = FALSE)
        } else {
          toplologyStart <- NULL # By setting them to null they are ignore from hereon out
          toplologyEnd   <- NULL
          
          # turn of visualization since no data
          plotTopology <- FALSE
        }
      } else {
        toplologyStart <- NULL # By setting them to null they are ignore from hereon out
        toplologyEnd   <- NULL
        
        # turn of visualization since no data
        plotTopology <- FALSE
        
      }
    }
    
    ### Trimmed sites
    if(TRUE) {
      if( inclTrimmedIsoforms ) {
        trimmedInfo <- orfInfo[,c('isoform_id','trimmedStartGenomic','orfEndGenomic')]
        trimmedInfo$orfEndGenomic[which(is.na(trimmedInfo$trimmedStartGenomic))] <- NA
        trimmedInfo$name <- NA
        
        trimmedStart <-
          split(
            trimmedInfo$trimmedStartGenomic,
            trimmedInfo$isoform_id,
            drop = FALSE
          )
        trimmedEnd   <-
          split(trimmedInfo$orfEndGenomic,
                trimmedInfo$isoform_id,
                drop = FALSE)
      } else {
        trimmedStart <- NULL # By setting them to null they are ignore from hereon out
        trimmedEnd   <- NULL
      }
      
    }
    
  }
  
  
  
  ### Align and massage regions
  if(TRUE) {
    ### Align everything and annoate each region
    if(TRUE) {
      ### Loop over each transcript and make the data.frame with all the annotation data (this is currently the rate limiting step)
      myTranscriptPlotDataList <- list()
      for (i in seq(along.with = exonInfoSplit)) {
        # extract data
        transcriptName <- names(exonInfoSplit)[i]
        localExons     <- exonInfoSplit[[transcriptName]]
        
        # extract local values of where to cut the transcript
        localOrfStart           <-
          orfStart      [[transcriptName]]
        localOrfEnd             <-
          orfEnd        [[transcriptName]] + 1
        localDomainStart        <-
          domainStart   [[transcriptName]]
        localDomainEnd          <-
          domainEnd     [[transcriptName]] + 1
        localIdrStart        <-
          idrStart      [[transcriptName]]
        localIdrEnd          <-
          idrEnd        [[transcriptName]] + 1
        localTopStart        <-
          toplologyStart[[transcriptName]]
        localTopEnd          <-
          toplologyEnd  [[transcriptName]] + 1
        localtrimmedStart    <-
          trimmedStart  [[transcriptName]]
        localtrimmedEnd      <-
          trimmedEnd    [[transcriptName]] + 1
        localPepticeCleaveage   <-
          cleaveageAfter[[transcriptName]]
        
        myCutValues <-
          unique(
            c(
              localOrfStart,
              localOrfEnd,
              localDomainStart,
              localDomainEnd,
              localIdrStart,
              localIdrEnd,
              localTopStart,
              localTopEnd,
              localtrimmedStart,
              localtrimmedEnd,
              localPepticeCleaveage
            )
          ) # NULLs are just removed
        
        # cut the exons into smaller part based on the ORF and domain coordinats (if needed)
        if (length(myCutValues) & !is.na(myCutValues[1])) {
          localExonsDevided <-
            cutGRanges(aGRange = localExons, cutValues = myCutValues)
        } else {
          localExonsDevided <- localExons[, 0]
        }
        
        ### Add annotation
        ## add standard annotation
        localExonsDevided$type <- 'utr'
        localExonsDevided$Domain <- ' transcript'
        localExonsDevided$topology <- 'none'
        
        ## modify if needed
        # ORF
        if(length(localOrfStart)) {
          if (!is.na(localOrfStart)) {
            coordinatPair <- c(localOrfStart, localOrfEnd)
            orfRange <-
              IRanges(min(coordinatPair), max(coordinatPair))
            localExonsDevided$type[queryHits(findOverlaps(
              subject = orfRange,
              query = ranges(localExonsDevided),
              type = 'within'
            ))] <- 'cds'
          }
        }
        
        # domain - loop over each domain
        if (length(localDomainStart)) {
          for (j in 1:length(localDomainStart)) {
            coordinatPair <- c(localDomainStart[j], localDomainEnd[j])
            if( all( !is.na(coordinatPair)) ) {
              domainRange <-
                IRanges(min(coordinatPair), max(coordinatPair))
              localExonsDevided$Domain[queryHits(findOverlaps(
                subject = domainRange,
                query = ranges(localExonsDevided),
                type = 'within'
              ))] <- domainName[[transcriptName]][j]
            }
          }
        }
        
        # IDR
        if (length(localIdrStart)) {
          for (j in 1:length(localIdrStart)) {
            coordinatPair <- c(localIdrStart[j], localIdrEnd[j])
            if( all( !is.na(coordinatPair)) ) {
              domainRange <-
                IRanges(min(coordinatPair), max(coordinatPair))
              localExonsDevided$Domain[queryHits(findOverlaps(
                subject = domainRange,
                query = ranges(localExonsDevided),
                type = 'within'
              ))] <- idrName[[transcriptName]][j]
            }
          }
        }
        
        # signal peptide
        if (length(localPepticeCleaveage)) {
          coordinatPair <- c(localOrfStart, localPepticeCleaveage)
          if( all( !is.na(coordinatPair)) ) {
            peptideRange <-
              IRanges(min(coordinatPair), max(coordinatPair))
            localExonsDevided$Domain[queryHits(findOverlaps(
              subject = peptideRange,
              query = ranges(localExonsDevided),
              type = 'within'
            ))] <- 'Signal Peptide'
          }
        }
        
        # trimmed
        if (length(localtrimmedStart)) {
          for (j in 1:length(localtrimmedStart)) {
            coordinatPair <- c(localtrimmedStart[j], localtrimmedEnd[j])
            if( all( !is.na(coordinatPair)) ) {
              domainRange <-
                IRanges(min(coordinatPair), max(coordinatPair))
              localExonsDevided$Domain[queryHits(findOverlaps(
                subject = domainRange,
                query = ranges(localExonsDevided),
                type = 'within'
              ))] <- 'Not Analyzed'
            }
          }
        }
        
        # topology - loop over each domain
        if (length(localTopStart)) {
          for (j in 1:length(localTopStart)) {
            coordinatPair <- c(localTopStart[j], localTopEnd[j])
            if( all( !is.na(coordinatPair)) ) {
              domainRange <-
                IRanges(min(coordinatPair), max(coordinatPair))
              localExonsDevided$topology[queryHits(findOverlaps(
                subject = domainRange,
                query = ranges(localExonsDevided),
                type = 'within'
              ))] <- toplologyType[[transcriptName]][j]
            }
          }
        }
        
        
        # convert to df
        localExonsDf <- as.data.frame(localExonsDevided)
        
        ### Massage
        localExonsDf$transcript <- transcriptName
        
        ### determine transcript type
        if (inclCodingPotential) {
          localCoding <- isCoding[[transcriptName]]
        } else {
          ### I no prediction default to annotation
          if( ! is.na(localOrfStart) ) {
            localCoding <- TRUE
          } else {
            localCoding <- NULL
          }
        }
        if (inclPTC) {
          localPTC <- isPTC[[transcriptName]]
        } else {
          localPTC <- NULL
        }
        localExonsDf$seqnames <-
          determineTranscriptClass(ptc =  localPTC , coding = localCoding)
        colnames(localExonsDf)[1] <- 'seqnames'
        
        # save
        myTranscriptPlotDataList[[transcriptName]] <- localExonsDf
      }
      myTranscriptPlotData <- do.call(rbind, myTranscriptPlotDataList)
      
    }
    # myTranscriptPlotData
    
    ### Correct names  !!! where Tx Names are changed !!!!
    if (TRUE) {
      ### correct transcript names
      # Create name annotation (this can be omitted when ggplots astetics mapping against type works)
      nameDF <- data.frame(
        oldTxName = unique(myTranscriptPlotData$transcript),
        newTxName = unique(myTranscriptPlotData$transcript),
        stringsAsFactors = FALSE
      )
      
      ### Modify if class code is defined
      if ('class_code' %in% colnames(isoInfo)) {
        nameDF$newTxName <- paste(
          nameDF$newTxName,
          ' (',
          isoInfo$class_code[match(nameDF$oldTxName, isoInfo$isoform_id)],
          ')',
          sep = ''
        )
      }
      
      if ( isConditional ) {
        ### Interpret direction
        isoInfo$direction                                      <- 'Unchanged usage'
        isoInfo$direction[which(isoInfo$dIF > dIFcutoff     )] <- 'Increased usage'
        isoInfo$direction[which(isoInfo$dIF < dIFcutoff * -1)] <- 'Decreased usage'
        
        if( ! optimizeForCombinedPlot ) {
          ### Add dIF
          if( ! all(isoInfo$dIF %in% c(0, Inf, -Inf)) ) {
            isoInfo$direction  <- paste(
              isoInfo$direction,
              ': dIF =',
              formatC(round(isoInfo$dIF, digits = 2),digits = 2, format='f') ,
              sep=' '
            )
          }
          
          ### Add q-values
          if(any( !is.na(isoInfo$isoform_switch_q_value))) {
            if( ! all(isoInfo$isoform_switch_q_value %in% c(1, -Inf)) ) {
              isoInfo$sig <- evalSig(isoInfo$isoform_switch_q_value, alphas = alphas)
              
              isoInfo$direction <- startCapitalLetter(
                paste0(
                  isoInfo$direction,
                  ' (',
                  isoInfo$sig,
                  ')'
                )
              )
            }
          }
        }
        
        ### Make new name
        nameDF$newTxName <- paste(
          nameDF$newTxName,
          '\n(',
          isoInfo$direction[match(nameDF$oldTxName, isoInfo$isoform_id)],
          ')',
          sep = ''
        )
        
      }
      
      
      ### Modify if sub-cell location is defined
      if( 'sub_cell_location' %in% colnames(isoInfo) ) {
        matchVec <- match(nameDF$oldTxName, isoInfo$isoform_id)
        
        if(abbreviateLocations) {
          isoInfo$sub_cell_location <- sapply(
            str_split(isoInfo$sub_cell_location, ','),
            function(x) {
              x <- gsub(pattern = 'Cell_membrane'        , replacement = 'Memb'   , x = x)
              x <- gsub(pattern = 'Cytoplasm'            , replacement = 'Cyto'   , x = x)
              x <- gsub(pattern = 'Endoplasmic_reticulum', replacement = 'ER'     , x = x)
              x <- gsub(pattern = 'Extracellular'        , replacement = 'ExtCell', x = x)
              x <- gsub(pattern = 'Golgi_apparatus'      , replacement = 'Golgi'  , x = x)
              x <- gsub(pattern = 'Lysosome_Vacuole'     , replacement = 'Lys'    , x = x)
              x <- gsub(pattern = 'Mitochondrion'        , replacement = 'Mito'   , x = x)
              x <- gsub(pattern = 'Nucleus'              , replacement = 'Nucl'   , x = x)
              x <- gsub(pattern = 'Peroxisome'           , replacement = 'Perox'  , x = x)
              x <- gsub(pattern = 'Plastid'              , replacement = 'Plastid', x = x)
              x <- paste(x, collapse = ', ')
              return(x)
            }
          )
          
        }
        
        nameDF$newTxName <- paste0(
          nameDF$newTxName,
          '\n(Location: ',
          isoInfo$sub_cell_location[match(
            nameDF$oldTxName,
            isoInfo$isoform_id
          )],
          ')'
        )
      }
      
      ### Modify if  solubility location is defined
      if( 'solubility_status' %in% colnames(isoInfo) ) {
        matchVec <- match(nameDF$oldTxName, isoInfo$isoform_id)
        
        nameDF$newTxName <- paste0(
          nameDF$newTxName,
          '\n(',
          gsub('_',' ', isoInfo$solubility_status[match(
            nameDF$oldTxName,
            isoInfo$isoform_id
          )]),
          ')'
        )
      }
      
      
      
      myTranscriptPlotData$transcript <-
        nameDF$newTxName[match(
          myTranscriptPlotData$transcript, nameDF$oldTxName
        )]
      
      ### Factorize order
      supposedOrder <-
        c('Coding',
          'Non-coding',
          'NMD Insensitive',
          'NMD Sensitive',
          'Transcripts')
      supposedOrder <-
        supposedOrder[which(
          supposedOrder %in% myTranscriptPlotData$seqnames
        )]
      
      myTranscriptPlotData$seqnames <-
        factor(myTranscriptPlotData$seqnames)
      newOrder <-
        match(levels(myTranscriptPlotData$seqnames),
              supposedOrder)
      myTranscriptPlotData$seqnames <-
        factor(myTranscriptPlotData$seqnames,
               levels = levels(myTranscriptPlotData$seqnames)[newOrder])
    }
    
    ### Rescale coordinats if nessesary
    if (rescaleTranscripts) {
      # Might as well work with smaller numbers
      myTranscriptPlotData[, c('start', 'end')] <-
        myTranscriptPlotData[, c('start', 'end')] -
        min(myTranscriptPlotData[, c('start', 'end')]) + 1
      
      ### create conversion table from origial coordiants to rescaled values
      allCoordinats <-
        sort(unique(
          c(
            myTranscriptPlotData$start,
            myTranscriptPlotData$end
          )
        ))
      myCoordinats <-
        data.frame(orgCoordinates = allCoordinats,
                   newCoordinates = allCoordinats)
      for (i in 2:nrow(myCoordinats)) {
        orgDistance <-
          myCoordinats$orgCoordinates[i] -
          myCoordinats$orgCoordinates[i - 1]
        
        newDistance <- max(c(1, round(
          orgDistance^(1/rescaleRoot)
        )))
        
        difference <- orgDistance - newDistance
        
        myCoordinats$newCoordinates [i:nrow(myCoordinats)] <-
          myCoordinats$newCoordinates [i:nrow(myCoordinats)] - difference
      }
      
      # replace original values
      myTranscriptPlotData$start <-
        myCoordinats$newCoordinates[match(
          myTranscriptPlotData$start,
          table = myCoordinats$orgCoordinates
        )]
      myTranscriptPlotData$end <-
        myCoordinats$newCoordinates[match(
          myTranscriptPlotData$end,
          table = myCoordinats$orgCoordinates
        )]
      
    }
    
    ### Revers coordinats if nesseary (and transcript is on minus strand)
    invertCoordinats <-
      reverseMinus & as.character(exonInfo@strand)[1] == '-' # exonInfo
    if (invertCoordinats) {
      # extract min coordinat
      minCoordinat <- min(myTranscriptPlotData[, c('start', 'end')])
      
      # transpose to
      myTranscriptPlotData[, c('start', 'end')] <-
        myTranscriptPlotData[, c('start', 'end')] - minCoordinat + 1
      
      # calculate how much to extract
      subtractNr <- max(myTranscriptPlotData[, c('start', 'end')])
      
      # calculate new coordinats (by subtracting everthing becomes negative, and abs inverts to postive = evertyhing is inversed)
      newCoordinats <-
        abs(myTranscriptPlotData[, c('start', 'end')] - subtractNr)
      
      # overwrite old coordinats
      myTranscriptPlotData$start <- newCoordinats$start
      myTranscriptPlotData$end <- newCoordinats$end
      myTranscriptPlotData$strand <- '+'
      
      # transpose back so min coordiant is the same
      myTranscriptPlotData[, c('start', 'end')] <-
        myTranscriptPlotData[, c('start', 'end')] + minCoordinat
    }
    
    ### order data
    if(TRUE) {
      ### Order by transcript category and name and add index (index ensures names and transcipts are properly plotted)
      myTranscriptPlotData <-
        myTranscriptPlotData[order(myTranscriptPlotData$seqnames,
                                   myTranscriptPlotData$transcript,
                                   decreasing = TRUE), ]
      myTranscriptPlotData$idNr <-
        match(myTranscriptPlotData$transcript ,
              unique(myTranscriptPlotData$transcript))
      
    }
    
  }
  
  ### Create objects for plot parts
  if(TRUE) {
    ### Create lines for toplogy
    if(plotTopology) {
      toplogyInfo <-
        myTranscriptPlotData %>%
        filter(topology %in% c(
          'signal',
          'outside',
          'TMhelix',
          'inside'
        )) %>%
        mutate(
          topGroup = 1:dplyr::n(),
          Topology = dplyr::case_when(
            topology == 'signal'  ~ 'Signal Peptide',
            topology == 'outside' ~ 'Extracellular',
            topology == 'inside'  ~ 'Intracellular',
            topology == 'TMhelix' ~ 'TM helix',
            TRUE ~ 'Something with toplogy assignment went wrong - contact developer'
          ),
          idNr = idNr + rectHegith * codingWidthFactor + rectHegith / 3 #lift them up
        ) %>%
        dplyr::select(idNr, seqnames, start, end, Topology, topGroup) %>%
        tidyr::pivot_longer(cols = c('start','end')) %>%
        as.data.frame()
    }
    
    ### Create coordiants to rectangels
    if (TRUE) {
      # ### Convert coordiants to rectangels coordinats (for plotting) and order them according to draw order
      
      ### calculate rectangle coordinates
      myTranscriptPlotData$ymin <- myTranscriptPlotData$start
      myTranscriptPlotData$ymax <- myTranscriptPlotData$end
      
      myTranscriptPlotData$xmin <-
        myTranscriptPlotData$idNr - rectHegith
      myTranscriptPlotData$xmax <-
        myTranscriptPlotData$idNr + rectHegith
      
      ### Change with of coding regions
      codingIndex <- which(myTranscriptPlotData$type == 'cds')
      myTranscriptPlotData$xmin[codingIndex] <-
        myTranscriptPlotData$idNr[codingIndex] -
        (rectHegith * codingWidthFactor)
      myTranscriptPlotData$xmax[codingIndex] <-
        myTranscriptPlotData$idNr[codingIndex] +
        (rectHegith * codingWidthFactor)
      
      
      ### Change order to reflect annotationImportance
      if(TRUE) {
        ### Create vector with supposed ordering
        annotNameList <- list(
          transcript = ' transcript',
          notAnalyzed = "Not Analyzed"
        )
        if(!is.null(domainStart)) {
          annotNameList$protein_domain <- unique(unlist(domainName))
        }
        if(!is.null(idrStart)) {
          annotNameList$idr <- unique(unlist(idrName))
        }
        if(!is.null(cleaveageAfter)) {
          annotNameList$signal_peptide <- 'Signal Peptide'
        }
        
        annotNameListOrdered <- unlist( annotNameList[c(
          'transcript',
          'notAnalyzed',
          rev(annotationImportance)
        )] )
        
        ### Reorder data
        myTranscriptPlotData$DomainRanking <- match(myTranscriptPlotData$Domain, annotNameListOrdered)
        
        myTranscriptPlotData <- myTranscriptPlotData[order(
          myTranscriptPlotData$DomainRanking,
          myTranscriptPlotData$seqnames,
          myTranscriptPlotData$transcript,
          decreasing = FALSE
        ),]
      }
      
    }
    
    ### Create transcript inton lines with arrows
    if (TRUE) {
      totalLengths <-
        diff(range(
          c(
            myTranscriptPlotData$ymin,
            myTranscriptPlotData$ymax
          )
        ))
      byFactor <- totalLengths / nrArrows
      
      myTranscriptPlotDataSplit <-
        split(myTranscriptPlotData, f = myTranscriptPlotData$idNr)
      arrowlineDataCombined <-
        do.call(rbind, lapply(myTranscriptPlotDataSplit, function(aDF) {
          # aDF <- myTranscriptPlotDataSplit[[1]]
          # extract introns (if coordinats are inverted start have the largest coordinat now)
          if (invertCoordinats) {
            localIntrons <-
              data.frame(gaps(IRanges(
                start = aDF$end, end = aDF$start
              )))
          } else {
            localIntrons <-
              data.frame(gaps(IRanges(
                start = aDF$start, end = aDF$end
              )))
          }
          
          if (nrow(localIntrons)) {
            # Determine munber of arrows based on total intronseize compare to all transcript
            totalIntronSize <- sum(localIntrons$width)
            localNrArrows <-
              totalIntronSize / totalLengths * nrArrows
            
            localIntrons$nrArrows <-
              floor(localIntrons$width /
                      sum(localIntrons$width) * localNrArrows)
            localIntrons$index <-
              seq(along.with = localIntrons$start)
            
            # for each intron make the calculate number of arrows (min 1 arrow pr intron since the seq(+2) will always give two coordiants)
            localArrowlineData <-
              do.call(rbind, lapply(split(
                localIntrons,
                f = seq(along.with = localIntrons$start)
              ), function(localDF) {
                # localDF <- localIntrons[1,]
                mySeq <-
                  seq(
                    min(localDF$start),
                    max(localDF$end) + 1,
                    length.out =  localDF$nrArrows + 2
                  )
                localArrowlineData <-
                  data.frame(
                    seqnames = aDF$seqnames[1],
                    x = aDF$idNr[1],
                    y = mySeq[-length(mySeq)] - 1,
                    yend = mySeq[-1],
                    nrArrows = localDF$nrArrows
                  )
                return(localArrowlineData)
              }))
            
            # reverse arrow direction if transcript is on minus strand (and reverse is off) - nessesary since calculations are done on + strand since IRanges cannot handle negative widths
            if (!reverseMinus &
                as.character(exonInfo@strand)[1] == '-') {
              localArrowlineData <-
                data.frame(
                  seqnames = localArrowlineData$seqnames,
                  x = localArrowlineData$x,
                  y = localArrowlineData$yend,
                  yend = localArrowlineData$y,
                  nrArrows = localArrowlineData$nrArrows
                )
            }
            
            return(localArrowlineData)
            
          } else {
            localArrowlineData <-
              data.frame(
                seqnames = aDF$seqnames[1],
                x = aDF$idNr[1],
                y = NA,
                yend = NA,
                nrArrows = c(0, 1)
              )
            return(localArrowlineData)
          }
          
          
        }))
      
      arrowlineDataArrows <-
        arrowlineDataCombined[which(arrowlineDataCombined$nrArrows != 0), ]
      arrowlineDataLines <-
        arrowlineDataCombined[which(arrowlineDataCombined$nrArrows == 0), ]
    }
    
    ### create data.frame for converting between index id and transcript name
    idData <-
      unique(myTranscriptPlotData[, c('seqnames', 'transcript', 'idNr')])
    
    ### create color code for domains
    if(TRUE) {
      domainsFound <-
        unique(myTranscriptPlotData$Domain [which(
          ! myTranscriptPlotData$Domain %in% c(' transcript', "Not Analyzed")
        )])
      
      ### Reorder
      domainsFound <- sort(domainsFound)
      domainsFound <- domainsFound[c(
        which( domainsFound == "Signal Peptide"),
        which( domainsFound != "Signal Peptide")
      )]
      
      ### Get colors
      if(TRUE) {
        if (length(domainsFound) == 0) {
          domainsColor <- NULL
        } else if (length(domainsFound) < 3) {
          #domainsColor <-
          #    RColorBrewer::brewer.pal(
          #        n = 3,
          #        name = 'Dark2'
          #    )[2:(length(domainsFound) + 1)]
          
          domainsColor <-
            RColorBrewer::brewer.pal(
              n = 12,
              name = 'Paired'
            )[c(2,4,6)[1:length(domainsFound)]]
          
        } else if (length(domainsFound) > 12) {
          gg_color_hue <- function(n) {
            hues <- seq(15, 375, length = n + 1)
            hcl(h = hues,
                l = 65,
                c = 100)[1:n]
          }
          domainsColor <- gg_color_hue(length(domainsFound))
        } else {
          domainsColor <-
            RColorBrewer::brewer.pal(n = length(domainsFound), name = 'Paired')
        }
      }
      
      ### Fix order
      if( "Not Analyzed" %in% myTranscriptPlotData$Domain ) {
        domainsFound <- c(domainsFound, "Not Analyzed")
        
        correspondingColors <- c(
          "#161616", # for ' transcript'
          domainsColor,
          '#595959' # For "not analyzed"
        )
        
        ### Move "not analyzed" color to its corresponding position
        apperanceInData <- sort(unique(myTranscriptPlotData$Domain))
        
        moveToIndex <- which(sort(apperanceInData) == "Not Analyzed")
        correspondingColors <- c(
          correspondingColors[1:(moveToIndex-1)],
          correspondingColors[length(correspondingColors)],
          correspondingColors[(moveToIndex):(length(correspondingColors)-1)]
        )
      } else {
        correspondingColors <- c(
          "#161616", # for ' transcript'
          domainsColor
        )
      }
      
    }
    
    ### Combined plotting
    if (optimizeForCombinedPlot) {
      
      ### Collect all labels
      allLabels <- c(
        condition1,
        condition2,
        unique(domainsFound)
      )
      if(plotTopology) {
        allLabels <- c(
          allLabels,
          unique(toplogyInfo$Topology)
        )
      }
      
      ### Find length of longest label
      analyzeStrandCompositionInWhiteSpaces <-
        function(aString) {
          # aString <- 'Ab1'
          
          round(sum(sapply(
            strsplit(aString, '')[[1]],
            function(aCharacter) {
              # Test if whitespace
              if (aCharacter == ' ') {
                return(1)
              }
              # Test if number
              if (!is.na(suppressWarnings(as.integer(aCharacter)))) {
                return(2) # whitespace pr number
              }
              # test symbols
              if (aCharacter == '_') {
                return(2) # whitespace pr character
              }
              # test symbols
              if (aCharacter == '.') {
                return(1) # whitespace pr numbers
              }
              # test upper
              if (aCharacter == toupper(aCharacter)) {
                return(2.4) # whitespace pr uppercase
              }
              # else it is probably lower
              return(1.8) # whitespace pr lowercase
              
            })))
          
        }
      
      maxCharacter <-
        max(c(
          sapply(allLabels, analyzeStrandCompositionInWhiteSpaces) ,
          50
        ))
      
      ### Modify names to match length
      modifyNames <- function(aVec, extendToLength) {
        tmp <- sapply(aVec, function(x) {
          currentLength <- analyzeStrandCompositionInWhiteSpaces(x)
          whitespacesToAdd <-
            round(extendToLength - currentLength - 1)
          if (whitespacesToAdd > 0) {
            x <-
              paste(x, paste(rep(
                ' ', whitespacesToAdd
              ), collapse = ''), collapse = '')
          }
          return(x)
        })
        names(tmp) <- NULL
        
        return(tmp)
      }
      
      
      ### Modify them
      domainsFound <- modifyNames(domainsFound, maxCharacter)
      
      modifiedNames <-
        data.frame(
          org = allLabels,
          new = modifyNames(allLabels, maxCharacter),
          stringsAsFactors = FALSE
        )
      indexToModify <-
        which(myTranscriptPlotData$Domain != ' transcript')
      myTranscriptPlotData$Domain[indexToModify] <-
        modifiedNames$new[match(
          myTranscriptPlotData$Domain[indexToModify] , modifiedNames$org
        )]
      if(plotTopology) {
        toplogyInfo$Topology <-
          modifiedNames$new[match(
            toplogyInfo$Topology, modifiedNames$org
          )]
        
      }
      
    }
    
    ### factorize seqnames for all 3 datasets to ensure correct facetting
    if(TRUE) {
      myTranscriptPlotData$seqnames <-
        factor(myTranscriptPlotData$seqnames, levels = supposedOrder)
      arrowlineDataArrows$seqnames  <-
        factor(arrowlineDataArrows$seqnames,  levels = supposedOrder)
      arrowlineDataLines$seqnames   <-
        factor(arrowlineDataLines$seqnames,   levels = supposedOrder)
      idData$seqnames               <-
        factor(idData$seqnames,               levels = supposedOrder)
      
      if(plotTopology) {
        toplogyInfo$seqnames               <-
          factor(toplogyInfo$seqnames,               levels = supposedOrder)
      }
    }
    
  }
  
  
  ### Build plot

    
    myPlot <- ggplot()
    
    ### Introns
    if (nrow(arrowlineDataLines)) {
      myPlot <-
        myPlot + geom_segment(data = arrowlineDataLines, aes(
          x = y,
          xend = yend,
          y = x,
          yend = x
        ))
      
    }
    if (nrow(arrowlineDataArrows)) {
      myPlot <-
        myPlot + geom_segment(
          data = arrowlineDataArrows,
          aes(
            x = y,
            xend = yend,
            y = x,
            yend = x
          ),
          arrow = arrow(length = unit(arrowSize, "cm"))
        )  # intron arrows
    }

###Manually change 5' and 3' UTRs
myTranscriptPlotData["PB.1207.309.1", 6] <- "5' UTR"
myTranscriptPlotData["PB.1207.309.33", 6] <- "3' UTR"

myTranscriptPlotData["PB.1207.317.1", 6] <- "5' UTR"
myTranscriptPlotData["PB.1207.317.33", 6] <- "3' UTR"

myTranscriptPlotData["PB.1207.336.1", 6] <- "5' UTR"
myTranscriptPlotData["PB.1207.336.2", 6] <- "5' UTR"
myTranscriptPlotData["PB.1207.336.3", 6] <- "5' UTR"
myTranscriptPlotData["PB.1207.336.30", 6] <- "3' UTR"

myTranscriptPlotData[which(grepl("cds", myTranscriptPlotData$type)), 6] <- "CDS"

myTranscriptPlotData[which(grepl("ArgoN", myTranscriptPlotData$Domain)), 6] <- "ArgoN"

myTranscriptPlotData$type <- factor(myTranscriptPlotData$type, levels=c("5' UTR", 'ArgoN', 'CDS', "3' UTR"))


    ### transcripts
      myPlot <- myPlot +
        geom_rect(
          data = myTranscriptPlotData,
          aes(
            xmin = ymin,
            ymin = xmin,
            xmax = ymax,
            ymax = xmax,
            fill = type
          )
        )
      
      ### Transcript color and theme
      myPlot <- myPlot +
        scale_fill_manual(
          values = c("#E25186", "#87CEFF", "#161616", "#91CE91")
        ) + # Correct domian color code so transcripts are black and not shown
        scale_y_continuous(breaks = idData$idNr, labels = idData$transcript) # change index numbers back to names
      
      
  
    
    ### Optimize apperance
    if(TRUE) {
      ### style and themes
      myPlot <-
        myPlot +
        localTheme + theme(strip.text.y = element_text(angle = 0)) + # change theme and rotate facette labes (ensures readability even though frew are pressent)
        theme(axis.title.y = element_blank()) + # remove y-axis label
        theme(strip.background = element_rect(fill = "white", linewidth = 0.5))
      
      # facette against transcript type if nessesary  <- CHANGING THIS!!!!
      if (all(supposedOrder == 'Transcripts')) {
        myPlot <-
          myPlot + facet_grid(seqnames ~ ., scales = 'free_y', space = 'free')
      }
      
      # Modify X axis
      if (!plotXaxis) {
        # remove axis if rescaled
        myPlot <-
          myPlot + theme(
            axis.line = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_blank()
          )
      } else {
        # add chr name
        myPlot <- myPlot + labs(x = chrName)
      }
      
      if( ! optimizeForCombinedPlot ) {
        if( isConditional ) {
          if( any(grepl('^placeholder_1$|^placeholder_2$', c(condition1, condition2)) ) ) {
            myPlot <- myPlot + labs(title = paste(
              'The isoform switch in',
              geneName,
              sep = ' '
            ))
          } else {
            myPlot <- myPlot + labs(title = paste(
              'The isoform switch in',
              geneName,
              paste0(' (', condition1, ' vs ', condition2,')'),
              sep = ' '
            ))
          }
        } else {
          myPlot <- myPlot + labs(title = paste(
            'The isoforms in',
            geneName,
            sep = ' '
          ))
        }
      }
    }
    
return(myPlot)

}
