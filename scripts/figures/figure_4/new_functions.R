# MODIFIED FUNCTIONS FROM ISOFORMSWITCH ANALYZER PACKAGE 
### For analyzing consequences
analyzeSwitchConsequences_new <- function(
    switchAnalyzeRlist,
    consequencesToAnalyze = c(
      'intron_retention',
      'coding_potential',
      'ORF_seq_similarity',
      'NMD_status',
      'domains_identified',
      'domain_isotype',
      'IDR_identified',
      'IDR_type',
      'signal_peptide_identified'
    ),
    alpha = 0.05,
    dIFcutoff = 0.1,
    onlySigIsoforms = FALSE,
    ntCutoff = 50,
    ntFracCutoff = NULL,
    ntJCsimCutoff = 0.8,
    AaCutoff = 10,
    AaFracCutoff = 0.8,
    AaJCsimCutoff = 0.9,
    removeNonConseqSwitches = TRUE,
    showProgress = TRUE,
    quiet = FALSE
) {
  ### Check input
  if (TRUE) {
    # check switchAnalyzeRlist
    if (class(switchAnalyzeRlist) != 'switchAnalyzeRlist') {
      stop(
        'The object supplied to \'switchAnalyzeRlist\' is not a \'switchAnalyzeRlist\''
      )
    }
    
    if (alpha < 0 |
        alpha > 1) {
      warning('The alpha parameter should usually be between 0 and 1 ([0,1]).')
    }
    if (alpha > 0.05) {
      warning(
        'Most journals and scientists consider an alpha larger than 0.05 untrustworthy. We therefore recommend using alpha values smaller than or queal to 0.05'
      )
    }
    
    # test wether switching have been analyzed
    if (!any(!is.na(
      switchAnalyzeRlist$isoformFeatures$gene_switch_q_value
    ))) {
      stop(
        'The analsis of isoform switching must be performed before functional consequences can be analyzed. Please run any of the isoformSwitchTest*() functions and try again.'
      )
    }
    
    acceptedTypes <- c(
      # Transcript
      'tss',
      'tts',
      'last_exon',
      'isoform_length',
      'exon_number',
      'intron_structure',
      'intron_retention',
      'isoform_class_code',
      # cpat
      'coding_potential',
      # ORF
      'ORF_genomic',
      'ORF_length',
      '5_utr_length',
      '3_utr_length',
      # seq similarity
      'isoform_seq_similarity',
      'ORF_seq_similarity',
      '5_utr_seq_similarity',
      '3_utr_seq_similarity',
      # ORF
      'NMD_status',
      # pfam
      'domains_identified',
      'genomic_domain_position',
      'domain_length',
      'domain_isotype',
      
      # SignalIP
      'signal_peptide_identified',
      
      # IDR
      'IDR_identified',
      'IDR_length',
      'IDR_type',
      
      # sub cell
      'sub_cell_location',
      'sub_cell_shift_to_cell_membrane',
      'sub_cell_shift_to_cytoplasm',
      'sub_cell_shift_to_nucleus',
      'sub_cell_shift_to_Extracellular',
      
      # topology
      'isoform_topology',
      'extracellular_region_count',
      'intracellular_region_count',
      'extracellular_region_length',
      'intracellular_region_length'
    )
    
    if (!all(consequencesToAnalyze %in% c('all', acceptedTypes))) {
      stop(
        paste(
          'The argument(s) supplied to \'typeOfconsequence\' are not accepted.',
          'Please see ?analyzeSwitchConsequences for description of which strings are allowed.',
          'The problem is:',
          paste(setdiff(
            consequencesToAnalyze , c('all', acceptedTypes)
          ), collapse = ', '),
          sep = ' '
        )
      )
    }
    
    if ('all' %in% consequencesToAnalyze) {
      consequencesToAnalyze <- acceptedTypes
    }
    
    ## Test whether annotation is advailable
    if ('intron_retention'  %in% consequencesToAnalyze) {
      if (is.null(switchAnalyzeRlist$intronRetentionAnalysis) & is.null( switchAnalyzeRlist$AlternativeSplicingAnalysis)) {
        stop(
          'To test for intron retention alternative splicing must first be classified. Please run analyzeIntronRetention() and try again.'
        )
      }
    }
    
    if (grepl('cufflinks', switchAnalyzeRlist$sourceId)) {
      if ('isoform_class_code'  %in% consequencesToAnalyze) {
        if (!'class_code' %in%
            colnames(switchAnalyzeRlist$isoformFeatures)
        ) {
          stop(
            'The switchAnalyzeRlist does not contail the calss_code information'
          )
        }
      }
    } else {
      consequencesToAnalyze <-
        setdiff(consequencesToAnalyze, 'isoform_class_code')
    }
    
    if (any(
      c(
        'ORF_length',
        '5_utr_length',
        '3_utr_length',
        'ORF_seq_similarity',
        '5_utr_seq_similarity',
        '3_utr_seq_similarity',
        'domains_identified',
        'genomic_domain_position',
        'domain_length',
        'signal_peptide_identified',
        'IDR_identified',
        'IDR_length',
        'IDR_type',
        'sub_cell_location'
      ) %in% consequencesToAnalyze
    )) {
      if (!'PTC' %in% colnames(switchAnalyzeRlist$isoformFeatures)) {
        stop(
          'To test differences in ORF or any annotation derived from these, ORF must be annotated. Please run \'addORFfromGTF()\' (and if nessesary \'analyzeNovelIsoformORF()\') and try again'
        )
        
      }
      
      if( ! is.null(switchAnalyzeRlist$orfAnalysis$orf_origin ) ) {
        if ( any( switchAnalyzeRlist$orfAnalysis$orf_origin == 'not_annotated_yet' )) {
          stop('Some ORFs have not been annotated yet. Please return to the analyzeNovelIsoformORF() step and start again.')
        }
      }
    }
    if ('coding_potential'  %in% consequencesToAnalyze) {
      if (
        !'codingPotential' %in%
        colnames(switchAnalyzeRlist$isoformFeatures))
      {
        stop(
          'To test differences in coding_potential, the result of the CPAT analysis must be advailable. Please run analyzeCPAT() or analyzeCPC2 and try again.'
        )
      }
    }
    if (any(
      c(
        'domains_identified',
        'domain_length',
        'genomic_domain_position'
      )  %in% consequencesToAnalyze
    )) {
      if (is.null(switchAnalyzeRlist$domainAnalysis)) {
        stop(
          'To test differences in protein domains, the result of the Pfam analysis must be advailable. Please run analyzePFAM() and try again.'
        )
      }
    }
    if ('signal_peptide_identified'  %in% consequencesToAnalyze) {
      if (is.null(switchAnalyzeRlist$signalPeptideAnalysis)) {
        stop(
          'To test differences in signal peptides, the result of the SignalP analysis must be advailable. Please run analyzeSignalP() and try again.'
        )
      }
    }
    if ( any(c('IDR_identified','IDR_type')  %in% consequencesToAnalyze)) {
      if (is.null(switchAnalyzeRlist$idrAnalysis)) {
        stop(
          'To test differences in IDR, the result of the NetSurfP2 analysis must be advailable. Please run analyzeNetSurfP2() and try again,'
        )
      }
    }
    if( 'IDR_type' %in% consequencesToAnalyze ) {
      if( ! 'idr_type' %in% colnames(switchAnalyzeRlist$idrAnalysis) ) {
        stop('To analyse IDR_type the IDR analysis must have been done using IUPred2A and imported with the analyzeIUPred2A() function.')
      }
    }
    
    if (!is.numeric(ntCutoff)) {
      stop('The \'ntCutoff\' arugment must be an numeric')
    }
    if (ntCutoff <= 0) {
      stop('The \'ntCutoff\' arugment must be an numeric > 0')
    }
    
    if (!is.null(ntFracCutoff)) {
      if (ntFracCutoff <= 0 | ntFracCutoff > 1) {
        stop(
          'The \'ntFracCutoff\' arugment must be a numeric in the interval (0,1]. Use NULL to disable.'
        )
      }
    }
    
    ### test sequence annotation
    if (any(
      consequencesToAnalyze %in% c(
        'isoform_seq_similarity',
        '5_utr_seq_similarity',
        '3_utr_seq_similarity'
      )
    )) {
      if (!any(names(switchAnalyzeRlist) == 'ntSequence')) {
        stop(
          'The transcrip nucleotide sequences must be added to the switchAnalyzeRlist before overlap analysis can be performed. These can be added by using the \'extractSequence()\' function.'
        )
      }
    }
    if (any(consequencesToAnalyze %in% c('ORF_seq_similarity'))) {
      if (!any(names(switchAnalyzeRlist) == 'aaSequence')) {
        stop(
          'The transcrip ORF amino acid sequences must be added to the switchAnalyzeRlist before ORF overlap analysis can be performed. These can be added by using the \'extractSequence()\' function.'
        )
      }
    }
    
    if ('sub_cell_location'  %in% consequencesToAnalyze) {
      if ( ! 'sub_cell_location' %in% colnames(switchAnalyzeRlist$isoformFeatures)) {
        stop(
          'Cannot test for differences in sub-cellular location as such results are not annotated. Run analyzeDeepLoc2() and try again.'
        )
      }
    }
    
    if ('isoform_topology'  %in% consequencesToAnalyze) {
      if ( ! 'topologyAnalysis' %in% names(switchAnalyzeRlist) ) {
        stop(
          'Cannot test for differences in topology as such results are not annotated. Run analyzeDeepTMHMM() and try again.'
        )
      }
    }
    
  }
  
  if (showProgress & !quiet) {
    progressBar <- 'text'
  } else {
    progressBar <- 'none'
  }
  
  ### Subset to relevant data
  if (TRUE) {
    if (!quiet) {
      message('Step 1 of 4: Extracting genes with isoform switches...')
    }
    
    ### Extract Iso pairs
    pairwiseIsoComparison <- extractSwitchPairs(
      switchAnalyzeRlist = switchAnalyzeRlist,
      alpha = alpha,
      dIFcutoff = dIFcutoff,
      onlySigIsoforms = onlySigIsoforms
    )
    
    ### Extract isoform_id pairs
    pairwiseIsoComparisonUniq <-
      unique(pairwiseIsoComparison[,c(
        'isoformUpregulated', 'isoformDownregulated'
      )])
    pairwiseIsoComparisonUniq$comparison <-
      1:nrow(pairwiseIsoComparisonUniq)
    
    ### Generate size reduced switchAnalyzeRList
    minimumSwitchList <- makeMinimumSwitchList(
      orgSwitchList = switchAnalyzeRlist,
      isoformsToKeep = unique(
        c(
          pairwiseIsoComparisonUniq$isoformUpregulated,
          pairwiseIsoComparisonUniq$isoformDownregulated
        )
      ))
    
    
  }
  
  ### Loop over all the the resulting genes and do a all pairwise comparison between up and down.
  if (TRUE) {
    if (!quiet) {
      message(
        paste(
          'Step 2 of 4: Analyzing',
          nrow(pairwiseIsoComparisonUniq),
          'pairwise isoforms comparisons...',
          sep = ' '
        )
      )
    }
    
    consequencesOfIsoformSwitching <- plyr::dlply(
      .data = pairwiseIsoComparisonUniq,
      .variables = 'comparison',
      .parallel = FALSE,
      .inform = TRUE,
      .progress = progressBar,
      .fun = function(aDF) {
        # aDF <- pairwiseIsoComparisonUniq[1,]
        compareAnnotationOfTwoIsoforms_new(
          switchAnalyzeRlist    = minimumSwitchList,
          consequencesToAnalyze = consequencesToAnalyze,
          upIso                 = aDF$isoformUpregulated,
          downIso               = aDF$isoformDownregulated,
          ntCutoff              = ntCutoff,
          ntFracCutoff          = ntFracCutoff,
          ntJCsimCutoff         = ntJCsimCutoff,
          AaCutoff              = AaCutoff,
          AaFracCutoff          = AaFracCutoff,
          AaJCsimCutoff         = AaJCsimCutoff,
          addDescription        = TRUE,
          testInput             = FALSE # already done by this function
        )
      }
    )
    
    ### Remove to those instances where there where no consequences
    if (removeNonConseqSwitches) {
      ## Remove to those instances where there where no consequences
      consequencesOfIsoformSwitching <-
        consequencesOfIsoformSwitching[which(
          sapply(consequencesOfIsoformSwitching, function(aDF) {
            any(aDF$isoformsDifferent)
          })
        )]
      if (!length(consequencesOfIsoformSwitching)) {
        stop('No isoform switches with the analyzed consequences were found.')
      }
    }
  }
  
  ### Massage result
  if (TRUE) {
    if (!quiet) {
      message(paste(
        'Step 3 of 4: Massaging isoforms comparisons results...',
        sep = ' '
      ))
    }
    ### Convert from list to df
    consequencesOfIsoformSwitchingDf <-
      myListToDf(consequencesOfIsoformSwitching)
    
    ### Add the origin info
    consequencesOfIsoformSwitchingDfcomplete <- dplyr::inner_join(
      consequencesOfIsoformSwitchingDf,
      pairwiseIsoComparison,
      by = c('isoformUpregulated', 'isoformDownregulated')
    )
    #### Add additional information
    #consequencesOfIsoformSwitchingDfcomplete <- dplyr::inner_join(
    #    consequencesOfIsoformSwitchingDfcomplete,
    #    switchAnalyzeRlist$isoformFeatures[match(
    #        unique(consequencesOfIsoformSwitchingDfcomplete$gene_ref),
    #        switchAnalyzeRlist$isoformFeatures$gene_ref
    #    ),
    #    c('gene_ref',
    #      'gene_id',
    #      'gene_name',
    #      'condition_1',
    #      'condition_2')],
    #    by = 'gene_ref'
    #)
    
    ### reorder
    newOrder <- na.omit(match(
      c(
        'gene_ref',
        'gene_id',
        'gene_name',
        'condition_1',
        'condition_2',
        'isoformUpregulated',
        'isoformDownregulated',
        'iso_ref_up',
        'iso_ref_down',
        'featureCompared',
        'isoformsDifferent',
        'switchConsequence'
      ),
      colnames(consequencesOfIsoformSwitchingDfcomplete)
    ))
    consequencesOfIsoformSwitchingDfcomplete <-
      consequencesOfIsoformSwitchingDfcomplete[, newOrder]
    
    consequencesOfIsoformSwitchingDfcomplete <-
      consequencesOfIsoformSwitchingDfcomplete[order(
        consequencesOfIsoformSwitchingDfcomplete$gene_ref,
        consequencesOfIsoformSwitchingDfcomplete$isoformUpregulated,
        consequencesOfIsoformSwitchingDfcomplete$isoformDownregulated
      ), ]
    
  }
  
  ### Add result to switchAnalyzeRlist
  if (TRUE) {
    if (!quiet) {
      message('Step 4 of 4: Preparing output...')
    }
    ### Add full analysis
    switchAnalyzeRlist$switchConsequence <-
      consequencesOfIsoformSwitchingDfcomplete
    
    # extract indexes of those analyzed
    indexesAnalyzed <-
      which(
        switchAnalyzeRlist$isoformFeatures$gene_ref %in%
          pairwiseIsoComparison$gene_ref
      )
    
    # Add indicator of switch consequence to those analyzed
    switchAnalyzeRlist$isoformFeatures$switchConsequencesGene <- NA
    
    switchAnalyzeRlist$isoformFeatures$switchConsequencesGene[
      indexesAnalyzed
    ] <- switchAnalyzeRlist$isoformFeatures$gene_ref[indexesAnalyzed] %in%
      consequencesOfIsoformSwitchingDfcomplete$gene_ref[which(
        consequencesOfIsoformSwitchingDfcomplete$isoformsDifferent
      )]
  }
  
  if (!quiet) {
    totalNrGene <-
      extractSwitchSummary(
        switchAnalyzeRlist,
        filterForConsequences = TRUE,
        includeCombined = TRUE
      )
    totalNrGene <- totalNrGene$nrGenes[which(
      totalNrGene$Comparison == 'combined'
    )]
    message(
      paste(
        'Identified',
        totalNrGene,
        'genes with containing isoforms switching with functional consequences...',
        sep = ' '
      )
    )
  }
  return(switchAnalyzeRlist)
}

compareAnnotationOfTwoIsoforms_new <- function (switchAnalyzeRlist, consequencesToAnalyze = "all", 
          upIso, downIso, addDescription = TRUE, onlyRepportDifferent = FALSE, 
          ntCutoff = 50, ntFracCutoff = NULL, ntJCsimCutoff = 0.8, 
          AaCutoff = 10, AaFracCutoff = 0.8, AaJCsimCutoff = 0.9, 
          testInput = TRUE) 
{
  if (testInput) {
    if (class(switchAnalyzeRlist) != "switchAnalyzeRlist") {
      stop("The object supplied to 'switchAnalyzeRlist' is not a 'switchAnalyzeRlist'")
    }
    if (!any(!is.na(switchAnalyzeRlist$isoformFeatures$gene_switch_q_value))) {
      stop("The analsis of isoform switching must be performed before functional consequences can be analyzed. Please run ?detectIsoformSwitching and try again.")
    }
    acceptedTypes <- c("tss", "tts", "last_exon", "isoform_length", 
                       "exon_number", "intron_structure", "intron_retention", 
                       "isoform_class_code", "coding_potential", "ORF_genomic", 
                       "ORF_length", "5_utr_length", "3_utr_length", "isoform_seq_similarity", 
                       "ORF_seq_similarity", "5_utr_seq_similarity", "3_utr_seq_similarity", 
                       "NMD_status", "domains_identified", "genomic_domain_position", 
                       "domain_length", "domain_isotype", "signal_peptide_identified", 
                       "IDR_identified", "IDR_length", "IDR_type", "sub_cell_location", 
                       "sub_cell_shift_to_cell_membrane", "sub_cell_shift_to_cytoplasm", 
                       "sub_cell_shift_to_nucleus", "sub_cell_shift_to_Extracellular", 
                       "isoform_topology", "extracellular_region_count", 
                       "intracellular_region_count", "extracellular_region_length", 
                       "intracellular_region_length")
    if (!all(consequencesToAnalyze %in% c("all", acceptedTypes))) {
      stop(paste("The argument(s) supplied to 'typeOfconsequence' are not accepted.", 
                 "Please see ?analyzeSwitchConsequences for description of which strings are allowed.", 
                 "The problem is:", paste(setdiff(consequencesToAnalyze, 
                                                  c("all", acceptedTypes)), collapse = ", "), 
                 sep = " "))
    }
    if ("all" %in% consequencesToAnalyze) {
      consequencesToAnalyze <- acceptedTypes
    }
    if ("intron_retention" %in% consequencesToAnalyze) {
      if (is.null(switchAnalyzeRlist$intronRetentionAnalysis) & 
          is.null(switchAnalyzeRlist$AlternativeSplicingAnalysis)) {
        stop("To test for intron retention alternative splicing must first be classified. Please run analyzeIntronRetention() and try again.")
      }
    }
    if (grepl("cufflinks", switchAnalyzeRlist$sourceId)) {
      if ("isoform_class_code" %in% consequencesToAnalyze) {
        if (!"class_code" %in% colnames(switchAnalyzeRlist$isoformFeatures)) {
          stop("The switchAnalyzeRlist does not contail the calss_code information")
        }
      }
    }
    else {
      consequencesToAnalyze <- setdiff(consequencesToAnalyze, 
                                       "isoform_class_code")
    }
    if (any(c("ORF_genomic", "ORF_length", "NMD_status") %in% 
            consequencesToAnalyze)) {
      if (!"PTC" %in% colnames(switchAnalyzeRlist$isoformFeatures)) {
        stop("To test differences in ORFs or PCT, ORF must be annotated. Please run 'addORFfromGTF()' (and if nessesary 'analyzeNovelIsoformORF()') and try again")
      }
    }
    if ("coding_potential" %in% consequencesToAnalyze) {
      if (!"codingPotential" %in% colnames(switchAnalyzeRlist$isoformFeatures)) {
        stop("To test differences in coding_potential, the result of the CPAT analysis must be advailable. Please run addCPATanalysis() and try again")
      }
    }
    if (any(c("domains_identified", "domain_length", "genomic_domain_position") %in% 
            consequencesToAnalyze)) {
      if (is.null(switchAnalyzeRlist$domainAnalysis)) {
        stop("To test differences in protein domains, the result of the Pfam analysis must be advailable. Please run addPFAManalysis() and try again")
      }
    }
    if ("signal_peptide_identified" %in% consequencesToAnalyze) {
      if (is.null(switchAnalyzeRlist$signalPeptideAnalysis)) {
        stop("To test differences in signal peptides, the result of the SignalP analysis must be advailable. Please run addSignalIPanalysis() and try again")
      }
    }
    if (any(c("sub_cell_location", "sub_cell_shift_to_cell_membrane", 
              "sub_cell_shift_to_cytoplasm", "sub_cell_shift_to_nucleus", 
              "sub_cell_shift_to_Extracellular") %in% consequencesToAnalyze)) {
      if (!"sub_cell_location" %in% colnames(switchAnalyzeRlist$isoformFeatures)) {
        stop("Cannot test for differences in sub-cellular location as such results are not annotated")
      }
    }
    if (any(c("isoform_topology", "intracellular_region_length", 
              "extracellular_region_length", "extracellular_region_count", 
              "intracellular_region_count") %in% consequencesToAnalyze) %in% 
        consequencesToAnalyze) {
      if (!"topologyAnalysis" %in% names(switchAnalyzeRlist)) {
        stop("Cannot test for differences in topology as such results are not annotated. Run analyzeDeepTMHMM() and try again.")
      }
    }
    if (!is.numeric(ntCutoff)) {
      stop("The ntCutoff arugment must be an numeric")
    }
    if (ntCutoff <= 0) {
      stop("The ntCutoff arugment must be an numeric > 0")
    }
    if (!is.null(ntFracCutoff)) {
      if (ntFracCutoff <= 0 | ntFracCutoff > 1) {
        stop("The ntFracCutoff arugment must be a numeric in the interval (0,1]. Use NULL to disable.")
      }
    }
    if (any(consequencesToAnalyze %in% c("isoform_seq_similarity", 
                                         "5_utr_seq_similarity", "3_utr_seq_similarity"))) {
      if (!any(names(switchAnalyzeRlist) == "ntSequence")) {
        stop("The transcrip nucleotide sequences must be added to the switchAnalyzeRlist before overlap analysis can be performed. Please run 'extractSequence' and try again.")
      }
    }
    if (any(consequencesToAnalyze %in% c("ORF_seq_similarity"))) {
      if (!any(names(switchAnalyzeRlist) == "aaSequence")) {
        stop("The transcrip ORF amino acid sequences must be added to the switchAnalyzeRlist before ORF overlap analysis can be performed. Please run 'extractSequence' and try again.")
      }
    }
  }
  if (TRUE) {
    if (!is.null(ntFracCutoff)) {
      fractionFilter <- TRUE
    }
    else {
      fractionFilter <- FALSE
    }
    isoformsToAnalyze <- c(upIso, downIso)
    names(isoformsToAnalyze) <- c("up", "down")
    exonData <- switchAnalyzeRlist$exons[which(switchAnalyzeRlist$exons$isoform_id %in% 
                                                 isoformsToAnalyze), "isoform_id"]
    exonDataList <- split(exonData, f = exonData$isoform_id)
    columnsToExtract <- "isoform_id"
    if ("isoform_class_code" %in% consequencesToAnalyze) {
      columnsToExtract <- c(columnsToExtract, "class_code")
    }
    if ("intron_retention" %in% consequencesToAnalyze) {
      columnsToExtract <- c(columnsToExtract, "IR")
    }
    if ("NMD_status" %in% consequencesToAnalyze) {
      columnsToExtract <- c(columnsToExtract, "PTC")
    }
    if ("coding_potential" %in% consequencesToAnalyze) {
      columnsToExtract <- c(columnsToExtract, "codingPotential")
    }
    if (any(c("sub_cell_location", "sub_cell_shift_to_cell_membrane", 
              "sub_cell_shift_to_cytoplasm", "sub_cell_shift_to_nucleus", 
              "sub_cell_shift_to_Extracellular") %in% consequencesToAnalyze)) {
      columnsToExtract <- c(columnsToExtract, "sub_cell_location")
    }
    transcriptData <- unique(switchAnalyzeRlist$isoformFeatures[which(switchAnalyzeRlist$isoformFeatures$isoform_id %in% 
                                                                        isoformsToAnalyze), columnsToExtract, drop = FALSE])
    columnsToExtract2 <- "isoform_id"
    if ("ORF_genomic" %in% consequencesToAnalyze) {
      columnsToExtract2 <- c(columnsToExtract2, c("orfStartGenomic", 
                                                  "orfEndGenomic"))
    }
    if (any(c("ORF_length", "5_utr_length", "3_utr_length", 
              "ORF_seq_similarity", "5_utr_seq_similarity", "3_utr_seq_similarity", 
              "domains_identified", "genomic_domain_position", 
              "domain_length", "signal_peptide_identified", "IDR_identified", 
              "IDR_length", "IDR_type", "sub_cell_location", "sub_cell_shift_to_cell_membrane", 
              "sub_cell_shift_to_cytoplasm", "sub_cell_shift_to_nucleus", 
              "sub_cell_shift_to_Extracellular", "soform_topology", 
              "intracellular_region_length", "extracellular_region_length", 
              "extracellular_region_count", "intracellular_region_count") %in% 
            consequencesToAnalyze)) {
      columnsToExtract2 <- c(columnsToExtract2, c("orfTransciptStart", 
                                                  "orfTransciptEnd", "orfTransciptLength"))
    }
    if (length(columnsToExtract2) > 1) {
      orfData <- unique(switchAnalyzeRlist$orfAnalysis[which(switchAnalyzeRlist$orfAnalysis$isoform_id %in% 
                                                               isoformsToAnalyze), columnsToExtract2])
      transcriptData <- dplyr::left_join(transcriptData, 
                                         orfData, by = "isoform_id", )
    }
    if ("intron_retention" %in% consequencesToAnalyze) {
      if (is.null(switchAnalyzeRlist$AlternativeSplicingAnalysis)) {
        localIRdata <- switchAnalyzeRlist$intronRetentionAnalysis[which(switchAnalyzeRlist$intronRetentionAnalysis$isoform_id %in% 
                                                                          isoformsToAnalyze), ]
      }
      else {
        localIRdata <- switchAnalyzeRlist$AlternativeSplicingAnalysis[which(switchAnalyzeRlist$AlternativeSplicingAnalysis$isoform_id %in% 
                                                                              isoformsToAnalyze), ]
      }
      if (nrow(localIRdata) != 2) {
        warning(paste("There was a problem with the extraction if intron retentions -", 
                      "please contact the developers with this example so they can fix it.", 
                      "For now they are ignored. The isoforms affected are:", 
                      paste(isoformsToAnalyze, collapse = ", "), 
                      sep = " "))
        consequencesToAnalyze <- consequencesToAnalyze[which(!consequencesToAnalyze %in% 
                                                               c("intron_retention"))]
      }
      else {
        localIRdata$irCoordinats <- paste(localIRdata$IR_genomic_start, 
                                          localIRdata$IR_genomic_end, sep = ":")
      }
    }
    onPlusStrand <- as.character(strand(exonData)[1]) == 
      "+"
    if ("wasTrimmed" %in% colnames(switchAnalyzeRlist$orfAnalysis)) {
      if (onPlusStrand) {
        localTrimmed <- switchAnalyzeRlist$orfAnalysis[which(switchAnalyzeRlist$orfAnalysis$isoform_id %in% 
                                                               isoformsToAnalyze), c("isoform_id", "wasTrimmed", 
                                                                                     "trimmedStartGenomic", "orfEndGenomic")]
      }
      else {
        localTrimmed <- switchAnalyzeRlist$orfAnalysis[which(switchAnalyzeRlist$orfAnalysis$isoform_id %in% 
                                                               isoformsToAnalyze), c("isoform_id", "wasTrimmed", 
                                                                                     "trimmedStartGenomic", "orfStartGenomic")]
        colnames(localTrimmed) <- c("isoform_id", "wasTrimmed", 
                                    "trimmedStartGenomic", "orfEndGenomic")
      }
      if (any(localTrimmed$wasTrimmed, na.rm = TRUE)) {
        localTrimmed <- localTrimmed[which(localTrimmed$wasTrimmed), 
        ]
        regionToOmmit <- GenomicRanges::reduce(IRanges(localTrimmed$trimmedStartGenomic, 
                                                       localTrimmed$orfEndGenomic))
      }
    }
    if (any(c("domains_identified", "genomic_domain_position", 
              "domain_length", "domain_isotype") %in% consequencesToAnalyze)) {
      domanData <- switchAnalyzeRlist$domainAnalysis[which(switchAnalyzeRlist$domainAnalysis$isoform_id %in% 
                                                             isoformsToAnalyze), ]
      domanData$isoform_id <- factor(domanData$isoform_id, 
                                     levels = isoformsToAnalyze)
      colIndex <- na.omit(match(c("hmm_name", "pfamStartGenomic", 
                                  "pfamEndGenomic", "domain_isotype_simple", "orf_aa_start", 
                                  "orf_aa_end"), colnames(domanData)))
      domanDataSplit <- split(domanData[, colIndex], f = domanData$isoform_id)
      if (exists("regionToOmmit")) {
        domanDataSplit <- lapply(domanDataSplit, function(aSet) {
          if (onPlusStrand) {
            aSet[which(!overlapsAny(IRanges(aSet$pfamStartGenomic, 
                                            aSet$pfamEndGenomic), regionToOmmit)), 
            ]
          }
          else {
            aSet[which(!overlapsAny(IRanges(aSet$pfamEndGenomic, 
                                            aSet$pfamStartGenomic), regionToOmmit)), 
            ]
          }
        })
      }
      isNAnames <- transcriptData$isoform_id[which(is.na(transcriptData$orfTransciptLength))]
      if (length(isNAnames)) {
        isNAindex <- which(names(domanDataSplit) %in% 
                             isNAnames)
        domanDataSplit[isNAindex] <- lapply(domanDataSplit[isNAindex], 
                                            function(aDF) {
                                              aDF[0, ]
                                            })
      }
    }
    if (any(c("IDR_identified", "IDR_type", "IDR_length") %in% 
            consequencesToAnalyze)) {
      idrData <- switchAnalyzeRlist$idrAnalysis[which(switchAnalyzeRlist$idrAnalysis$isoform_id %in% 
                                                        isoformsToAnalyze), ]
      idrData$isoform_id <- factor(idrData$isoform_id, 
                                   levels = isoformsToAnalyze)
      idrDataSplit <- split(idrData[, c("idrStartGenomic", 
                                        "idrEndGenomic", "orf_aa_start", "orf_aa_end", 
                                        "idr_type")], f = idrData$isoform_id)
      if (exists("regionToOmmit")) {
        idrDataSplit <- lapply(idrDataSplit, function(aSet) {
          if (onPlusStrand) {
            aSet[which(!overlapsAny(IRanges(aSet$idrStartGenomic, 
                                            aSet$idrEndGenomic), regionToOmmit)), 
            ]
          }
          else {
            aSet[which(!overlapsAny(IRanges(aSet$idrEndGenomic, 
                                            aSet$idrStartGenomic), regionToOmmit)), 
            ]
          }
        })
      }
      isNAnames <- transcriptData$isoform_id[which(is.na(transcriptData$orfTransciptLength))]
      if (length(isNAnames)) {
        isNAindex <- which(names(idrDataSplit) %in% 
                             isNAnames)
        idrDataSplit[isNAindex] <- lapply(idrDataSplit[isNAindex], 
                                          function(aDF) {
                                            aDF[0, ]
                                          })
      }
    }
    if ("signal_peptide_identified" %in% consequencesToAnalyze) {
      peptideData <- switchAnalyzeRlist$signalPeptideAnalysis[which(switchAnalyzeRlist$signalPeptideAnalysis$isoform_id %in% 
                                                                      isoformsToAnalyze), ]
      peptideData$isoform_id <- factor(peptideData$isoform_id, 
                                       levels = isoformsToAnalyze)
      peptideDataSplit <- split(peptideData, f = peptideData$isoform_id)
      if (exists("regionToOmmit")) {
        peptideDataSplit <- lapply(peptideDataSplit, 
                                   function(aSet) {
                                     aSet[which(!overlapsAny(IRanges(aSet$genomicClevageAfter, 
                                                                     aSet$genomicClevageAfter), regionToOmmit)), 
                                     ]
                                   })
      }
      isNAnames <- transcriptData$isoform_id[which(is.na(transcriptData$orfTransciptLength))]
      if (length(isNAnames)) {
        isNAindex <- which(names(peptideDataSplit) %in% 
                             isNAnames)
        peptideDataSplit[isNAindex] <- lapply(peptideDataSplit[isNAindex], 
                                              function(aDF) {
                                                aDF[0, ]
                                              })
      }
    }
    if (any(c("isoform_topology", "intracellular_region_length", 
              "extracellular_region_length", "extracellular_region_count", 
              "intracellular_region_count") %in% consequencesToAnalyze)) {
      topData <- switchAnalyzeRlist$topologyAnalysis[which(switchAnalyzeRlist$topologyAnalysis$isoform_id %in% 
                                                             isoformsToAnalyze), ]
      topData$isoform_id <- factor(topData$isoform_id, 
                                   levels = isoformsToAnalyze)
      colIndex <- na.omit(match(c("region_type", "regionStartGenomic", 
                                  "regionEndGenomic"), colnames(topData)))
      topDataSplit <- split(topData[, colIndex], f = topData$isoform_id)
      if (exists("regionToOmmit")) {
        topDataSplit <- lapply(topDataSplit, function(aSet) {
          if (onPlusStrand) {
            aSet[which(!overlapsAny(IRanges(aSet$regionStartGenomic, 
                                            aSet$regionEndGenomic), regionToOmmit)), 
            ]
          }
          else {
            aSet[which(!overlapsAny(IRanges(aSet$regionEndGenomic, 
                                            aSet$regionStartGenomic), regionToOmmit)), 
            ]
          }
        })
      }
      isNAnames <- transcriptData$isoform_id[which(is.na(transcriptData$orfTransciptLength))]
      if (length(isNAnames)) {
        isNAindex <- which(names(topDataSplit) %in% 
                             isNAnames)
        topDataSplit[isNAindex] <- lapply(topDataSplit[isNAindex], 
                                          function(aDF) {
                                            aDF[0, ]
                                          })
      }
      topDataSplit <- lapply(topDataSplit, function(x) {
        x$regionLength <- abs(x$regionEndGenomic - x$regionStartGenomic)
        return(x)
      })
    }
    if (any(c("isoform_length", "5_utr_length", "3_utr_length", 
              "ORF_length", "isoform_seq_similarity", "5_utr_seq_similarity", 
              "3_utr_seq_similarity", "ORF_seq_similarity") %in% 
            consequencesToAnalyze)) {
      isoform_length <- sapply(exonDataList, function(x) sum(width(x)))
      transcriptData$length <- isoform_length[match(transcriptData$isoform_id, 
                                                    names(isoform_length))]
    }
    if (any(c("isoform_seq_similarity", "5_utr_seq_similarity", 
              "3_utr_seq_similarity") %in% consequencesToAnalyze)) {
      ntSeq <- switchAnalyzeRlist$ntSequence[which(names(switchAnalyzeRlist$ntSequence) %in% 
                                                     c(upIso, downIso))]
      if (length(ntSeq) == 2) {
        upNtSeq <- ntSeq[upIso]
        downNtSeq <- ntSeq[downIso]
      }
      else {
        consequencesToAnalyze <- consequencesToAnalyze[which(!consequencesToAnalyze %in% 
                                                               c("isoform_seq_similarity", "5_utr_seq_similarity", 
                                                                 "3_utr_seq_similarity"))]
      }
    }
    if ("ORF_seq_similarity" %in% consequencesToAnalyze) {
      aaSeq <- switchAnalyzeRlist$aaSequence[which(names(switchAnalyzeRlist$aaSequence) %in% 
                                                     transcriptData$isoform_id[which(!is.na(transcriptData$orfTransciptLength))])]
      if (length(aaSeq) == 2) {
        upAAseq <- aaSeq[upIso]
        downAAseq <- aaSeq[downIso]
      }
    }
    if (length(consequencesToAnalyze) == 0) {
      return(NULL)
    }
    if (nrow(transcriptData) != 2) {
      return(NULL)
    }
  }
  if (TRUE) {
    isoComparison <- data.frame(isoformUpregulated = upIso, 
                                isoformDownregulated = downIso, featureCompared = consequencesToAnalyze, 
                                isoformsDifferent = NA)
    if (addDescription) {
      isoComparison$switchConsequence <- NA
    }
    if ("tss" %in% consequencesToAnalyze) {
      localExonData <- unlist(range(exonDataList))
      if (as.character(localExonData@strand[1]) == "+") {
        tssCoordinats <- start(localExonData)
        tssDifferent <- abs(tssCoordinats[1] - tssCoordinats[2]) > 
          ntCutoff
        mostUpstream <- names(localExonData)[which.min(tssCoordinats)]
      }
      else {
        tssCoordinats <- end(localExonData)
        tssDifferent <- abs(tssCoordinats[1] - tssCoordinats[2]) > 
          ntCutoff
        mostUpstream <- names(localExonData)[which.max(tssCoordinats)]
      }
      localIndex <- which(isoComparison$featureCompared == 
                            "tss")
      isoComparison$isoformsDifferent[localIndex] <- tssDifferent
      if (tssDifferent & addDescription) {
        switchMoreUpstram <- mostUpstream == upIso
        if (switchMoreUpstram) {
          isoComparison$switchConsequence[localIndex] <- "Tss more upstream"
        }
        else {
          isoComparison$switchConsequence[localIndex] <- "Tss more downstream"
        }
      }
    }
    if ("tts" %in% consequencesToAnalyze) {
      localExonData <- unlist(range(exonDataList))
      if (as.character(localExonData@strand[1]) == "+") {
        ttsCoordinats <- end(localExonData)
        ttsDifferent <- abs(ttsCoordinats[1] - ttsCoordinats[2]) > 
          ntCutoff
        mostDownstream <- names(localExonData)[which.max(ttsCoordinats)]
      }
      else {
        ttsCoordinats <- start(localExonData)
        ttsDifferent <- abs(ttsCoordinats[1] - ttsCoordinats[2]) > 
          ntCutoff
        mostDownstream <- names(localExonData)[which.min(ttsCoordinats)]
      }
      localIndex <- which(isoComparison$featureCompared == 
                            "tts")
      isoComparison$isoformsDifferent[localIndex] <- ttsDifferent
      if (ttsDifferent & addDescription) {
        switchMoreDownstream <- mostDownstream == upIso
        if (switchMoreDownstream) {
          isoComparison$switchConsequence[localIndex] <- "Tts more downstream"
        }
        else {
          isoComparison$switchConsequence[localIndex] <- "Tts more upstream"
        }
      }
    }
    if ("last_exon" %in% consequencesToAnalyze) {
      isPlusStrand <- as.logical(exonData@strand[1] == 
                                   "+")
      if (isPlusStrand) {
        lastExons <- lapply(exonDataList, function(x) x[length(x), 
                                                        0])
      }
      else {
        lastExons <- lapply(exonDataList, function(x) x[1, 
                                                        0])
      }
      lastExonDifferent <- !overlapsAny(lastExons[[1]], 
                                        lastExons[[2]])
      localIndex <- which(isoComparison$featureCompared == 
                            "last_exon")
      isoComparison$isoformsDifferent[localIndex] <- lastExonDifferent
      if (lastExonDifferent & addDescription) {
        if (isPlusStrand) {
          localEndCoordinats <- sapply(lastExons, function(aGRange) end(aGRange))
          mostDownstream <- names(localEndCoordinats)[which.max(localEndCoordinats)]
        }
        else {
          localEndCoordinats <- sapply(lastExons, function(aGRange) start(aGRange))
          mostDownstream <- names(localEndCoordinats)[which.min(localEndCoordinats)]
        }
        switchMoreDownstream <- mostDownstream == upIso
        if (switchMoreDownstream) {
          isoComparison$switchConsequence[localIndex] <- "Last exon more downstream"
        }
        else {
          isoComparison$switchConsequence[localIndex] <- "Last exon more upstream"
        }
      }
    }
    if ("isoform_length" %in% consequencesToAnalyze) {
      differentLength <- abs(diff(isoform_length)) > ntCutoff
      if (fractionFilter) {
        downLength <- isoform_length[which(names(isoform_length) == 
                                             downIso)]
        fractionDifference <- abs(diff(isoform_length))/downLength > 
          ntFracCutoff
        differentLength <- differentLength & fractionDifference
      }
      localIndex <- which(isoComparison$featureCompared == 
                            "isoform_length")
      isoComparison$isoformsDifferent[localIndex] <- differentLength
      if (differentLength & addDescription) {
        lengthGain <- names(isoform_length)[which.max(isoform_length)] == 
          upIso
        if (lengthGain) {
          isoComparison$switchConsequence[localIndex] <- "Length gain"
        }
        else {
          isoComparison$switchConsequence[localIndex] <- "Length loss"
        }
      }
    }
    if ("isoform_seq_similarity" %in% consequencesToAnalyze) {
      localAlignment <- Biostrings::pairwiseAlignment(pattern = upNtSeq, 
                                                      subject = downNtSeq, type = "global")
      overlapSize <- min(c(nchar(gsub("-", "", as.character(Biostrings::alignedSubject(localAlignment)))), 
                           nchar(gsub("-", "", as.character(Biostrings::alignedPattern(localAlignment))))))
      totalWidth <- width(localAlignment@subject@unaligned) + 
        width(localAlignment@pattern@unaligned) - overlapSize + 
        1
      jcDist <- overlapSize/totalWidth
      differentIsoformOverlap <- jcDist < ntJCsimCutoff & 
        totalWidth - overlapSize > ntCutoff
      localIndex <- which(isoComparison$featureCompared == 
                            "isoform_seq_similarity")
      isoComparison$isoformsDifferent[localIndex] <- differentIsoformOverlap
      if (differentIsoformOverlap & addDescription) {
        lengthGain <- names(isoform_length)[which.max(isoform_length)] == 
          upIso
        if (lengthGain) {
          isoComparison$switchConsequence[localIndex] <- "Length gain"
        }
        else {
          isoComparison$switchConsequence[localIndex] <- "Length loss"
        }
      }
    }
    if ("exon_number" %in% consequencesToAnalyze) {
      localNrExons <- sapply(exonDataList, length)
      exonDifferent <- localNrExons[1] != localNrExons[2]
      localIndex <- which(isoComparison$featureCompared == 
                            "exon_number")
      isoComparison$isoformsDifferent[localIndex] <- exonDifferent
      if (exonDifferent & addDescription) {
        switchGainsExons <- names(localNrExons)[which.max(localNrExons)] == 
          upIso
        if (switchGainsExons) {
          isoComparison$switchConsequence[localIndex] <- "Exon gain"
        }
        else {
          isoComparison$switchConsequence[localIndex] <- "Exon loss"
        }
      }
    }
    if ("intron_structure" %in% consequencesToAnalyze) {
      localIntrons <- lapply(exonDataList, function(aGR) {
        gaps(ranges(aGR))
      })
      differentintron_structure <- !any(all(localIntrons[[1]] %in% 
                                              localIntrons[[2]]), all(localIntrons[[2]] %in% 
                                                                        localIntrons[[1]]))
      localIndex <- which(isoComparison$featureCompared == 
                            "intron_structure")
      isoComparison$isoformsDifferent[localIndex] <- differentintron_structure
    }
    if ("intron_retention" %in% consequencesToAnalyze) {
      if (all(!is.na(transcriptData$IR))) {
        differentNrIR <- localIRdata$irCoordinats[1] != 
          localIRdata$irCoordinats[2]
        localIndex <- which(isoComparison$featureCompared == 
                              "intron_retention")
        isoComparison$isoformsDifferent[localIndex] <- differentNrIR
        if (differentNrIR & addDescription) {
          if (transcriptData$IR[1] == transcriptData$IR[2]) {
            isoComparison$switchConsequence[localIndex] <- "Intron retention switch"
          }
          else if (transcriptData$isoform_id[which.max(transcriptData$IR)] == 
                   upIso) {
            isoComparison$switchConsequence[localIndex] <- "Intron retention gain"
          }
          else {
            isoComparison$switchConsequence[localIndex] <- "Intron retention loss"
          }
        }
      }
    }
    if ("isoform_class_code" %in% consequencesToAnalyze) {
      differentAnnotation <- transcriptData$class_code[1] != 
        transcriptData$class_code[2]
      localIndex <- which(isoComparison$featureCompared == 
                            "isoform_class_code")
      isoComparison$isoformsDifferent[localIndex] <- differentAnnotation
    }
    if ("ORF_genomic" %in% consequencesToAnalyze) {
      localIndex <- which(isoComparison$featureCompared == 
                            "ORF_genomic")
      if (all(!is.na(transcriptData$orfStartGenomic))) {
        genomicOrfDifferent <- transcriptData$orfStartGenomic[1] != 
          transcriptData$orfStartGenomic[2] | transcriptData$orfEndGenomic[1] != 
          transcriptData$orfEndGenomic[2]
        isoComparison$isoformsDifferent[localIndex] <- genomicOrfDifferent
      }
      else if (sum(!is.na(transcriptData$orfStartGenomic)) == 
               1) {
        isoComparison$isoformsDifferent[localIndex] <- TRUE
      }
      else {
        isoComparison$isoformsDifferent[localIndex] <- FALSE
      }
    }
    if ("ORF_length" %in% consequencesToAnalyze) {
      localIndex <- which(isoComparison$featureCompared == 
                            "ORF_length")
      if (all(!is.na(transcriptData$orfTransciptLength))) {
        orfDifferent <- abs(diff(transcriptData$orfTransciptLength)) > 
          ntCutoff
        if (fractionFilter) {
          downLength <- isoform_length[which(names(isoform_length) == 
                                               downIso)]
          fractionDifference <- abs(diff(transcriptData$orfTransciptLength))/downLength > 
            ntFracCutoff
          orfDifferent <- orfDifferent & fractionDifference
        }
        isoComparison$isoformsDifferent[localIndex] <- orfDifferent
        if (orfDifferent & addDescription) {
          orfGain <- transcriptData$isoform_id[which.max(transcriptData$orfTransciptLength)] == 
            upIso
          if (orfGain) {
            isoComparison$switchConsequence[localIndex] <- "ORF is longer"
          }
          else {
            isoComparison$switchConsequence[localIndex] <- "ORF is shorter"
          }
        }
      }
      else if (sum(!is.na(transcriptData$orfTransciptLength)) == 
               1) {
        isoComparison$isoformsDifferent[localIndex] <- TRUE
        if (addDescription) {
          orfLoss <- is.na(transcriptData$orfTransciptLength[which(transcriptData$isoform_id == 
                                                                     upIso)])
          if (orfLoss) {
            isoComparison$switchConsequence[localIndex] <- "Complete ORF loss"
          }
          else {
            isoComparison$switchConsequence[localIndex] <- "Complete ORF gain"
          }
        }
      }
      else {
        isoComparison$isoformsDifferent[localIndex] <- FALSE
      }
    }
    if ("ORF_seq_similarity" %in% consequencesToAnalyze) {
      localIndex <- which(isoComparison$featureCompared == 
                            "ORF_seq_similarity")
      if (all(!is.na(transcriptData$orfTransciptLength))) {
        if (length(aaSeq) == 2) {
          localAlignment <- Biostrings::pairwiseAlignment(pattern = upAAseq, 
                                                          subject = downAAseq, type = "global")
          overlapSize <- min(c(nchar(gsub("-", "", as.character(Biostrings::alignedSubject(localAlignment)))), 
                               nchar(gsub("-", "", as.character(Biostrings::alignedPattern(localAlignment))))))
          totalWidth <- width(localAlignment@subject@unaligned) + 
            width(localAlignment@pattern@unaligned) - 
            overlapSize + 1
          jcDist <- overlapSize/totalWidth
          differentORFoverlap <- jcDist < AaJCsimCutoff & 
            (totalWidth - overlapSize) > AaCutoff
          isoComparison$isoformsDifferent[localIndex] <- differentORFoverlap
          if (differentORFoverlap & addDescription) {
            lengthGain <- transcriptData$isoform_id[which.max(transcriptData$orfTransciptLength)] == 
              upIso
            if (lengthGain) {
              isoComparison$switchConsequence[localIndex] <- "ORF is longer"
            }
            else {
              isoComparison$switchConsequence[localIndex] <- "ORF is shorter"
            }
          }
        }
      }
      else if (sum(!is.na(transcriptData$orfTransciptLength)) == 
               1) {
        isoComparison$isoformsDifferent[localIndex] <- TRUE
        if (addDescription) {
          orfLoss <- is.na(transcriptData$orfTransciptLength[which(transcriptData$isoform_id == 
                                                                     upIso)])
          if (orfLoss) {
            isoComparison$switchConsequence[localIndex] <- "Complete ORF loss"
          }
          else {
            isoComparison$switchConsequence[localIndex] <- "Complete ORF gain"
          }
        }
      }
      else {
        isoComparison$isoformsDifferent[localIndex] <- FALSE
      }
    }
    if ("5_utr_length" %in% consequencesToAnalyze) {
      if (all(!is.na(transcriptData$orfTransciptStart))) {
        utr5Different <- abs(diff(transcriptData$orfTransciptStart)) > 
          ntCutoff
        if (fractionFilter) {
          downLength <- isoform_length[which(names(isoform_length) == 
                                               downIso)]
          fractionDifference <- abs(diff(transcriptData$orfTransciptStart))/downLength > 
            ntFracCutoff
          utr5Different <- utr5Different & fractionDifference
        }
        localIndex <- which(isoComparison$featureCompared == 
                              "5_utr_length")
        isoComparison$isoformsDifferent[localIndex] <- utr5Different
        if (utr5Different & addDescription) {
          utr5Gain <- transcriptData$isoform_id[which.max(transcriptData$orfTransciptStart)] == 
            upIso
          if (utr5Gain) {
            isoComparison$switchConsequence[localIndex] <- "5UTR is longer"
          }
          else {
            isoComparison$switchConsequence[localIndex] <- "5UTR is shorter"
          }
        }
      }
    }
    if ("5_utr_seq_similarity" %in% consequencesToAnalyze) {
      if (all(!is.na(transcriptData$orfTransciptStart))) {
        localUpNt <- Biostrings::subseq(upNtSeq, 1, 
                                        transcriptData$orfTransciptStart[which(transcriptData$isoform_id == 
                                                                                 upIso)] - 1)
        localDownNt <- Biostrings::subseq(downNtSeq, 
                                          1, transcriptData$orfTransciptStart[which(transcriptData$isoform_id == 
                                                                                      downIso)] - 1)
        if (width(localUpNt) > 0 & width(localDownNt) > 
            0) {
          localAlignment <- Biostrings::pairwiseAlignment(pattern = localUpNt, 
                                                          subject = localDownNt, type = "overlap")
          overlapSize <- min(c(nchar(gsub("-", "", as.character(Biostrings::alignedSubject(localAlignment)))), 
                               nchar(gsub("-", "", as.character(Biostrings::alignedPattern(localAlignment))))))
          totalWidth <- width(localAlignment@subject@unaligned) + 
            width(localAlignment@pattern@unaligned) - 
            overlapSize + 1
          jcDist <- overlapSize/totalWidth
        }
        else {
          overlapSize <- 0
          totalWidth <- abs(width(localUpNt) - width(localDownNt))
          jcDist <- 0
        }
        differenttUTRoverlap <- jcDist < ntJCsimCutoff & 
          totalWidth - overlapSize > ntCutoff
        localIndex <- which(isoComparison$featureCompared == 
                              "5_utr_seq_similarity")
        isoComparison$isoformsDifferent[localIndex] <- differenttUTRoverlap
        if (differenttUTRoverlap & addDescription) {
          utr5Gain <- transcriptData$isoform_id[which.max(transcriptData$orfTransciptStart)] == 
            upIso
          if (utr5Gain) {
            isoComparison$switchConsequence[localIndex] <- "5UTR is longer"
          }
          else {
            isoComparison$switchConsequence[localIndex] <- "5UTR is shorter"
          }
        }
      }
    }
    if ("3_utr_length" %in% consequencesToAnalyze) {
      if (all(!is.na(transcriptData$orfTransciptEnd))) {
        transcriptData$utr3length <- transcriptData$length - 
          (transcriptData$orfTransciptEnd + 1)
        utr3Different <- abs(diff(transcriptData$utr3length)) > 
          ntCutoff
        if (fractionFilter) {
          downLength <- isoform_length[which(names(isoform_length) == 
                                               downIso)]
          fractionDifference <- abs(diff(transcriptData$utr3length))/downLength > 
            ntFracCutoff
          utr3Different <- utr3Different & fractionDifference
        }
        localIndex <- which(isoComparison$featureCompared == 
                              "3_utr_length")
        isoComparison$isoformsDifferent[localIndex] <- utr3Different
        if (utr3Different & addDescription) {
          utr3Gain <- transcriptData$isoform_id[which.max(transcriptData$utr3length)] == 
            upIso
          if (utr3Gain) {
            isoComparison$switchConsequence[localIndex] <- "3UTR is longer"
          }
          else {
            isoComparison$switchConsequence[localIndex] <- "3UTR is shorter"
          }
        }
      }
    }
    if ("3_utr_seq_similarity" %in% consequencesToAnalyze) {
      if (all(!is.na(transcriptData$orfTransciptEnd))) {
        transcriptData$utr3length <- transcriptData$length - 
          (transcriptData$orfTransciptEnd + 1)
        up3UTRstart <- transcriptData$orfTransciptEnd[which(transcriptData$isoform_id == 
                                                              upIso)] + 4
        down3UTRstart <- transcriptData$orfTransciptEnd[which(transcriptData$isoform_id == 
                                                                downIso)] + 4
        if (up3UTRstart > width(upNtSeq) + 1) {
          up3UTRstart <- width(upNtSeq) + 1
        }
        if (down3UTRstart > width(downNtSeq) + 1) {
          down3UTRstart <- width(downNtSeq) + 1
        }
        localUpNt <- subseq(upNtSeq, up3UTRstart, width(upNtSeq))
        localDownNt <- subseq(downNtSeq, down3UTRstart, 
                              width(downNtSeq))
        if (width(localUpNt) > 0 & width(localDownNt) > 
            0) {
          localAlignment <- Biostrings::pairwiseAlignment(pattern = localUpNt, 
                                                          subject = localDownNt, type = "overlap")
          overlapSize <- min(c(nchar(gsub("-", "", as.character(Biostrings::alignedSubject(localAlignment)))), 
                               nchar(gsub("-", "", as.character(Biostrings::alignedPattern(localAlignment))))))
          totalWidth <- width(localAlignment@subject@unaligned) + 
            width(localAlignment@pattern@unaligned) - 
            overlapSize + 1
          jcDist <- overlapSize/totalWidth
        }
        else {
          overlapSize <- 0
          totalWidth <- abs(width(localUpNt) - width(localDownNt))
          jcDist <- 0
        }
        different3UTRoverlap <- jcDist < ntJCsimCutoff & 
          totalWidth - overlapSize > ntCutoff
        localIndex <- which(isoComparison$featureCompared == 
                              "3_utr_seq_similarity")
        isoComparison$isoformsDifferent[localIndex] <- different3UTRoverlap
        if (different3UTRoverlap & addDescription) {
          utr3Gain <- transcriptData$isoform_id[which.max(transcriptData$utr3length)] == 
            upIso
          if (utr3Gain) {
            isoComparison$switchConsequence[localIndex] <- "3UTR is longer"
          }
          else {
            isoComparison$switchConsequence[localIndex] <- "3UTR is shorter"
          }
        }
      }
    }
    if ("NMD_status" %in% consequencesToAnalyze) {
      if (all(!is.na(transcriptData$PTC))) {
        ptcDifferent <- transcriptData$PTC[1] != transcriptData$PTC[2]
        localIndex <- which(isoComparison$featureCompared == 
                              "NMD_status")
        isoComparison$isoformsDifferent[localIndex] <- ptcDifferent
        if (ptcDifferent & addDescription) {
          upIsSensitive <- transcriptData$PTC[which(transcriptData$isoform_id == 
                                                      upIso)]
          if (upIsSensitive) {
            isoComparison$switchConsequence[localIndex] <- "NMD sensitive"
          }
          else {
            isoComparison$switchConsequence[localIndex] <- "NMD insensitive"
          }
        }
      }
    }
    if ("coding_potential" %in% consequencesToAnalyze) {
      if (all(!is.na(transcriptData$codingPotential))) {
        cpDifferent <- transcriptData$codingPotential[1] != 
          transcriptData$codingPotential[2]
        localIndex <- which(isoComparison$featureCompared == 
                              "coding_potential")
        isoComparison$isoformsDifferent[localIndex] <- cpDifferent
        if (cpDifferent & addDescription) {
          upIsCoding <- transcriptData$codingPotential[which(transcriptData$isoform_id == 
                                                               upIso)]
          if (upIsCoding) {
            isoComparison$switchConsequence[localIndex] <- "Transcript is coding"
          }
          else {
            isoComparison$switchConsequence[localIndex] <- "Transcript is Noncoding"
          }
        }
      }
    }
    if ("domains_identified" %in% consequencesToAnalyze) {
      if (sum(!is.na(transcriptData$orfTransciptLength)) > 
          0) {
        domianNames <- lapply(domanDataSplit, function(x) x$hmm_name)
        t1 <- rle(domianNames[[1]])
        t2 <- rle(domianNames[[2]])
        differentDomains <- !identical(t1, t2)
        localIndex <- which(isoComparison$featureCompared == 
                              "domains_identified")
        isoComparison$isoformsDifferent[localIndex] <- differentDomains
        if (differentDomains & addDescription) {
          if (sum(t1$lengths) != sum(t2$lengths)) {
            nrDomains <- sapply(domanDataSplit, nrow)
            upHasMoreDomains <- names(nrDomains)[which.max(nrDomains)] == 
              upIso
            if (upHasMoreDomains) {
              isoComparison$switchConsequence[localIndex] <- "Domain gain"
            }
            else {
              isoComparison$switchConsequence[localIndex] <- "Domain loss"
            }
          }
          else {
            isoComparison$switchConsequence[localIndex] <- "Domain switch"
          }
        }
      }
    }
    if ("genomic_domain_position" %in% consequencesToAnalyze) {
      if (sum(!is.na(transcriptData$orfTransciptLength)) > 
          0) {
        localDomanDataSplit <- lapply(domanDataSplit, 
                                      function(aDF) {
                                        aDF[sort.list(aDF$pfamStartGenomic), c("hmm_name", 
                                                                               "pfamStartGenomic", "pfamEndGenomic")]
                                      })
        genomicDomainDifferent <- !identical(localDomanDataSplit[[1]], 
                                             localDomanDataSplit[[2]])
        localIndex <- which(isoComparison$featureCompared == 
                              "genomic_domain_position")
        isoComparison$isoformsDifferent[localIndex] <- genomicDomainDifferent
      }
    }
    if ("domain_length" %in% consequencesToAnalyze) {
      if (sum(!is.na(transcriptData$orfTransciptLength)) > 
          0) {
        domanDataSplit <- lapply(domanDataSplit, function(x) {
          x$length <- x$orf_aa_end - x$orf_aa_start + 
            1
          return(x)
        })
        nDom <- sapply(domanDataSplit, nrow)
        if (all(nDom)) {
          domainRanges <- lapply(domanDataSplit, function(x) {
            if (x$pfamStartGenomic[1] < x$pfamEndGenomic[1]) {
              GRanges(seqnames = "artificial", IRanges(start = x$pfamStartGenomic, 
                                                       end = x$pfamEndGenomic), hmm_name = x$hmm_name, 
                      orfLength = x$length)
            }
            else {
              GRanges(seqnames = "artificial", IRanges(start = x$pfamEndGenomic, 
                                                       end = x$pfamStartGenomic), hmm_name = x$hmm_name, 
                      orfLength = x$length)
            }
          })
          localOverlap <- findOverlaps(query = domainRanges[[upIso]], 
                                       subject = domainRanges[[downIso]])
          if (length(localOverlap)) {
            localOverlapDf <- as.data.frame(localOverlap)
            localOverlapDf$upName <- domainRanges[[upIso]]$hmm_name[queryHits(localOverlap)]
            localOverlapDf$dnName <- domainRanges[[downIso]]$hmm_name[subjectHits(localOverlap)]
            localOverlapDf$upLength <- domainRanges[[upIso]]$orfLength[queryHits(localOverlap)]
            localOverlapDf$dnLength <- domainRanges[[downIso]]$orfLength[subjectHits(localOverlap)]
            localOverlapDf <- localOverlapDf[which(localOverlapDf$upName == 
                                                     localOverlapDf$dnName), ]
            localOverlapDf$lengthDiff <- abs(localOverlapDf$upLength - 
                                               localOverlapDf$dnLength) > AaCutoff
            if (!is.null(AaFracCutoff)) {
              localOverlapDf$minLength <- apply(localOverlapDf[, 
                                                               c("upLength", "dnLength")], 1, min)
              localOverlapDf$maxLength <- apply(localOverlapDf[, 
                                                               c("upLength", "dnLength")], 1, max)
              localOverlapDf$lengthDiff <- localOverlapDf$minLength/localOverlapDf$maxLength < 
                AaFracCutoff & localOverlapDf$lengthDiff
            }
            localOverlapDfDiff <- localOverlapDf[which(localOverlapDf$lengthDiff), 
            ]
            differentDomainLength <- nrow(localOverlapDfDiff) > 
              0
            localIndex <- which(isoComparison$featureCompared == 
                                  "domain_length")
            isoComparison$isoformsDifferent[localIndex] <- differentDomainLength
          }
          else {
            differentDomainLength <- FALSE
          }
        }
        else {
          differentDomainLength <- FALSE
        }
        if (differentDomainLength & addDescription) {
          localOverlapDfDiff$maxIsUp <- localOverlapDfDiff$maxLength == 
            localOverlapDfDiff$upLength
          if (all(c(TRUE, FALSE) %in% localOverlapDfDiff$maxIsUp)) {
            isoComparison$switchConsequence[localIndex] <- "Mixed Domain length differences"
          }
          else if (all(localOverlapDfDiff$maxIsUp)) {
            isoComparison$switchConsequence[localIndex] <- "Domain length gain"
          }
          else if (all(!localOverlapDfDiff$maxIsUp)) {
            isoComparison$switchConsequence[localIndex] <- "Domain length loss"
          }
          else {
            stop("Something with Domain length analysis went wrong")
          }
        }
      }
    }
    if ("domain_isotype" %in% consequencesToAnalyze) {
      if (sum(!is.na(transcriptData$orfTransciptLength)) > 
          0) {
        domanDataSplit <- lapply(domanDataSplit, function(x) {
          x$length <- x$orf_aa_end - x$orf_aa_start + 
            1
          return(x)
        })
        nDom <- sapply(domanDataSplit, nrow)
        if (all(nDom)) {
          domainRanges <- lapply(domanDataSplit, function(x) {
            if (x$pfamStartGenomic[1] < x$pfamEndGenomic[1]) {
              GRanges(seqnames = "artificial", IRanges(start = x$pfamStartGenomic, 
                                                       end = x$pfamEndGenomic), hmm_name = x$hmm_name, 
                      domain_isotype = x$domain_isotype_simple)
            }
            else {
              GRanges(seqnames = "artificial", IRanges(start = x$pfamEndGenomic, 
                                                       end = x$pfamStartGenomic), hmm_name = x$hmm_name, 
                      domain_isotype = x$domain_isotype_simple)
            }
          })
          localOverlap <- findOverlaps(query = domainRanges[[upIso]], 
                                       subject = domainRanges[[downIso]])
          if (length(localOverlap)) {
            localOverlapDf <- as.data.frame(localOverlap)
            localOverlapDf$upName <- domainRanges[[upIso]]$hmm_name[queryHits(localOverlap)]
            localOverlapDf$dnName <- domainRanges[[downIso]]$hmm_name[subjectHits(localOverlap)]
            localOverlapDf$upStructure <- domainRanges[[upIso]]$domain_isotype[queryHits(localOverlap)]
            localOverlapDf$dnStructure <- domainRanges[[downIso]]$domain_isotype[subjectHits(localOverlap)]
            localOverlapDf <- localOverlapDf[which(localOverlapDf$upName == 
                                                     localOverlapDf$dnName), ]
            if (nrow(localOverlapDf)) {
              localOverlapDf$upCombined <- stringr::str_c(localOverlapDf$upName, 
                                                          "_", localOverlapDf$upStructure)
              localOverlapDf$dnCombined <- stringr::str_c(localOverlapDf$dnName, 
                                                          "_", localOverlapDf$dnStructure)
              differentDomainStructure <- any(localOverlapDf$upCombined != 
                                                localOverlapDf$dnCombined)
            }
            else {
              differentDomainStructure <- FALSE
            }
          }
          else {
            differentDomainStructure <- FALSE
          }
        }
        else {
          differentDomainStructure <- FALSE
        }
        localIndex <- which(isoComparison$featureCompared == 
                              "domain_isotype")
        isoComparison$isoformsDifferent[localIndex] <- differentDomainStructure
        if (differentDomainStructure & addDescription) {
          upCount <- sum(localOverlapDf$upStructure != 
                           "Reference")
          dnCount <- sum(localOverlapDf$dnStructure != 
                           "Reference")
          if (all(c(upCount, dnCount) > 0) & upCount == 
              dnCount) {
            isoComparison$switchConsequence[localIndex] <- "Mixed domain isotype changes"
          }
          else if (upCount > dnCount) {
            isoComparison$switchConsequence[localIndex] <- "Domain non-reference isotype gain"
          }
          else if (dnCount > upCount) {
            isoComparison$switchConsequence[localIndex] <- "Domain non-reference isotype loss"
          }
        }
      }
    }
    if ("IDR_identified" %in% consequencesToAnalyze) {
      if (sum(!is.na(transcriptData$orfTransciptLength)) > 
          0) {
        nIdr <- sapply(idrDataSplit, nrow)
        if (all(nIdr)) {
          idrRanges <- lapply(idrDataSplit, function(x) {
            if (x$idrStartGenomic[1] < x$idrEndGenomic[1]) {
              IRanges(start = x$idrStartGenomic, end = x$idrEndGenomic)
            }
            else {
              IRanges(start = x$idrEndGenomic, end = x$idrStartGenomic)
            }
          })
          overlap1 <- overlapsAny(idrRanges[[1]], idrRanges[[2]])
          overlap2 <- overlapsAny(idrRanges[[2]], idrRanges[[1]])
        }
        else if (nIdr[1] == 0 & nIdr[2] == 0) {
          overlap1 <- NA
          overlap2 <- NA
        }
        else if (nIdr[1] == 0) {
          overlap1 <- TRUE
          overlap2 <- FALSE
        }
        else {
          overlap1 <- FALSE
          overlap2 <- TRUE
        }
        differentIdr <- any(c(!overlap1, !overlap2), 
                            na.rm = TRUE)
        localIndex <- which(isoComparison$featureCompared == 
                              "IDR_identified")
        isoComparison$isoformsDifferent[localIndex] <- differentIdr
        if (differentIdr & addDescription) {
          if (any(!overlap1) & any(!overlap2)) {
            isoComparison$switchConsequence[localIndex] <- "IDR switch"
          }
          else {
            if (any(!overlap1)) {
              theDiff <- 1
            }
            else {
              theDiff <- 2
            }
            upHasMoreIDR <- names(idrDataSplit)[theDiff] == 
              upIso
            if (upHasMoreIDR) {
              isoComparison$switchConsequence[localIndex] <- "IDR gain"
            }
            else {
              isoComparison$switchConsequence[localIndex] <- "IDR loss"
            }
          }
        }
      }
    }
    if ("IDR_type" %in% consequencesToAnalyze) {
      if (sum(!is.na(transcriptData$orfTransciptLength)) > 
          0) {
        nIdr <- sapply(idrDataSplit, nrow)
        if (all(nIdr)) {
          idrRanges <- lapply(idrDataSplit, function(x) {
            if (x$idrStartGenomic[1] < x$idrEndGenomic[1]) {
              GRanges(seqnames = "artificial", IRanges(start = x$idrStartGenomic, 
                                                       end = x$idrEndGenomic), type = x$idr_type)
            }
            else {
              GRanges(seqnames = "artificial", IRanges(start = x$idrEndGenomic, 
                                                       end = x$idrStartGenomic), type = x$idr_type)
            }
          })
          localOverlap <- findOverlaps(query = idrRanges[[upIso]], 
                                       subject = idrRanges[[downIso]])
          if (length(localOverlap)) {
            localOverlapDf <- as.data.frame(localOverlap)
            localOverlapDf$upType <- idrRanges[[upIso]]$type[queryHits(localOverlap)]
            localOverlapDf$dnType <- idrRanges[[downIso]]$type[subjectHits(localOverlap)]
            differentIdrType <- any(na.omit(localOverlapDf$upType != 
                                              localOverlapDf$dnType))
            localIndex <- which(isoComparison$featureCompared == 
                                  "IDR_type")
            isoComparison$isoformsDifferent[localIndex] <- differentIdrType
          }
          else {
            differentIdrType <- FALSE
          }
        }
        else {
          differentIdrType <- FALSE
        }
        if (differentIdrType & addDescription) {
          localOverlapDf <- localOverlapDf[which(localOverlapDf$upType != 
                                                   localOverlapDf$dnType), ]
          localOverlapDf$bindingGain <- localOverlapDf$dnType == 
            "IDR" & localOverlapDf$upType == "IDR_w_binding_region"
          localOverlapDf$bindingLoss <- localOverlapDf$dnType == 
            "IDR_w_binding_region" & localOverlapDf$upType == 
            "IDR"
          if (any(localOverlapDf$bindingGain) & any(localOverlapDf$bindingLoss)) {
            isoComparison$switchConsequence[localIndex] <- "IDR w binding region switch"
          }
          else if (any(localOverlapDf$bindingGain) & 
                   !any(localOverlapDf$bindingLoss)) {
            isoComparison$switchConsequence[localIndex] <- "IDR w binding region gain"
          }
          else if (!any(localOverlapDf$bindingGain) & 
                   any(localOverlapDf$bindingLoss)) {
            isoComparison$switchConsequence[localIndex] <- "IDR w binding region loss"
          }
          else {
            stop("Something with idr binding regions went wrong")
          }
        }
      }
    }
    if ("IDR_length" %in% consequencesToAnalyze) {
      if (sum(!is.na(transcriptData$orfTransciptLength)) > 
          0) {
        idrDataSplit <- lapply(idrDataSplit, function(x) {
          x$length <- x$orf_aa_end - x$orf_aa_start + 
            1
          return(x)
        })
        nIdr <- sapply(idrDataSplit, nrow)
        if (all(nIdr)) {
          idrRanges <- lapply(idrDataSplit, function(x) {
            if (x$idrStartGenomic[1] < x$idrEndGenomic[1]) {
              GRanges(seqnames = "artificial", IRanges(start = x$idrStartGenomic, 
                                                       end = x$idrEndGenomic), type = x$idr_type, 
                      orfLength = x$length)
            }
            else {
              GRanges(seqnames = "artificial", IRanges(start = x$idrEndGenomic, 
                                                       end = x$idrStartGenomic), type = x$idr_type, 
                      orfLength = x$length)
            }
          })
          localOverlap <- findOverlaps(query = idrRanges[[upIso]], 
                                       subject = idrRanges[[downIso]])
          if (length(localOverlap)) {
            localOverlapDf <- as.data.frame(localOverlap)
            localOverlapDf$upLength <- idrRanges[[upIso]]$orfLength[queryHits(localOverlap)]
            localOverlapDf$dnLength <- idrRanges[[downIso]]$orfLength[subjectHits(localOverlap)]
            localOverlapDf$lengthDiff <- abs(localOverlapDf$upLength - 
                                               localOverlapDf$dnLength) > AaCutoff
            if (!is.null(AaFracCutoff)) {
              localOverlapDf$minLength <- apply(localOverlapDf[, 
                                                               c("upLength", "dnLength")], 1, min)
              localOverlapDf$maxLength <- apply(localOverlapDf[, 
                                                               c("upLength", "dnLength")], 1, max)
              localOverlapDf$lengthDiff <- localOverlapDf$minLength/localOverlapDf$maxLength < 
                AaFracCutoff & localOverlapDf$lengthDiff
            }
            localOverlapDfDiff <- localOverlapDf[which(localOverlapDf$lengthDiff), 
            ]
            differentIdrLength <- nrow(localOverlapDfDiff) > 
              0
            localIndex <- which(isoComparison$featureCompared == 
                                  "IDR_length")
            isoComparison$isoformsDifferent[localIndex] <- differentIdrLength
          }
          else {
            differentIdrType <- FALSE
          }
        }
        else {
          differentIdrType <- FALSE
        }
        if (differentIdrType & addDescription) {
          localOverlapDfDiff$maxIsUp <- localOverlapDfDiff$maxLength == 
            localOverlapDfDiff$upLength
          if (all(c(TRUE, FALSE) %in% localOverlapDfDiff$maxIsUp)) {
            isoComparison$switchConsequence[localIndex] <- "Mixed IDR length differences"
          }
          else if (all(localOverlapDfDiff$maxIsUp)) {
            isoComparison$switchConsequence[localIndex] <- "IDR length gain"
          }
          else if (all(!localOverlapDfDiff$maxIsUp)) {
            isoComparison$switchConsequence[localIndex] <- "IDR length loss"
          }
          else {
            stop("Something with idr length went wrong")
          }
        }
      }
    }
    if ("signal_peptide_identified" %in% consequencesToAnalyze) {
      if (sum(!is.na(transcriptData$orfTransciptLength)) > 
          0) {
        nrSignalPeptides <- sapply(peptideDataSplit, 
                                   nrow)
        differentSignaPeptied <- nrSignalPeptides[1] != 
          nrSignalPeptides[2]
        localIndex <- which(isoComparison$featureCompared == 
                              "signal_peptide_identified")
        isoComparison$isoformsDifferent[localIndex] <- differentSignaPeptied
        if (differentSignaPeptied & addDescription) {
          upHasSignal <- names(nrSignalPeptides)[which.max(nrSignalPeptides)] == 
            upIso
          if (upHasSignal) {
            isoComparison$switchConsequence[localIndex] <- "Signal peptide gain"
          }
          else {
            isoComparison$switchConsequence[localIndex] <- "Signal peptide loss"
          }
        }
      }
    }
    if ("sub_cell_location" %in% consequencesToAnalyze) {
      if (sum(!is.na(transcriptData$orfTransciptLength)) == 
          2) {
        if (sum(!is.na(transcriptData$sub_cell_location)) == 
            2) {
          upIsoIndex <- which(transcriptData$isoform_id == 
                                upIso)
          dnIsoIndex <- which(transcriptData$isoform_id != 
                                upIso)
          locUp <- sort(unlist(strsplit(transcriptData$sub_cell_location[upIsoIndex], 
                                        ",")))
          locDn <- sort(unlist(strsplit(transcriptData$sub_cell_location[dnIsoIndex], 
                                        ",")))
          differentLoc <- !identical(locUp, locDn)
          localIndex <- which(isoComparison$featureCompared == 
                                "sub_cell_location")
          isoComparison$isoformsDifferent[localIndex] <- differentLoc
          if (differentLoc & addDescription) {
            upUniqueLength <- length(setdiff(locUp, 
                                             locDn))
            dnUniqueLength <- length(setdiff(locDn, 
                                             locUp))
            if (upUniqueLength > 0 & dnUniqueLength > 
                0) {
              isoComparison$switchConsequence[localIndex] <- "SubCell location switch"
            }
            if (upUniqueLength > 0 & dnUniqueLength == 
                0) {
              isoComparison$switchConsequence[localIndex] <- "SubCell location gain"
            }
            if (upUniqueLength == 0 & dnUniqueLength > 
                0) {
              isoComparison$switchConsequence[localIndex] <- "SubCell location loss"
            }
          }
        }
      }
    }
    if ("sub_cell_shift_to_cell_membrane" %in% consequencesToAnalyze) {
      if (sum(!is.na(transcriptData$orfTransciptLength)) == 
          2) {
        if (sum(!is.na(transcriptData$sub_cell_location)) == 
            2) {
          upIsoIndex <- which(transcriptData$isoform_id == 
                                upIso)
          dnIsoIndex <- which(transcriptData$isoform_id != 
                                upIso)
          upLength <- length(intersect(unlist(strsplit(transcriptData$sub_cell_location[upIsoIndex], 
                                                       ",")), "Cell_membrane"))
          dnLength <- length(intersect(unlist(strsplit(transcriptData$sub_cell_location[dnIsoIndex], 
                                                       ",")), "Cell_membrane"))
          differentLoc <- upLength != dnLength
          localIndex <- which(isoComparison$featureCompared == 
                                "sub_cell_shift_to_cell_membrane")
          isoComparison$isoformsDifferent[localIndex] <- differentLoc
          if (differentLoc & addDescription) {
            if (upLength > 0 & dnLength == 0) {
              isoComparison$switchConsequence[localIndex] <- "SubCell location memb gain"
            }
            if (upLength == 0 & dnLength > 0) {
              isoComparison$switchConsequence[localIndex] <- "SubCell location memb loss"
            }
          }
        }
      }
    }
    if ("sub_cell_shift_to_cytoplasm" %in% consequencesToAnalyze) {
      if (sum(!is.na(transcriptData$orfTransciptLength)) == 
          2) {
        if (sum(!is.na(transcriptData$sub_cell_location)) == 
            2) {
          upIsoIndex <- which(transcriptData$isoform_id == 
                                upIso)
          dnIsoIndex <- which(transcriptData$isoform_id != 
                                upIso)
          upLength <- length(intersect(unlist(strsplit(transcriptData$sub_cell_location[upIsoIndex], 
                                                       ",")), "Cytoplasm"))
          dnLength <- length(intersect(unlist(strsplit(transcriptData$sub_cell_location[dnIsoIndex], 
                                                       ",")), "Cytoplasm"))
          differentLoc <- upLength != dnLength
          localIndex <- which(isoComparison$featureCompared == 
                                "sub_cell_shift_to_cytoplasm")
          isoComparison$isoformsDifferent[localIndex] <- differentLoc
          if (differentLoc & addDescription) {
            if (upLength > 0 & dnLength == 0) {
              isoComparison$switchConsequence[localIndex] <- "SubCell location cyto gain"
            }
            if (upLength == 0 & dnLength > 0) {
              isoComparison$switchConsequence[localIndex] <- "SubCell location cyto loss"
            }
          }
        }
      }
    }
    if ("sub_cell_shift_to_nucleus" %in% consequencesToAnalyze) {
      if (sum(!is.na(transcriptData$orfTransciptLength)) == 
          2) {
        if (sum(!is.na(transcriptData$sub_cell_location)) == 
            2) {
          upIsoIndex <- which(transcriptData$isoform_id == 
                                upIso)
          dnIsoIndex <- which(transcriptData$isoform_id != 
                                upIso)
          upLength <- length(intersect(unlist(strsplit(transcriptData$sub_cell_location[upIsoIndex], 
                                                       ",")), "Nucleus"))
          dnLength <- length(intersect(unlist(strsplit(transcriptData$sub_cell_location[dnIsoIndex], 
                                                       ",")), "Nucleus"))
          differentLoc <- upLength != dnLength
          localIndex <- which(isoComparison$featureCompared == 
                                "sub_cell_shift_to_nucleus")
          isoComparison$isoformsDifferent[localIndex] <- differentLoc
          if (differentLoc & addDescription) {
            if (upLength > 0 & dnLength == 0) {
              isoComparison$switchConsequence[localIndex] <- "SubCell location nucl gain"
            }
            if (upLength == 0 & dnLength > 0) {
              isoComparison$switchConsequence[localIndex] <- "SubCell location nucl loss"
            }
          }
        }
      }
    }
    if ("sub_cell_shift_to_Extracellular" %in% consequencesToAnalyze) {
      if (sum(!is.na(transcriptData$orfTransciptLength)) == 
          2) {
        if (sum(!is.na(transcriptData$sub_cell_location)) == 
            2) {
          upIsoIndex <- which(transcriptData$isoform_id == 
                                upIso)
          dnIsoIndex <- which(transcriptData$isoform_id != 
                                upIso)
          upLength <- length(intersect(unlist(strsplit(transcriptData$sub_cell_location[upIsoIndex], 
                                                       ",")), "Extracellular"))
          dnLength <- length(intersect(unlist(strsplit(transcriptData$sub_cell_location[dnIsoIndex], 
                                                       ",")), "Extracellular"))
          differentLoc <- upLength != dnLength
          localIndex <- which(isoComparison$featureCompared == 
                                "sub_cell_shift_to_Extracellular")
          isoComparison$isoformsDifferent[localIndex] <- differentLoc
          if (differentLoc & addDescription) {
            if (upLength > 0 & dnLength == 0) {
              isoComparison$switchConsequence[localIndex] <- "SubCell location ext cell gain"
            }
            if (upLength == 0 & dnLength > 0) {
              isoComparison$switchConsequence[localIndex] <- "SubCell location ext cell loss"
            }
          }
        }
      }
    }
    if ("isoform_topology" %in% consequencesToAnalyze) {
      if (sum(!is.na(transcriptData$orfTransciptLength)) > 
          0) {
        toplogy <- lapply(topDataSplit, function(x) x$region_type)
        t1 <- rle(toplogy[[1]])
        t2 <- rle(toplogy[[2]])
        differentTopology <- !identical(t1, t2)
        localIndex <- which(isoComparison$featureCompared == 
                              "isoform_topology")
        isoComparison$isoformsDifferent[localIndex] <- differentTopology
        if (differentTopology & addDescription) {
          if (sum(t1$lengths) != sum(t2$lengths)) {
            nrTop <- sapply(topDataSplit, nrow)
            upHasMoreTop <- names(nrTop)[which.max(nrTop)] == 
              upIso
            if (upHasMoreTop) {
              isoComparison$switchConsequence[localIndex] <- "Topology complexity gain"
            }
            else {
              isoComparison$switchConsequence[localIndex] <- "Topology complexity loss"
            }
          }
          else {
            isoComparison$switchConsequence[localIndex] <- NA
          }
        }
      }
    }
    if ("extracellular_region_count" %in% consequencesToAnalyze) {
      if (sum(!is.na(transcriptData$orfTransciptLength)) > 
          0) {
        upCount <- sum(topDataSplit[[upIso]]$region_type == 
                         "outside")
        dnCount <- sum(topDataSplit[[downIso]]$region_type == 
                         "outside")
        differentTopology <- upCount != dnCount
        localIndex <- which(isoComparison$featureCompared == 
                              "extracellular_region_count")
        isoComparison$isoformsDifferent[localIndex] <- differentTopology
        if (differentTopology & addDescription) {
          if (upCount > dnCount) {
            isoComparison$switchConsequence[localIndex] <- "Extracellular region gain"
          }
          else {
            isoComparison$switchConsequence[localIndex] <- "Extracellular region loss"
          }
        }
      }
    }
    if ("intracellular_region_count" %in% consequencesToAnalyze) {
      if (sum(!is.na(transcriptData$orfTransciptLength)) > 
          0) {
        upCount <- sum(topDataSplit[[upIso]]$region_type == 
                         "inside")
        dnCount <- sum(topDataSplit[[downIso]]$region_type == 
                         "inside")
        differentTopology <- upCount != dnCount
        localIndex <- which(isoComparison$featureCompared == 
                              "intracellular_region_count")
        isoComparison$isoformsDifferent[localIndex] <- differentTopology
        if (differentTopology & addDescription) {
          if (upCount > dnCount) {
            isoComparison$switchConsequence[localIndex] <- "Intracellular region gain"
          }
          else {
            isoComparison$switchConsequence[localIndex] <- "Intracellular region loss"
          }
        }
      }
    }
    if ("extracellular_region_length" %in% consequencesToAnalyze) {
      if (sum(!is.na(transcriptData$orfTransciptLength)) > 
          0) {
        upLength <- sum(topDataSplit[[upIso]]$regionLength[which(topDataSplit[[upIso]]$region_type == 
                                                                   "outside")])
        dnLength <- sum(topDataSplit[[downIso]]$regionLength[which(topDataSplit[[downIso]]$region_type == 
                                                                     "outside")])
        lengthGain <- upLength - dnLength
        minLength <- min(c(upLength, dnLength))
        maxLength <- max(c(upLength, dnLength))
        differentLength <- abs(lengthGain) > AaCutoff & 
          (minLength/maxLength) < AaFracCutoff
        localIndex <- which(isoComparison$featureCompared == 
                              "extracellular_region_length")
        isoComparison$isoformsDifferent[localIndex] <- differentLength
        if (differentLength & addDescription) {
          if (lengthGain > 0) {
            isoComparison$switchConsequence[localIndex] <- "Extracellular length gain"
          }
          else {
            isoComparison$switchConsequence[localIndex] <- "Extracellular length loss"
          }
        }
      }
    }
    if ("intracellular_region_length" %in% consequencesToAnalyze) {
      if (sum(!is.na(transcriptData$orfTransciptLength)) > 
          0) {
        upLength <- sum(topDataSplit[[upIso]]$regionLength[which(topDataSplit[[upIso]]$region_type == 
                                                                   "inside")])
        dnLength <- sum(topDataSplit[[downIso]]$regionLength[which(topDataSplit[[downIso]]$region_type == 
                                                                     "inside")])
        lengthGain <- upLength - dnLength
        minLength <- min(c(upLength, dnLength))
        maxLength <- max(c(upLength, dnLength))
        differentLength <- abs(lengthGain) > AaCutoff & 
          (minLength/maxLength) < AaFracCutoff
        localIndex <- which(isoComparison$featureCompared == 
                              "intracellular_region_length")
        isoComparison$isoformsDifferent[localIndex] <- differentLength
        if (differentLength & addDescription) {
          if (lengthGain > 0) {
            isoComparison$switchConsequence[localIndex] <- "Intracellular length gain"
          }
          else {
            isoComparison$switchConsequence[localIndex] <- "Intracellular length loss"
          }
        }
      }
    }
  }
  if (onlyRepportDifferent) {
    isoComparison <- isoComparison[which(isoComparison$isoformsDifferent), 
    ]
  }
  return(isoComparison)
}

extractSwitchPairs <- function (switchAnalyzeRlist, alpha = 0.05, dIFcutoff = 0.1, 
          onlySigIsoforms = FALSE) 
{
  if (TRUE) {
    if (class(switchAnalyzeRlist) != "switchAnalyzeRlist") {
      stop("The object supplied to 'switchAnalyzeRlist' is not a 'switchAnalyzeRlist'")
    }
    if (alpha < 0 | alpha > 1) {
      warning("The alpha parameter should usually be between 0 and 1 ([0,1]).")
    }
    if (alpha > 0.05) {
      warning("Most journals and scientists consider an alpha larger than 0.05 untrustworthy. We therefore recommend using alpha values smaller than or queal to 0.05")
    }
    if (!any(!is.na(switchAnalyzeRlist$isoformFeatures$gene_switch_q_value))) {
      stop("The analsis of isoform switching must be performed before functional consequences can be analyzed. Please run ?detectIsoformSwitching and try again.")
    }
  }
  if (TRUE) {
    localData <- switchAnalyzeRlist$isoformFeatures[which(switchAnalyzeRlist$isoformFeatures$gene_switch_q_value < 
                                                            alpha & abs(switchAnalyzeRlist$isoformFeatures$dIF) > 
                                                            dIFcutoff), c("iso_ref", "gene_ref", "isoform_switch_q_value", 
                                                                          "gene_switch_q_value", "dIF")]
    if (!nrow(localData)) {
      stop("No genes were considered switching with the used cutoff values")
    }
    localData$switchDirection <- NA
    localData$switchDirection[which(sign(localData$dIF) == 
                                      1)] <- "up"
    localData$switchDirection[which(sign(localData$dIF) == 
                                      -1)] <- "down"
    isoResTest <- any(!is.na(switchAnalyzeRlist$isoformFeatures$isoform_switch_q_value))
    if (isoResTest) {
      localData$isoSig <- localData$isoform_switch_q_value < 
        alpha & abs(localData$dIF) > dIFcutoff
    }
    else {
      localData$isoSig <- localData$gene_switch_q_value < 
        alpha & abs(localData$dIF) > dIFcutoff
    }
    if (onlySigIsoforms) {
      localData <- localData[which(localData$isoSig), 
      ]
    }
  }
  if (TRUE) {
    sigUpData <- localData[which(localData$isoSig & localData$switchDirection == 
                                   "up"), c("iso_ref", "gene_ref")]
    sigDnData <- localData[which(localData$isoSig & localData$switchDirection == 
                                   "down"), c("iso_ref", "gene_ref")]
    colnames(sigUpData)[1] <- c("iso_ref_up")
    colnames(sigDnData)[1] <- c("iso_ref_down")
    if (!onlySigIsoforms) {
      justUpData <- localData[which(localData$switchDirection == 
                                      "up"), c("iso_ref", "gene_ref")]
      justDnData <- localData[which(localData$switchDirection == 
                                      "down"), c("iso_ref", "gene_ref")]
      colnames(justUpData)[1] <- c("iso_ref_up")
      colnames(justDnData)[1] <- c("iso_ref_down")
    }
  }
  if (TRUE) {
    if (onlySigIsoforms) {
      pairwiseIsoComparison <- dplyr::inner_join(sigUpData, 
                                                 sigDnData, by = "gene_ref", multiple = "all")
    }
    else {
      upPairs <- dplyr::inner_join(sigUpData, justDnData, 
                                   by = "gene_ref", multiple = "all")
      dnPairs <- dplyr::inner_join(justUpData, sigDnData, 
                                   by = "gene_ref", multiple = "all")
      pairwiseIsoComparison <- unique(rbind(upPairs, dnPairs))
      pairwiseIsoComparison <- pairwiseIsoComparison[, 
                                                     c("gene_ref", "iso_ref_up", "iso_ref_down")]
      pairwiseIsoComparison <- pairwiseIsoComparison[order(pairwiseIsoComparison$gene_ref, 
                                                           pairwiseIsoComparison$iso_ref_up, pairwiseIsoComparison$iso_ref_down), 
      ]
    }
  }
  if (TRUE) {
    matchVectorUp <- match(pairwiseIsoComparison$iso_ref_up, 
                           switchAnalyzeRlist$isoformFeatures$iso_ref)
    matchVectorDn <- match(pairwiseIsoComparison$iso_ref_down, 
                           switchAnalyzeRlist$isoformFeatures$iso_ref)
    pairwiseIsoComparison$isoformUpregulated <- switchAnalyzeRlist$isoformFeatures$isoform_id[matchVectorUp]
    pairwiseIsoComparison$isoformDownregulated <- switchAnalyzeRlist$isoformFeatures$isoform_id[matchVectorDn]
    pairwiseIsoComparison$gene_id <- switchAnalyzeRlist$isoformFeatures$gene_id[matchVectorUp]
    pairwiseIsoComparison$gene_name <- switchAnalyzeRlist$isoformFeatures$gene_name[matchVectorUp]
    pairwiseIsoComparison$condition_1 <- switchAnalyzeRlist$isoformFeatures$condition_1[matchVectorUp]
    pairwiseIsoComparison$condition_2 <- switchAnalyzeRlist$isoformFeatures$condition_2[matchVectorUp]
  }
  return(pairwiseIsoComparison)
}

myListToDf <- function (aList, ignoreColNames = FALSE, addOrignAsRowNames = FALSE, 
          addOrignAsColumn = FALSE, addOrgRownames = FALSE) 
{
  if (class(aList) != "list") {
    stop("Input is not a list")
  }
  aList <- aList[which(!sapply(aList, is.null))]
  if (class(aList[[1]]) != "data.frame") {
    aList <- lapply(aList, function(x) as.data.frame(t(x)))
  }
  nCol <- unique(sapply(aList, ncol))
  if (length(nCol) != 1) {
    stop("Interies in the list does not have the same number of collums/")
  }
  if (!ignoreColNames) {
    if (length(unique(as.vector(sapply(aList, names)))) != 
        nCol) {
      stop("Interies in the list does not have the collum names")
    }
  }
  df <- data.frame(matrix(NA, ncol = nCol, nrow = sum(sapply(aList, 
                                                             nrow))))
  for (i in 1:nCol) {
    df[, i] <- as.vector(unlist(sapply(aList, function(x) x[, 
                                                            i])))
  }
  colnames(df) <- colnames(aList[[1]])
  if (addOrignAsColumn) {
    df$orign <- rep(names(aList), sapply(aList, nrow))
  }
  if (addOrignAsRowNames) {
    rownames(df) <- rep(names(aList), sapply(aList, nrow))
  }
  if (addOrgRownames) {
    rownames(df) <- rep(sapply(aList, rownames), sapply(aList, 
                                                        nrow))
  }
  return(df)
}


makeMinimumSwitchList <- function (orgSwitchList, isoformsToKeep) 
{
  if (class(orgSwitchList) != "switchAnalyzeRlist") {
    stop("The object supplied to 'orgSwitchList' must be a 'switchAnalyzeRlist'")
  }
  orgSwitchList <- subsetSwitchAnalyzeRlist(orgSwitchList, 
                                            orgSwitchList$isoformFeatures$isoform_id %in% isoformsToKeep)
  orgSwitchList$isoformSwitchAnalysis <- NULL
  orgSwitchList$switchConsequence <- NULL
  colsToAnnulate <- c("gene_name", "gene_value_1", "gene_value_2", 
                      "gene_stderr_1", "gene_stderr_2", "gene_log2_fold_change", 
                      "gene_q_value", "iso_value_1", "iso_value_2", "iso_stderr_1", 
                      "iso_stderr_2", "iso_log2_fold_change", "iso_q_value", 
                      "IF1", "IF2", "dIF", "isoform_switch_q_value")
  orgSwitchList$isoformFeatures <- orgSwitchList$isoformFeatures[, 
                                                                 which(!colnames(orgSwitchList$isoformFeatures) %in% 
                                                                         colsToAnnulate)]
  orgSwitchList$isoformFeatures$condition_1 <- 1
  orgSwitchList$isoformFeatures$condition_2 <- 2
  orgSwitchList$isoformFeatures$gene_switch_q_value <- 1
  orgSwitchList$isoformFeatures$isoform_switch_q_value <- 1
  orgSwitchList$isoformFeatures <- unique(orgSwitchList$isoformFeatures)
  return(orgSwitchList)
}
