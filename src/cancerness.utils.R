library(dplyr)
library(stringr)
library(pROC)
library(ggembl)
library(tidyverse)
library(ggforce)
library(ggbeeswarm)



quantiles <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1)
evaluated.color <- '#21520e'  

groupColors <- c('#9D0208', evaluated.color, '#0353A4')

colorVec <- c(colorRampPalette(c(groupColors[1], 'white'))(length(quantiles)),
              colorRampPalette(c(groupColors[3], 'white'))(length(quantiles)),
              colorRampPalette(c(groupColors[2], 'white'))(length(quantiles)))

colorVecEval <- colorRampPalette(c(evaluated.color, 'white'))(length(quantiles))
colorVecRed <- colorRampPalette(c('#9D0208', 'white'))(length(quantiles))
colorVecBlue <- colorRampPalette(c('#0353A4', 'white'))(length(quantiles))

pc <- 1E-5



# This is a functional quantile plot
get_quantile_plot <- function(inputData, axisColumn, labelColumn, valueColumn, expectedNumLevels = 3, xlab = 'Genera', ylab = 'Relative Abundances (log10)') {
  
  colnames(inputData)[colnames(inputData) == axisColumn] <- 'genus'
  colnames(inputData)[colnames(inputData) == labelColumn] <- 'label'
  colnames(inputData)[colnames(inputData) == valueColumn] <- 'relAb'
  
  inputData$plotGroup <- 1
  print(head(inputData))
  quantileData <- list()
  for (quantile in quantiles) {
    tmp <- inputData %>% 
      #select(cyl, wt) %>% 
      group_by(genus, label, plotGroup) %>%
      #mutate(cyl) %>%
      summarize(value_max = quantile(relAb, probs = quantile),
                value_min = quantile(relAb, probs = 1-quantile)) %>%
      mutate(quantile = quantile)
    quantileData[[length(quantileData) + 1]] <- tmp
  }
  
  quantileData <- do.call('rbind', quantileData) %>%

    mutate(genus = as.factor(genus)) %>%
    arrange(desc(quantile))  %>%
    mutate(Quantiles = map2_chr(quantile, label, function(x, y) {
      return(str_c((x * 100), '% - ', y))
    }))
  
  l <- quantileData %>%
    ungroup() %>%
    select(Quantiles, quantile, label) %>%
    distinct() %>%
    group_by(label) %>%
    nest() %>%
    mutate(data = map(data, function(x) return(x %>% arrange(quantile)))) %>%
    unnest() %>%
    # arrange(quantile) %>%
    ungroup() %>%
    pull(Quantiles)
  
  names(colorVec) <- l
  
  quantileData <- quantileData %>%
    mutate(Quantiles = factor(Quantiles, levels = (l))) %>%
    arrange(desc(quantile))
  
  
  labelMap <- levels(quantileData$genus)  
  names(labelMap) <- 1:length(labelMap)
  
  #print(levels(quantileData$Quantiles))
  print(length(unique(quantileData$label)))
  print(expectedNumLevels)
  if (length(unique(quantileData$label)) != expectedNumLevels) {
    asdaddads
  }
  
  groupLevels <- unique(quantileData$label)
  
  colorVec <- c(colorRampPalette(c(groupColors[1], 'white'))(length(quantiles)),
                colorRampPalette(c(groupColors[2], 'white'))(length(quantiles)))
  names(colorVec) <- l
  # print(head(quantileData))
  
  p <- ggplot(data = quantileData %>%
                mutate(Quantiles = factor(Quantiles, levels = l))) +
    geom_rect(data = quantileData %>%
                filter(label == groupLevels[1]) %>%
                mutate(Quantiles = factor(Quantiles, levels = l)),
              aes(xmin = as.integer(genus) -0.1 - 0.15, xmax = as.integer(genus) + 0.1 - 0.15, ymin = value_max, ymax = value_min, fill = Quantiles), color = 'black') +
    geom_rect(data = quantileData %>%
                filter(label == groupLevels[2]) %>%
                mutate(Quantiles = factor(Quantiles, levels = l)),
              aes(xmin = as.integer(genus) -0.1 + 0.15, xmax = as.integer(genus) + 0.1 + 0.15, ymin = value_max, ymax = value_min, fill = Quantiles), color = 'black') +
    geom_point(data = quantileData %>%
                 filter(label == groupLevels[1]) %>%
                 filter(str_detect(Quantiles, '50'))%>%
                 mutate(Quantiles = factor(Quantiles, levels = l)),
               aes(x = as.integer(genus)-0.15, y = value_max), fill = evaluated.color, size = 3, pch = 23, color = 'white') +
    geom_point(data = quantileData %>%
                 filter(label == groupLevels[2]) %>%
                 filter(str_detect(Quantiles, '50'))%>%
                 mutate(Quantiles = factor(Quantiles, levels = l)),
               aes(x = as.integer(genus)+0.15, y = value_max), fill = evaluated.color, size = 3, pch = 23, color = 'white') +            
    scale_x_continuous(breaks = as.integer(names(labelMap)),  labels = labelMap) +
    scale_fill_manual(values = colorVec[!str_detect(names(colorVec), '50')], breaks = names(colorVec[!str_detect(names(colorVec), '50')])) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    #coord_flip() +
    facet_wrap(~plotGroup, nrow = 2, scales = 'free') +
    xlab(xlab) +
    ylab(ylab) +
    NULL 
  return(p)
}

get_quantile_data <- function(i, cont = 'relAb', ...) {
  args <- as.character(list(...))
  quantileData <- list()
  #SSS <- inputData
  # inputData <- inputData %>%
  #  filter(type == 'reference (CV)')
  # medians <- i %>%
  # group_by(...) %>%
  # summarize(median = median(relAb + pc))
  for (quantile in quantiles) {
    tmp <- i %>%
      #select(cyl, wt) %>% 
      group_by(across(all_of(args))) %>%
      #mutate(cyl) %>%
      summarize(value_max = quantile(.data[[cont]]+1E-5, probs = quantile),
                value_min = quantile(.data[[cont]]+1E-5, probs = 1-quantile)) %>%
      mutate(quantile = quantile)
    quantileData[[length(quantileData) + 1]] <- tmp
  }
  var <- args[1]
  quantileData <- do.call('rbind', quantileData) %>%
    #left_join(medians, by = c('genus', 'label')) %>%
    mutate({{var}} := as.factor(.data[[var]])) %>%
    arrange(desc(quantile)) %>%
    mutate(Quantiles = map2_chr(quantile, .data[[args[length(args)]]], function(x, y) return(str_c((1-x) * 100, '% - ', x * 100, '% - ', y))))
  return(quantileData)
}


get_cancerness <- function(siam, pro, meta,
                           nudge_x = 0.5, pc = 1E-5, singleSample,
                           evaluated.color = '#21520e',
                           evaluated.label = 'New') {

  require(tidyverse)
  require(SIAMCAT)
  require(ggrepel)
  require(scales)
  require(purrr)
  
  # Recompute (CV) predictions from siamcat model	
  siam <- make.predictions(siam, siamcat.holdout = NULL)
  
  # Generate siamcat profile from profiles and meta
  label <- create.label(meta=meta,
                        label='label',
                        case='CRC',
                        control = 'Control')
  
  
  siamcatTest <- siamcat(feat=pro,
                         meta=meta,
                         label =label,
                         case='CRC')
  
  print(mean(rownames(pro) %in% rownames(siam@filt_feat[[1]])))
  
  siamcatTest <- normalize.features(siamcatTest,
                                    norm.param=norm_params(siam),
                                    feature.type='original',
                                    verbose = 2)
  
  # Predict
  siamcatTest <- make.predictions(siam, siamcatTest)
  
  testPredictions <- siamcatTest
  
  CVModelPredictions <- siam
  
  predMatrix <- testPredictions@pred_matrix %>%
    as.data.frame() %>%
    rownames_to_column('sampleID') %>%
    as_tibble() %>%
    filter(sampleID != 'dummySample') %>%
    pivot_longer(-sampleID) %>%
    group_by(sampleID) %>%
    summarize(medianPredictionProb = median(value)) %>%
    mutate(label = 'CRC') %>%
    mutate(type = 'evaluated') %>%
    mutate(alpha = 1)
  
  cvPredMatrix <- CVModelPredictions@pred_matrix %>%
    as.data.frame() %>%
    rownames_to_column('sampleID') %>%
    as_tibble() %>%
    filter(sampleID != 'dummySample') %>%
    pivot_longer(-sampleID) %>%
    group_by(sampleID) %>%
    summarize(medianPredictionProb = median(value)) %>%
    left_join(CVModelPredictions@label$label %>% as.data.frame() %>% dplyr::rename(label = '.') %>% rownames_to_column('sampleID'), by = 'sampleID') %>%
    mutate(label = ifelse(label == '1', 'CRC', 'Control')) %>%
    mutate(type = 'reference (CV)') %>%
    mutate(alpha = 0.3)
  
  dataFin <- rbind(predMatrix, cvPredMatrix)
  qTPR <- c(0.75, 0.5, 0.25)
  TPRVals <- map_dbl(qTPR, function(x) dataFin %>%
                       filter(type == 'reference (CV)') %>%
                       filter(label == 'CRC') %>%
                       dplyr::pull(medianPredictionProb) %>%
                       quantile(probs= x))
  #print(TPRVals)
  fracsTPR <- map_dbl(TPRVals, function(x) mean((dataFin %>%
                                                   filter(type == 'reference (CV)') %>%
                                                   filter(label == 'CRC') %>%
                                                   dplyr::pull(medianPredictionProb)) > x))
  #print(fracsTPR)
  
  dataFin <- rbind(predMatrix, cvPredMatrix)
  qFPR <- c(0.90, 0.95, 0.98, 0.99)
  FPRVals <- map_dbl(qFPR, function(x) dataFin %>%
                       filter(type == 'reference (CV)') %>%
                       filter(label == 'Control') %>%
                       dplyr::pull(medianPredictionProb) %>%
                       quantile(probs = x))
  
  
  fracsTPR <- map_dbl(FPRVals, function(x) mean((dataFin %>%
                                                   filter(type == 'reference (CV)') %>%
                                                   filter(label == 'CRC') %>%
                                                   dplyr::pull(medianPredictionProb)) > x))
  
  #print(FPRVals)
  
  dataFin_FPR <- lapply(FPRVals, function(threshold) {
    above_threshold <- dataFin %>% 
      filter(type=='reference (CV)') %>% 
      filter(label=='CRC') %>% 
      filter(medianPredictionProb > threshold)
    
    total_samples <- nrow(dataFin %>%filter(type=='reference (CV)') %>% 
                            filter(label=='CRC') )
    above_threshold_samples <- nrow(above_threshold)
    
    percentage_above_threshold <- (above_threshold_samples / total_samples) * 100
    
    return(data.frame(Threshold = threshold, PercentageAbove = percentage_above_threshold, totalSample=above_threshold_samples))
  })
  
  
  threshold_names <- c('Threshold 10%', 'Threshold 5%', 'Threshold 2%', 'Threshold 1%')
  dataFin_FPR <- bind_rows(dataFin_FPR) %>%
    mutate(Threshold = factor(Threshold, levels = FPRVals, labels = threshold_names))
  
  dataFin_new <- lapply(FPRVals, function(threshold) {
    above_threshold <- dataFin %>%
      filter(type=='evaluated') %>% 
      filter(label=='CRC') %>% 
      filter(medianPredictionProb > threshold)
    
    total_samples <- nrow(dataFin %>% filter(type=='evaluated') %>% 
                            filter(label=='CRC'))
    above_threshold_samples <- nrow(above_threshold)
    
    percentage_above_threshold <- (above_threshold_samples / total_samples) * 100
    
    return(data.frame(Threshold = threshold, PercentageAbove = percentage_above_threshold, totalSample=above_threshold_samples))
  })

  dataFin_new <- bind_rows(dataFin_new) %>%
    mutate(Threshold = factor(Threshold, levels = FPRVals, labels = threshold_names))
  

  
  dataFin <- dataFin %>%
      mutate(label = ifelse(type == 'evaluated', evaluated.label, label)) %>%
      mutate(label = factor(label, levels = c('CRC', evaluated.label, 'Control'))) %>% 
      mutate(label2=label)
    
  
  eval_data <- dataFin %>% filter(type=='evaluated') %>%  
    mutate(dataset=evaluated.label)  
  
  
  result_all_eval <- get_quantile_data(
    eval_data,
    cont = 'medianPredictionProb',
    args = c('dataset')
  )
  
  # Define quantile levels
  quantile_levels <- c('0% – 100%', '10% – 90%', '20% – 80%', '30% – 70%', '40% – 60%', '50% – 50%')
  
  # Generate gradient using the provided function and color '#FFBA08'
  generate_gradient <- function(color) {
    (colorRampPalette(c('white', color))(length(quantile_levels)))
  }
  
  result_all_eval <- result_all_eval %>%
    mutate(QuantileGroup = case_when(
      quantile == 1    ~ '0% – 100%',
      quantile == 0.9  ~ '10% – 90%',
      quantile == 0.8  ~ '20% – 80%',
      quantile == 0.7  ~ '30% – 70%',
      quantile == 0.6  ~ '40% – 60%',
      quantile == 0.5  ~ '50% – 50%',
      TRUE ~ NA_character_
    ))
  
  eval_gradient <- (generate_gradient(`evaluated.color`))
  
  # Assign fill colors to each quantile group
  result_all_eval <- result_all_eval %>%
    mutate(
      QuantileGroup = factor(QuantileGroup, levels = quantile_levels),
      FillColor = eval_gradient[as.integer(QuantileGroup)]
    )
  
  # Calculate median point for the diamond marker
  median_point_eval<- result_all_eval %>% 
    filter(QuantileGroup == '50% – 50%')

  # Plotting
  
  p1 <- ggplot() +
    # DRAW RECTANGLES FIRST (background)
    geom_rect(
      data = result_all_eval %>% arrange(desc(quantile)),
      aes(
        xmin = 2 - 0.32,
        xmax = 2 + 0.32,
        ymin = value_min,
        ymax = value_max,
        fill = FillColor
      ),
      color = 'black',
      linewidth = 0.5
    ) +
    
    # Draw lines on the rect
    geom_segment(
      data = result_all_eval,
      aes(x = 2 - 0.2, xend = 2 + 0.2, y = value_min, yend = value_min),
      color = 'black',
      linewidth = 0.3
    ) +
    geom_segment(
      data = result_all_eval,
      aes(x = 2 - 0.2, xend = 2 + 0.2, y = value_max, yend = value_max),
      color = 'black',
      linewidth = 0.3
    ) +
    
    # THEN POINTS SO THEY SHOW ON TOP
    geom_point(
      data = dataFin,
      aes(x = label, y = medianPredictionProb, color = label),
      size = 1.5,
      alpha = 0.4,
      position = position_jitternormal(sd_x = 0.1, sd_y = 0),
      show.legend = FALSE
    ) +
    geom_point(
      data = median_point_eval,
      aes(x = `evaluated.label`, y = value_max),
      shape = 23,
      size = 3.5,
      fill = 'white',
      color = 'black',
      stroke = 1
    ) +
    
    # TPR/FPR lines and text
    geom_segment(
      data = data.frame(label = 'CRC', TPRVals = TPRVals),
      aes(x = 1 - 0.5, xend = 1 + 0.5, y = TPRVals, yend = TPRVals),
      alpha = 0.8
    ) +
    geom_text(
      data = data.frame(label = 'CRC', TPRVals = rev(TPRVals), TPRLabel = qTPR * 100),
      aes(x = 1 - 0.5, y = TPRVals + 0.025, label = str_c(TPRLabel, '%'))
    ) +
    annotate(geom = 'text', x = 0.5, y = 0.95, label = 'TPR', fontface = 2) +
    geom_segment(
      data = data.frame(label = 'Control', FPRVals = FPRVals),
      aes(x = 3 - 0.5, xend = 3 + 0.5, y = FPRVals, yend = FPRVals),
      alpha = 0.8
    ) +
    geom_text(
      data = data.frame(label = 'Control', FPRVals = FPRVals, FPRLabel = (1 - qFPR) * 100),
      aes(x = 3 + 0.45, y = FPRVals + 0.025, label = str_c(FPRLabel, '%'))
    ) +
    annotate(geom = 'text', x = 3.5, y = 0.95, label = 'FPR', fontface = 2) +
    
    # scales, labels, themes
    scale_fill_identity() +
    scale_color_manual(
      values = c('CRC' = '#852020', 'Control' = 'dodgerblue4') %>%
        c(setNames(evaluated.color, evaluated.label))
    )+
    scale_x_discrete(limits = c('CRC', evaluated.label, 'Control')) +
    ylim(0, 1) +
    labs(
      y = 'Colorectal cancer microbiome signature score',
      x = NULL
    ) +
    theme_presentation() +
    theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    )
  
  
  
  
  return(list(p1,dataFin %>% select(c(sampleID, medianPredictionProb, label,   type  ))))
  


}
      
      
      
      
      
      
      
      
      
  
  