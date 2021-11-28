#' @title Fit Pseudotime Model
#'
#' @description Fit slingshot pseudotime model on training data
#' @param train_meth a data frame of training methylation data
#' @param train_meta a data frame of meta data for the training samples
#' @param clust_num (optional) a numeric value of cluster number
#' @param tissue_type (optional) a string of tissue type for plot title
#' @param model_label (optional) a string of model label for plot file name
#' @param choose_young_clust (optional) a boolean whether to use cluster with
#' youngest mean age as starting cluster of the trajectory
#' @param plot (optional) a boolean of whether to plot pseudotime trajectory
#' and pseudotime vs chronological age and save plots to disk
#' @param plot_trend (optional) a boolean of whether to fit and plot
#' nls best fit sum of two exponential trend line for pseudotime vs age
#' @param start_plot_fit (optional) a named vector of starting value for nls
#' fitting sum of two exponential trend to pseudotime vs age

#' @return a list of 2
#' \describe{
#'   \itemize{
#'     \item "reduced_dim_model" prcomp model object used in model training
#'     \item "pstime_model" trained slingshot model object
#'   }
#' }
#' @export
#' @examples
#' \dontrun{
#' fit_pseudotime(train_meth, train_meta, tissue_type = "Blood",
#'                model_label = "train_ver1", plot = TRUE)
#' }

fit_pseudotime <- function(train_meth,
                           train_meta,
                           clust_num,
                           tissue_type = "",
                           model_label = "",
                           choose_young_clust = FALSE,
                           plot = TRUE,
                           plot_trend = FALSE,
                           start_plot_fit = c(a0=9.79,b1=2.0,c1=0.025,c2=0.412)){

  reduced_dim_model <- stats::prcomp(t(train_meth), scale. = FALSE)
  reduced_dim <- reduced_dim_model$x[,1:2]
  if(missing(clust_num)){
    clust_model <- mclust::Mclust(reduced_dim)
    clusters <- clust_model$classification
  }else{
    clust_model <- mclust::Mclust(reduced_dim, G = clust_num)
    clusters <- clust_model$classification
  }


  #calculate mean age of each cluster
  if (choose_young_clust == TRUE){
    mean_clust_age <- vector()
    clust_labels <- paste(unique(clusters))
    for (i in clust_labels){
      curr_clust <- names(clusters)[clusters == i]
      mean_clust_age[i] <- mean(subset(train_meta, train_meta$sample_labels %in% curr_clust)$chrono_age)
    }
    youngest_clust <- names(which.min(mean_clust_age))
    oldest_clust <- names(which.max(mean_clust_age))

    #fit model
    pstime_model <- slingshot::slingshot(reduced_dim, clusters,
                                         start.clus = youngest_clust,
                                         end.clus = oldest_clust)
    #fit linage
    lin_num <- length(pstime_model@metadata$lineages)
    if (lin_num != 1){
      pstime_model <- slingshot::slingshot(reduced_dim, clusters,
                                           start.clus = youngest_clust)
      lin_num <- length(pstime_model@metadata$lineages)
      if (lin_num != 1){
        message("Warning: There are more than 1 linage")
      }
    }
  }else{
    #fit model
    pstime_model <- slingshot::slingshot(reduced_dim, clusters)
    lin_num <- length(pstime_model@metadata$lineages)
    if (lin_num != 1){
      message("Warning: There are more than 1 linage")
    }
  }

  #pseudotime value
  pstime <- get_single_pseudotime(pstime_model)

  if (plot == TRUE){
    #plot data points in reduced dimension
    if (length(unique(train_meta$exp_labels))>1){
      plot_scatter(x = reduced_dim[,1],
                   y = reduced_dim[,2],
                   group = train_meta$exp_labels,
                   tissue_type = tissue_type,
                   model_label = model_label,
                   xlab = "PC1", ylab = "PC2",
                   title = "Principal Components",
                   filename = "PCs",
                   grouplab = "group",
                   width = 2100,
                   height = 1500)

      #plot trajectory on clustered dimension
      plot_pseudotime_cluster(pstime_model,
                              reduced_dim,
                              clusters,
                              tissue_type,
                              model_label)
      #plot pseudotime trajectory
      plot_pseudotime_value(pstime_model,
                            reduced_dim,
                            tissue_type,
                            model_label)

      if(plot_trend == TRUE){
        df <- data.frame(x = train_meta$chrono_age+0.5,
                         y = pstime)
        doubleExp <- nlsr::wrapnlsr(y ~ a0*((1+b1)-exp(-c1*x)-b1*exp(-c2*x)), data = df,
                                    start=start_plot_fit)
        plot_scatter(x = train_meta$chrono_age,
                     y = pstime,
                     title = "Pseudotime VS Chronological Age",
                     trend_model = doubleExp,
                     func_form = "doubleExp",
                     group = train_meta$exp_labels,
                     tissue_type = tissue_type,
                     model_label = model_label,
                     xlab = "Chronological Age (yrs)", ylab = "Pseudotime",
                     grouplab = "group",
                     filename = "chronoage",
                     width = 2600,
                     height = 1800)
      }else{
        plot_scatter(x = train_meta$chrono_age,
                     y = pstime,
                     title = "Pseudotime VS Chronological Age",
                     group = train_meta$exp_labels,
                     tissue_type = tissue_type,
                     model_label = model_label,
                     xlab = "Chronological Age (yrs)", ylab = "Pseudotime",
                     grouplab = "group",
                     filename = "chronoage",
                     width = 2600,
                     height = 1800)
      }

    }
    else{
      plot_scatter(x = reduced_dim[,1],
                   y = reduced_dim[,2],
                   #group = train_meta$exp_labels,
                   tissue_type = tissue_type,
                   model_label = model_label,
                   xlab = "PC1", ylab = "PC2",
                   title = "Principal Components",
                   filename = "PCs",
                   #grouplab = "group",
                   width = 2100,
                   height = 1500)

      #plot trajectory on clustered dimension
      plot_pseudotime_cluster(pstime_model,
                              reduced_dim,
                              clusters,
                              tissue_type,
                              model_label)
      #plot pseudotime trajectory
      plot_pseudotime_value(pstime_model,
                            reduced_dim,
                            tissue_type,
                            model_label)

      if(plot_trend == TRUE){
        df <- data.frame(x = train_meta$chrono_age+0.5,
                         y = pstime)
        doubleExp <- nlsr::wrapnlsr(y ~ a0*((1+b1)-exp(-c1*x)-b1*exp(-c2*x)), data = df,
                                    start=start_plot_fit)
        plot_scatter(x = train_meta$chrono_age,
                   y = pstime,
                   title = "Pseudotime VS Chronological Age",
                   trend_model = doubleExp,
                   func_form = "doubleExp",
                   tissue_type = tissue_type,
                   model_label = model_label,
                   xlab = "Chronological Age (yrs)", ylab = "Pseudotime",
                   filename = "chronoage",
                   width = 2400,
                   height = 1800)
      }else{
        plot_scatter(x = train_meta$chrono_age,
                     y = pstime,
                     title = "Pseudotime VS Chronological Age",
                     tissue_type = tissue_type,
                     model_label = model_label,
                     xlab = "Chronological Age (yrs)", ylab = "Pseudotime",
                     filename = "chronoage",
                     width = 2400,
                     height = 1800)
      }
    }
  }

  return(list("reduced_dim_model" = reduced_dim_model,
              "pstime_model" = pstime_model))
}

#' @title Get single pseudotime value
#'
#' @description Get single pseudotime value for each data points.
#' If there are more than 1 linages, the pseudotime is calculated
#' as a weighted average across all linage.
#' @param pseudotime_obj trained slingshot model object
#' @return an array of pseudotime value for each data points
#' @export
#' @examples
#' \dontrun{get_single_pseudotime(slingshot_model)}
get_single_pseudotime <- function(pseudotime_obj){
  ps_time <- slingshot::slingPseudotime(pseudotime_obj)
  weight <- slingshot::slingCurveWeights(pseudotime_obj)
  denominator <- ceiling(rowSums(weight))
  weighted_ps_time <- ps_time*weight
  weighted_ps_time[is.na(weighted_ps_time)] <- 0
  weighted_avg <- rowSums(weighted_ps_time)/ceiling(rowSums(weight))
}

#' @title Mean absolute error
#'
#' @description Calculate mean absolute error
#' @param y a numeric vector of predicted data
#' @param x a numeric vector of observed data
#' @return Mean absolute error between x and y
#' @export
#' @examples
#' \dontrun{mae(c(1,2,3),c(4,5,6))}
mae <- function(y,x){
  return(sum(abs(x-y))/length(x))
}

#' @title Root mean square error
#'
#' @description Calculate root mean square error
#' @param y a numeric vector of predicted data
#' @param x a numeric vector of observed data
#' @return root mean square error between x and y
#' @export
#' @examples
#' \dontrun{rmse(c(1,2,3),c(4,5,6))}
rmse <- function(y,x){
  return(sqrt(sum((x-y)**2)/length(x)))
}

#' @title R squared
#'
#' @description Calculate R squared value
#' @param preds a numeric vector of predicted data
#' @param actual a numeric vector of observed data
#' @return R squared value between preds and actual
#' @export
#' @examples
#' \dontrun{r2(c(1,2,3),c(4,5,6))}
r2 <- function(preds,actual){
  return(1- sum((preds - actual) ** 2)/sum((actual - mean(actual))**2))
}
