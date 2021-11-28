#' @title Predict Pseudotime value
#'
#' @description Predict pseudotime value of a dataset given a pseudotime model
#' @param test_meth a data frame of methylation data
#' @param test_meta a data frame of meta data
#' @param model_list a list of prcomp model object used in model training and
#' a slingshot object returned from CV_pseudotime or fit_pseudotime
#' @param tissue_type (optional) a string of tissue type for plot title
#' @param model_label (optional) a string of model label for plot file name
#' @param plot (optional) a boolean of whether to plot pseudotime trajectory
#' and pseudotime vs chronological age and save plots to disk
#' @param plot_trend (optional) a boolean of whether to fit and plot
#' nls best fit sum of two exponential trend line for pseudotime vs age
#' @param start_plot_fit (optional) a named vector of starting value for nls
#' fitting sum of two exponential trend to pseudotime vs age
#' @return a numeric vector of predicted pseudotime value
#' @export
#' @examples
#' \dontrun{predict_pseudotime(test_meth, test_meta,
#'                    list("reduced_dim_model" = prcomp_model,
#'                         "pstime_model" = slingshot_model),
#'                    tissue_type = "Blood",
#'                    model_label = "test_ver1",
#'                    plot = TRUE)}

predict_pseudotime <- function(test_meth,
                               test_meta,
                               model_list,
                               tissue_type = "",
                               model_label = "",
                               plot = TRUE,
                               plot_trend = FALSE,
                               start_plot_fit = c(a0=9.79,b1=2.0,c1=0.025,c2=0.412)){

  reduced_dim_model <- model_list[["reduced_dim_model"]]
  pstime_model <- model_list[["pstime_model"]]

  test_reduced_dim <- stats::predict(reduced_dim_model,
                                     t(test_meth))[,1:2]
  test_slg <- slingshot::predict(pstime_model,test_reduced_dim)
  test_pstime <- get_single_pseudotime(test_slg)
  if (plot == TRUE){

    if(length(unique(test_meta$exp_labels))>1){
      plot_scatter(x = test_reduced_dim[,1],
                   y = test_reduced_dim[,2],
                   title = "Principal Components",
                   group = test_meta$exp_labels,
                   tissue_type = tissue_type,
                   model_label = model_label,
                   xlab = "PC1", ylab = "PC2",
                   grouplab = "group",
                   filename = "PCs",
                   width = 2600,
                   height = 1800)

      plot_pseudotime_value(test_slg,
                            test_reduced_dim,
                            tissue_type = tissue_type,
                            model_label = model_label)

      if(plot_trend == TRUE){
        df <- data.frame(x = test_meta$chrono_age+0.5,
                         y = test_pstime)
        doubleExp <- nlsr::wrapnlsr(y ~ a0*((1+b1)-exp(-c1*x)-b1*exp(-c2*x)), data = df,
                                    start=start_plot_fit)

        plot_scatter(x = test_meta$chrono_age,
                     y = test_pstime,
                     title = "Pseudotime VS Chronological Age",
                     trend_model = doubleExp,
                     func_form = "doubleExp",
                     group = test_meta$exp_labels,
                     tissue_type = tissue_type,
                     model_label = model_label,
                     xlab = "Chronological Age (yrs)", ylab = "Pseudotime",
                     grouplab = "group",
                     filename = "chronoage",
                     width = 2600,
                     height = 1800)

      }
      else{
        plot_scatter(x = test_meta$chrono_age,
                     y = test_pstime,
                     title = "Pseudotime VS Chronological Age",
                     group = test_meta$exp_labels,
                     tissue_type = tissue_type,
                     model_label = model_label,
                     xlab = "Chronological Age (yrs)", ylab = "Pseudotime",
                     grouplab = "group",
                     filename = "chronoage",
                     width = 2600,
                     height = 1800)
      }


      plot_scatter(x = test_meta$epigen_age,
                   y = test_pstime,
                   title = "Pseudotime VS Epigenetic State",
                   group = test_meta$exp_labels,
                   tissue_type = tissue_type,
                   model_label = model_label,
                   xlab = "Epigenetic State", ylab = "Pseudotime",
                   grouplab = "group",
                   filename = "epigenage",
                   width = 2600,
                   height = 1800)

    }else{
      plot_scatter(x = test_reduced_dim[,1],
                   y = test_reduced_dim[,2],
                   title = "Principal Components",
                   #group = test_meta$exp_labels,
                   tissue_type = tissue_type,
                   model_label = model_label,
                   xlab = "PC1", ylab = "PC2",
                   #grouplab = "group",
                   filename = "PCs",
                   width = 2100,
                   height = 1800)

      plot_pseudotime_value(test_slg,
                            test_reduced_dim,
                            tissue_type = tissue_type,
                            model_label = model_label)
      if(plot_trend==TRUE){
        df <- data.frame(x = test_meta$chrono_age+0.5,
                         y = test_pstime)
        doubleExp <- nlsr::wrapnlsr(y ~ a0*((1+b1)-exp(-c1*x)-b1*exp(-c2*x)), data = df,
                                    start=start_plot_fit)

        plot_scatter(x = test_meta$chrono_age,
                     y = test_pstime,
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
        plot_scatter(x = test_meta$chrono_age,
                     y = test_pstime,
                     title = "Pseudotime VS Chronological Age",
                     tissue_type = tissue_type,
                     model_label = model_label,
                     xlab = "Chronological Age (yrs)", ylab = "Pseudotime",
                     filename = "chronoage",
                     width = 2400,
                     height = 1800)
      }


      plot_scatter(x = test_meta$epigen_age,
                   y = test_pstime,
                   title = "Pseudotime VS Epigenetic State",
                   #group = test_meta$exp_labels,
                   tissue_type = tissue_type,
                   model_label = model_label,
                   xlab = "Epigenetic State", ylab = "Pseudotime",
                   #grouplab = "group",
                   filename = "epigenage",
                   width = 2100,
                   height = 1800)
    }

  }


  return(test_pstime)

}

