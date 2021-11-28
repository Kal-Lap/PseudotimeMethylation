#' @title Plot data points
#'
#' @description Make scatter plot of data points
#' @param x a numeric vector of x value
#' @param y a numeric vector of y value
#' @param title a string of plot title
#' @param group (optional) a vector of group labels of data
#' @param trend_model (optional) fitted (lm/nls) model of the trend line
#' @param func_form (optional) functional form of the trend line
#' @param tissue_type (optional) a string of tissue type for plot title
#' @param model_label (optional) a string of label for plot file name
#' @param xlab (optional) a string of x label
#' @param ylab (optional) a string of y label
#' @param grouplab (optional) a string of label for group labels
#' @param filename (optional) a string of plot file name
#' @param width (optional) a numeric value of width of the plot in pixel
#' @param height (optional) a numeric value of height of the plot in pixel

#' @export
#' @examples
#' \dontrun{plot_scatter(x = chrono_age,
#'              y = pseudotime_value,
#'              title = "Pseudotime VS Chronological Age",
#'              group = exp_labels,
#'              tissue_type = "Blood",
#'              model_label = "testv.1",
#'              xlab = "Chronological Age (yrs)", ylab = "Pseudotime",
#'              grouplab = "experiment",
#'              filename = "epigenage",
#'              width = 2600,
#'              height = 1350)}

plot_scatter <- function(x,y,
                         title,
                         group,
                         trend_model,
                         func_form = c("exp","quad","log","sqrt","doubleExp"),
                         tissue_type = "",
                         model_label = "",
                         xlab = "x", ylab = "y",
                         grouplab = "group",
                         filename = "scatter_plot",
                         width = 2100, height = 1500){
  if(tissue_type == ""){
    filename_tissue_type <- ""
    title_tissue_type <- ""
  }
  else{
    filename_tissue_type <- paste0(tissue_type,"_")
    title_tissue_type <- paste0(tissue_type,": ")
  }
  if(model_label != ""){
    model_label <- paste0("model_",model_label,"_")
  }

  theme <- ggplot2::theme(panel.border = ggplot2::element_blank(),
                          panel.background = ggplot2::element_blank(),
                          panel.grid.major = ggplot2::element_blank(),
                          panel.grid.minor = ggplot2::element_blank(),
                          axis.line = ggplot2::element_line(colour = "black"),
                          plot.title = ggplot2::element_text(hjust = 0.5,
                                                             face="bold"))


  if (missing(group)){
    plot_df <-data.frame(x, y)
    colnames(plot_df)<- c("x","y")
    aes <- ggplot2::aes_string(x="x", y = "y")

    p <- ggplot2::ggplot(plot_df,aes) +
      ggplot2::geom_point(colour = "grey") +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab) +
      #ggplot2::ggtitle(paste0(title_tissue_type,title)) +
      theme

  }else{
    plot_df <-data.frame(x, y, group)
    colnames(plot_df)<- c("x","y",grouplab)
    aes <- ggplot2::aes_string(x="x", y = "y", color = grouplab)

    legend_cols = ceiling(length(unique(group))/12)

    p <- ggplot2::ggplot(plot_df,aes) +
      ggplot2::geom_point() +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab) +
      #ggplot2::ggtitle(paste0(title_tissue_type,title)) +
      ggplot2::guides(color = ggplot2::guide_legend(ncol = legend_cols)) +
      theme
  }

  if (!missing(trend_model)){
    if(func_form == "exp"){
      trend_func <- function(x, coef){
        exp_term <- paste0("-exp(-",coef[-1],"*x)")
        eval(parse(text = paste0(coef[1],"*(",length(coef)-1, paste0(exp_term, collapse = ""),")")))
      }
      trend_label <- function(coef){
        exp_term <- paste0("-exp(-",round(coef[-1],digits=3),"x)")
        paste0(round(coef[1],digits=3),"(",length(coef)-1, paste0(exp_term, collapse = ""),")")
      }
    }

    if(func_form == "doubleExp"){
      trend_func <- function(x, coef){
        first_term <- paste0("-exp(-",coef[3],"*x)")
        second_term <- paste0("-",coef[2],"*exp(-",coef[4],"*x)")
        eval(parse(text = paste0(coef[1],"*((1+",coef[2],")",first_term,second_term,")")))
      }
      trend_label <- function(coef){
        first_term <- paste0("-exp(-",round(coef[3],digits=3),"x)")
        second_term <- paste0("-",round(coef[2],digits=3),"exp(-",round(coef[4],digits=3),"x)")
        paste0(round(coef[1],digits=3),"((1+",round(coef[2],digits=3),")",first_term,second_term,")")
      }
    }

    if(func_form == "quad"){
      trend_func <- function(x, coef){
        eval(parse(text = paste0(coef[1], "+", coef[2],"*x" ,"+", coef[3],"*x^2")))
      }
      trend_label <- function(coef){
        paste0(round(coef[1],digits = 3), "+", round(coef[2],digits = 3),"x","+", round(coef[3],digits = 3),"x^2")
      }
    }

    if(func_form == "linear"){
      trend_func <- function(x, coef){
        eval(parse(text = paste0(coef[1], "+", coef[2],"*x")))
      }
      trend_label <- function(coef){
        paste0(round(coef[1],digits = 3), "+", round(coef[2],digits = 3),"x")
      }
    }

    if(func_form == "log"){
      trend_func <- function(x, coef){
        eval(parse(text = paste0(coef[1], "+", coef[2],"*log(x)")))
      }
      trend_label <- function(coef){
        paste0(round(coef[1],digits = 3), "+", round(coef[2],digits = 3),"log(x)")
      }
    }

    if(func_form == "sqrt"){
      trend_func <- function(x, coef){
        eval(parse(text = paste0(coef[1], "+", coef[2],"*x^0.5")))
      }
      trend_label <- function(coef){
        paste0(round(coef[1],digits = 3), "+", round(coef[2],digits = 3),"x^0.5")
      }
    }


    p <- p + ggplot2::stat_function(mapping = ggplot2::aes_string(x="x"),
                                    fun = trend_func,
                                    args = list(coef= stats::coef(trend_model)),
                                    inherit.aes = FALSE,
                                    colour = "black",
                                    size = 1) +
              ggplot2::annotate("text",
                                x = min(plot_df$x),
                                y = max(plot_df$y),
                                label = trend_label(stats::coef(trend_model)),
                                hjust = 0,
                                size = 3.5)
  }

  plot_file_name <- paste0(filename_tissue_type, model_label,filename,".JPEG")
  ggplot2::ggsave(plot_file_name,
                  plot = p,
                  width = width,
                  height = height,
                  units = "px")
}

#' @title Plot pseudotime trajectory over cluster
#'
#' @description Plot pseudotime trajectory over clustered train data points
#' @param pstime_model a slingshot pseudotime object
#' @param reduced_dim a dataframe of training data in low dimensional space
#' @param clusters a vector of cluster labels of reduced_dim
#' @param tissue_type (optional) a string of tissue type for plot title
#' @param model_label (optional) a string of label for plot file name
#' @export
#' @examples
#' \dontrun{plot_pseudotime_value(slingshot_model,
#'                       reduced_dim,
#'                       clusters,
#'                       tissue_type = "Blood",
#'                       model_label = "train_v1")}

plot_pseudotime_cluster <- function(pstime_model,
                                    reduced_dim,
                                    clusters,
                                    tissue_type = "",
                                    model_label = ""){
  if(tissue_type == ""){
    filename_tissue_type <- tissue_type
    title_tissue_type <- tissue_type
  }
  else{
    filename_tissue_type <- paste0(tissue_type,"_")
    title_tissue_type <- paste0(tissue_type,": ")
  }

  if (model_label != ""){
    model_label <- paste0("model_",model_label,"_")
  }

  plot_file_name <- paste0(filename_tissue_type, model_label,"pstime_traj_clust.JPEG")

  grDevices::jpeg(plot_file_name,
                  width = 3300,
                  height = 2800,
                  res   = 450,
                  pointsize = 12)
  plot(reduced_dim,
       col = RColorBrewer::brewer.pal(9,"Set2")[clusters],
       pch=16,
       asp = 1,
       bty = "l",
       #main= paste0(title_tissue_type, "Pseudotime Trajectory"),
       xlab="PC1",
       ylab="PC2",
       #cex.main = 1.5,
       cex.lab = 1.2)
  graphics::lines(slingshot::SlingshotDataSet(pstime_model),
                  lwd = 3,
                  col = 'black')
  graphics::legend("topright",
                   legend = paste("Cluster",sort(unique(clusters))),
                   pch = 16,
                   bty="n",
                   col = RColorBrewer::brewer.pal(9,"Set2"))
  grDevices::dev.off()

}

#' @title Plot pseudotime trajectory with pseudotime value
#'
#' @description Plot pseudotime trajectory with pseudotime value of each data point
#' @param pstime_model a slingshot pseudotime object
#' @param reduced_dim a dataframe of training data in low dimensional space
#' @param tissue_type (optional) a string of tissue type for plot title
#' @param model_label (optional) a string of label for plot file name
#' @export
#' @examples
#' \dontrun{plot_pseudotime_value(slingshot_model,
#'                       reduced_dim,
#'                       tissue_type = "Blood",
#'                       model_label = "train_v1")}

plot_pseudotime_value <- function(pstime_model,
                                  reduced_dim,
                                  tissue_type = "",
                                  model_label = ""){
  if(tissue_type == ""){
    filename_tissue_type <- tissue_type
    title_tissue_type <- tissue_type
  }
  else{
    filename_tissue_type <- paste0(tissue_type,"_")
    title_tissue_type <- paste0(tissue_type,": ")
  }
  if (model_label != ""){
    model_label <- paste0("model_",model_label,"_")
  }

  plot_file_name <- paste0(filename_tissue_type, model_label,"pstime_traj_value.JPEG")
  pstime <- get_single_pseudotime(pstime_model)
  pal <- viridis::viridis(100, end = 0.95)
  colors <- pal[cut(pstime, breaks = 100)]


  grDevices::jpeg(plot_file_name,
                  width = 4000,
                  height = 2800,
                  res   = 450,
                  pointsize = 12)
  graphics::par(mar=c(6.1, 5.1, 4.1, 9.1))
  plot(reduced_dim,
       col = colors,
       pch=16,
       asp = 1,
       bty = "l",
       #main= paste0(title_tissue_type, "Pseudotime Trajectory"),
       xlab="PC1",
       ylab="PC2",
       #cex.main = 1.5,
       cex.lab = 1.2)
  graphics::lines(slingshot::SlingshotDataSet(pstime_model),
                  lwd = 3,
                  col = 'black')
  graphics::axis(side = 1)
  add_legend(pstime)
  grDevices::dev.off()
}

#' @title Add legend of pseudotime value
#'
#' @description Add continuous legend to pseudotime value plot
#' @param pstime a numeric vector of pseudotime value
#' @export
#' @examples
#' \dontrun{add_legend(c(32.3,45.3,23.5,1.2,55.6))}
add_legend <- function(pstime) {
  pal <- viridis::viridis(100, end = 0.95)
  lgd = rep(NA,length(pal))
  min = min(pstime)
  max = ceiling(max(pstime))
  mid = mean(c(min,max))
  lgd[c(1,50,100)] = c(min, mid,max)

  #graphics::par(mar=c(6.1, 5.1, 5.1, 9.1), xpd=TRUE)
  graphics::legend("topright",
                   #inset = c(0,0.04),
                   inset = c(-0.13,0.05),
                   legend = lgd,
                   fill = grDevices::colorRampPalette(colors = pal)(length(pal)),
                   border = NA,
                   y.intersp = 0.17,
                   bty = "n",
                   cex = 1,
                   text.font = 0.8,
                   xpd = TRUE)
}
