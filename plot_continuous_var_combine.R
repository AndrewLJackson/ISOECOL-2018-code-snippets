#' Plot proportions by a continuous covariate
#'
#' \code{plot_continuous_var} creates a plot of how the mixture proportions
#' change according to a continuous covariate, as well as plots of the mixture
#' proportions for the individuals with minimum, median, and maximum covariate
#' values. Called by \code{\link{output_JAGS}} if any continuous effects are in
#' the model.
#'
#' MixSIAR fits a continuous covariate as a linear regression in ILR/transform-space.
#' Two terms are fit for the proportion of each source: an intercept and a slope.
#' The plotted line uses the posterior median estimates of the intercept and slope, and
#' the lines are curved because of the ILR-transform back into proportion-space. The 
#' 95\% credible intervals are shaded.
#' 
#' If the model contains both a continuous AND a categorical (factor) covariate, MixSIAR
#' fits a different intercept term for each factor level and all levels share the
#' same slope term.
#'
#' @param jags.1 output from \code{\link{run_model}}
#' @param mix output from \code{\link{load_mix_data}}
#' @param source output from \code{\link{load_source_data}}
#' @param output_options list containing options for plots and saving,
#' passed from \code{\link{output_JAGS}}
#'
#' @seealso Francis et al. 2011
plot_continuous_var_combine <- function(jags.1, mix, source, output_options){
  # added only to pass R CMD check
  # ilr.global <- x <- p.global <- p.ind <- sources <- ..scaled.. <- NULL
  R2jags::attach.jags(jags.1)
  n.sources <- source$n.sources
  source_names <- source$source_names
  
  for(ce in 1:mix$n.ce){
    
    if(mix$n.effects == 1){ # if there is a FE/RE in addition to continuous effect
      for(f1 in 1:mix$FAC[[1]]$levels){
        fac.lab <- mix$FAC[[1]]$labels[f1]
        label <- mix$cont_effects[ce]
        cont <- mix$CE[[ce]]
        ilr.cont <- get(paste("ilr.cont",ce,sep=""))
        
        get_high <- function(x){return(quantile(x,.975))}
        get_low <- function(x){return(quantile(x,.025))}
        n.plot = 200
        chain.len = dim(p.global)[1]
        Cont1.plot <- seq(from=round(min(cont),1), to=round(max(cont),1), length.out=n.plot)
        ilr.plot <- array(NA,dim=c(n.plot, n.sources-1, chain.len))
        ilr.median <- array(NA,dim=c(n.plot, n.sources-1))
        ilr.low <- array(NA,dim=c(n.plot, n.sources-1))
        ilr.high <- array(NA,dim=c(n.plot, n.sources-1))
        for(src in 1:n.sources-1){
          for(i in 1:n.plot){
            ilr.plot[i,src,] <- ilr.global[,src] + ilr.cont[,src]*Cont1.plot[i] + ilr.fac1[,f1,src]
            ilr.low[i,src] <- get_low(ilr.plot[i,src,])
            ilr.median[i,src] <- median(ilr.plot[i,src,])
            ilr.high[i,src] <- get_high(ilr.plot[i,src,])
          }
        }
        
        # Transform regression lines from ILR-space to p-space
        e <- matrix(rep(0,n.sources*(n.sources-1)),nrow=n.sources,ncol=(n.sources-1))
        for(i in 1:(n.sources-1)){
          e[,i] <- exp(c(rep(sqrt(1/(i*(i+1))),i),-sqrt(i/(i+1)),rep(0,n.sources-i-1)))
          e[,i] <- e[,i]/sum(e[,i])
        }
        cross.med <- array(data=NA,dim=c(n.plot, n.sources, n.sources-1))  # dummy variable for inverse ILR calculation
        tmp.p.med <- array(data=NA,dim=c(n.plot, n.sources))              # dummy variable for inverse ILR calculation
        p.median <- array(data=NA,dim=c(n.plot, n.sources))
        cross.low <- array(data=NA,dim=c(n.plot, n.sources, n.sources-1))  # dummy variable for inverse ILR calculation
        tmp.p.low <- array(data=NA,dim=c(n.plot, n.sources))              # dummy variable for inverse ILR calculation
        p.low <- array(data=NA,dim=c(n.plot, n.sources))
        cross.high <- array(data=NA,dim=c(n.plot, n.sources, n.sources-1))  # dummy variable for inverse ILR calculation
        tmp.p.high <- array(data=NA,dim=c(n.plot, n.sources))              # dummy variable for inverse ILR calculation
        p.high <- array(data=NA,dim=c(n.plot, n.sources))
        for(i in 1:n.plot){
          for(j in 1:(n.sources-1)){
            cross.med[i,,j] <- (e[,j]^ilr.median[i,j])/sum(e[,j]^ilr.median[i,j]);
            cross.low[i,,j] <- (e[,j]^ilr.low[i,j])/sum(e[,j]^ilr.low[i,j]);
            cross.high[i,,j] <- (e[,j]^ilr.high[i,j])/sum(e[,j]^ilr.high[i,j]);
          }
          for(src in 1:n.sources){
            tmp.p.med[i,src] <- prod(cross.med[i,src,]);
            tmp.p.low[i,src] <- prod(cross.low[i,src,]);
            tmp.p.high[i,src] <- prod(cross.high[i,src,]);
          }
          for(src in 1:n.sources){
            p.median[i,src] <- tmp.p.med[i,src]/sum(tmp.p.med[i,]);
            p.low[i,src] <- tmp.p.low[i,src]/sum(tmp.p.low[i,]);
            p.high[i,src] <- tmp.p.high[i,src]/sum(tmp.p.high[i,]);
          }
        }
        colnames(p.median) <- source_names
        
        Cont1.plot <- Cont1.plot*mix$CE_scale + mix$CE_center # transform Cont1.plot (x-axis) back to the original scale
        df <- data.frame(reshape2::melt(p.median)[,2:3], rep(Cont1.plot,n.sources), reshape2::melt(p.low)[,3], reshape2::melt(p.high)[,3])
        colnames(df) <- c("source","median","x","low","high")
        
        # medians <- data.frame(cont,apply(p.ind,c(2,3),median))
        # colnames(medians) <- c("cont",source_names)
        # medians <- melt(medians,id="cont")
        
        # Plot of Diet vs. Cont effect
        # Page 370 in Francis et al (2011)
        dev.new()
        print(ggplot2::ggplot(data=df,ggplot2::aes(x=x,y=median)) +
                ggplot2::geom_line(ggplot2::aes(x=x, y=median,group=source,colour=source),size=1.5) +
                ggplot2::geom_ribbon(ggplot2::aes(ymin=low, ymax=high, group=source, fill=source), alpha=0.35) +
                ggplot2::labs(title = fac.lab) +
                ggplot2::ylab("Diet Proportion") +
                ggplot2::xlab(label) +
                ggplot2::scale_y_continuous(expand = c(0, 0), limits=c(0,1)) +
                ggplot2::theme_bw() +
                ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(), 
                               panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), 
                               axis.line = ggplot2::element_line(colour = "black"), axis.title=ggplot2::element_text(size=16), 
                               axis.text=ggplot2::element_text(size=14), legend.text=ggplot2::element_text(size=14), legend.position=c(.02,1), 
                               legend.justification=c(0,1), legend.title=ggplot2::element_blank()))
        

      }
    } # end if YES FE/RE
    
    if(mix$n.effects == 0){
      label <- mix$cont_effects[ce]
      cont <- mix$CE[[ce]]
      ilr.cont <- get(paste("ilr.cont",ce,sep=""))
      
      get_high <- function(x){return(quantile(x,.975))}
      get_low <- function(x){return(quantile(x,.025))}
      n.plot = 200
      chain.len = dim(p.global)[1]
      Cont1.plot <- seq(from=round(min(cont),1), to=round(max(cont),1), length.out=n.plot)
      ilr.plot <- array(NA,dim=c(n.plot, n.sources-1, chain.len))
      ilr.median <- array(NA,dim=c(n.plot, n.sources-1))
      ilr.low <- array(NA,dim=c(n.plot, n.sources-1))
      ilr.high <- array(NA,dim=c(n.plot, n.sources-1))
      for(src in 1:n.sources-1){
        for(i in 1:n.plot){
          ilr.plot[i,src,] <- ilr.global[,src] + ilr.cont[,src]*Cont1.plot[i]
          ilr.low[i,src] <- get_low(ilr.plot[i,src,])
          ilr.median[i,src] <- median(ilr.plot[i,src,])
          ilr.high[i,src] <- get_high(ilr.plot[i,src,])
        }
      }
      
      # Transform regression lines from ILR-space to p-space
      e <- matrix(rep(0,n.sources*(n.sources-1)),nrow=n.sources,ncol=(n.sources-1))
      for(i in 1:(n.sources-1)){
        e[,i] <- exp(c(rep(sqrt(1/(i*(i+1))),i),-sqrt(i/(i+1)),rep(0,n.sources-i-1)))
        e[,i] <- e[,i]/sum(e[,i])
      }
      cross.med <- array(data=NA,dim=c(n.plot, n.sources, n.sources-1))  # dummy variable for inverse ILR calculation
      tmp.p.med <- array(data=NA,dim=c(n.plot, n.sources))              # dummy variable for inverse ILR calculation
      p.median <- array(data=NA,dim=c(n.plot, n.sources))
      cross.low <- array(data=NA,dim=c(n.plot, n.sources, n.sources-1))  # dummy variable for inverse ILR calculation
      tmp.p.low <- array(data=NA,dim=c(n.plot, n.sources))              # dummy variable for inverse ILR calculation
      p.low <- array(data=NA,dim=c(n.plot, n.sources))
      cross.high <- array(data=NA,dim=c(n.plot, n.sources, n.sources-1))  # dummy variable for inverse ILR calculation
      tmp.p.high <- array(data=NA,dim=c(n.plot, n.sources))              # dummy variable for inverse ILR calculation
      p.high <- array(data=NA,dim=c(n.plot, n.sources))
      for(i in 1:n.plot){
        for(j in 1:(n.sources-1)){
          cross.med[i,,j] <- (e[,j]^ilr.median[i,j])/sum(e[,j]^ilr.median[i,j]);
          cross.low[i,,j] <- (e[,j]^ilr.low[i,j])/sum(e[,j]^ilr.low[i,j]);
          cross.high[i,,j] <- (e[,j]^ilr.high[i,j])/sum(e[,j]^ilr.high[i,j]);
        }
        for(src in 1:n.sources){
          tmp.p.med[i,src] <- prod(cross.med[i,src,]);
          tmp.p.low[i,src] <- prod(cross.low[i,src,]);
          tmp.p.high[i,src] <- prod(cross.high[i,src,]);
        }
        for(src in 1:n.sources){
          p.median[i,src] <- tmp.p.med[i,src]/sum(tmp.p.med[i,]);
          p.low[i,src] <- tmp.p.low[i,src]/sum(tmp.p.low[i,]);
          p.high[i,src] <- tmp.p.high[i,src]/sum(tmp.p.high[i,]);
        }
      }
      colnames(p.median) <- source_names
      
      Cont1.plot <- Cont1.plot*mix$CE_scale + mix$CE_center # transform Cont1.plot (x-axis) back to the original scale
      df <- data.frame(reshape2::melt(p.median)[,2:3], rep(Cont1.plot,n.sources), reshape2::melt(p.low)[,3], reshape2::melt(p.high)[,3])
      colnames(df) <- c("source","median","x","low","high")
      
      # medians <- data.frame(cont,apply(p.ind,c(2,3),median))
      # colnames(medians) <- c("cont",source_names)
      # medians <- melt(medians,id="cont")
      
      # Plot of Diet vs. Cont effect
      # Page 370 in Francis et al (2011)
      # dev.new()
      print(ggplot2::ggplot(data=df,ggplot2::aes(x=x,y=median)) +
              ggplot2::geom_line(ggplot2::aes(x=x, y=median,group=source,colour=source),size=1.5) +
              ggplot2::geom_ribbon(ggplot2::aes(ymin=low, ymax=high, group=source, fill=source), alpha=0.35) +
              ggplot2::ylab("Diet Proportion") +
              ggplot2::xlab(label) +
              ggplot2::scale_y_continuous(expand = c(0, 0), limits=c(0,1)) +
              ggplot2::theme_bw() +
              ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(), 
                             panel.grid.minor = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), 
                             axis.line = ggplot2::element_line(colour = "black"), axis.title=ggplot2::element_text(size=16), 
                             axis.text=ggplot2::element_text(size=14), legend.text=ggplot2::element_text(size=14), legend.position=c(.02,1), 
                             legend.justification=c(0,1), legend.title=ggplot2::element_blank()))
      
     
    } # end if NO FE/RE
    
    
    
  } # end loop over ce
} #end plot_continuous_var function
