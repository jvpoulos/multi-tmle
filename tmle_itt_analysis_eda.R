#################################################################
# Descriptive tables and plots for ITT analysis                 #
# (not longitudinal)                                            #
#################################################################

print(paste0("Summarizing results from output directory: ", output_dir))

if(est.LMTP){
  n.estimators <- 3
}else{
  n.estimators <- 2
}

blind.drugs <- TRUE
slides <- TRUE

print(paste0("Generating tables and plots for outcome:", outcome, "; Condition: ", condition))

outcome.labels <-list("combined"=c("diabetes diagnosis or death within 36 months"), 
                      "diabetes"=c("diabetes diagnosis within 36 months"),
                      "death"=c("death within 36 months"))

condition.labels <- list("none"=c("all patients"),
                         "schizophrenia"=c("patients with a Schizophrenia diagnosis at baseline"),
                         "mdd"=c("patients with a MDD diagnosis at baseline"),
                         "black"=c("black patients"),
                         "white"=c("white patients"),
                         "latino"=c("latino patients"))

blinded.labels <- c("Reference","A","B","C","D","E")

if(use.SL){
  estimator.labels <- c("Multinomial (super learner)", "Binomial (super learner)")
}else{
  estimator.labels <- c("Multinomial (Logistic)", "Binomial (Logistic)")
}

if(outcome=="combined" & condition=="none" & use.SL==TRUE){
  
  ## Density plot comparing individual probs. of treatment (multinomial- unbounded)
  
  print("unbounded treatment probs - multinomial")
  print(summary(g_preds))
  
  if(blind.drugs){
    colnames(g_preds) <- blinded.labels
  }else{
    colnames(g_preds) <- proper(colnames(g_preds))
  }
  
  treatment.probs.m <- melt(g_preds)
  
  propensity.plot.multinomial <- ggplot(treatment.probs.m,aes(x=value, fill=variable)) + 
    geom_density(alpha=0.35) +
    scale_x_continuous(breaks = seq(0,0.8,0.2)) +
    ylab("Density") +
    xlab(paste0("Estimated treatment probability"))+ 
    labs(fill = "Treatment") +
    theme(legend.position="bottom") +   theme(plot.title = element_text(hjust = 0.5, family="serif", size=16)) +
    theme(axis.title=element_text(family="serif", size=16)) +
    theme(axis.text.y=element_text(family="serif", size=14)) +
    theme(axis.text.x=element_text(family="serif", size=14)) +
    theme(legend.text=element_text(family="serif", size = 14)) +
    theme(legend.title=element_text(family="serif", size = 14)) +
    theme(strip.text.x = element_text(family="serif", size = 14)) +
    theme(strip.text.y = element_text(family="serif", size = 14)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l =0))) +
    theme(axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l =0))) +
    theme(panel.spacing = unit(1, "lines"))
  
  if(blind.drugs){
    if(slides){
      ggsave(paste0(output_dir, outcome, "_",condition, "_use_SL_",use.SL, "_propensity_plot_multinomial_blinded_slides.png"),propensity.plot.multinomial+ggtitle(paste0("Implementation: ", estimator.labels[1])), scale=1.25)
    } else{
      ggsave(paste0(output_dir, outcome, "_",condition, "_use_SL_",use.SL, "_propensity_plot_multinomial_blinded.png"),propensity.plot.multinomial, scale=1.25)
    }
  }
  ggsave(paste0(output_dir, outcome, "_",condition, "_use_SL_",use.SL, "_propensity_plot_multinomial.png"),propensity.plot.multinomial, scale=1.25)
  
  if(est.binomial==TRUE){
    # multiple binary treatment probs. 
    
    print("unbounded treatment probs - multiple binary")
    print(summary(g_preds_bin))
    
    if(blind.drugs){
      colnames(g_preds_bin) <- blinded.labels
    }else{
      colnames(g_preds_bin) <- proper(colnames(g_preds_bin))
    }
    
    treatment.probs.bin.m <- melt(g_preds_bin)
    propensity.plot.bin <- ggplot(treatment.probs.bin.m,aes(x=value, fill=Var2)) + 
      geom_density(alpha=0.35) +
      scale_x_continuous(breaks = seq(0,0.8,0.2)) +
      ylab("Density") +
      xlab(paste0("Estimated treatment probability"))+ 
      labs(fill = "Treatment") +
      theme(legend.position="bottom") +   theme(plot.title = element_text(hjust = 0.5, family="serif", size=16)) +
      theme(axis.title=element_text(family="serif", size=16)) +
      theme(axis.text.y=element_text(family="serif", size=14)) +
      theme(axis.text.x=element_text(family="serif", size=14)) +
      theme(legend.text=element_text(family="serif", size = 14)) +
      theme(legend.title=element_text(family="serif", size = 14)) +
      theme(strip.text.x = element_text(family="serif", size = 14)) +
      theme(strip.text.y = element_text(family="serif", size = 14)) +
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l =0))) +
      theme(axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l =0))) +
      theme(panel.spacing = unit(1, "lines"))
    
    if(blind.drugs){
      if(slides){
        ggsave(paste0(output_dir, outcome, "_",condition, "_use_SL_",use.SL, "_propensity_plot_bin_blinded_slides.png"),propensity.plot.bin+ggtitle(paste0("Implementation: ", estimator.labels[2])), scale=1.25)
      } else{
        ggsave(paste0(output_dir, outcome, "_",condition, "_use_SL_",use.SL, "_propensity_plot_bin_blinded.png"),propensity.plot.bin, scale=1.25)
      }
    }
    
    ggsave(paste0(output_dir, outcome, "_",condition, "_use_SL_",use.SL, "_propensity_plot_binomial.png"),propensity.plot.bin, scale=1.25)
  }
  
  ## Plot comparing estimated vs. observed treatment (overlap)
  obs.treatment <- dummify(A) # dummy matrix
  for(estimator in estimator.labels){
    for(j in 1:ncol(obs.treatment)){
      if(estimator %in% c("Multinomial (super learner)","Multinomial (Logistic)")){
        obs.treatment.probs <- g_preds*obs.treatment[,j]
      }else if (estimator %in% c("Binomial (super learner)","Binomial (Logistic)")){
        obs.treatment.probs <- g_preds_bin*obs.treatment[,j]
      }
      
      est.obs.treat <- data.frame("est"=c(c(t(obs.treatment.probs[rowSums(obs.treatment.probs[, -1])>0, ]))),
                                  "Estimator"=c(rep(estimator,length.out=length(c(t(obs.treatment.probs[rowSums(obs.treatment.probs[, -1])>0, ]))))),
                                  "Treatment"=c(rep(colnames(g_preds),length.out=length(c(t(obs.treatment.probs[rowSums(obs.treatment.probs[, -1])>0, ]))))))
      est.obs.treat.m  <- melt(est.obs.treat, id.vars=c("Estimator","Treatment"))
      
      ggplot(est.obs.treat.m[which(est.obs.treat.m$Estimator==estimator),],
             aes(x=value, fill=Treatment)) +
        geom_density(alpha=0.35) +
        xlab("Estimated treatment probability") + ylab("Density") + ggtitle(paste0("Patients assigned to ", colnames(g_preds)[j]))  +
        labs(fill = "Treatment") +
        theme(legend.position="bottom") +   theme(plot.title = element_text(hjust = 0.5, family="serif", size=16), plot.subtitle = element_text(hjust = 0.5, family="serif", size=12)) +
        theme(axis.title=element_text(family="serif", size=16)) +
        theme(axis.text.y=element_text(family="serif", size=14)) +
        theme(axis.text.x=element_text(family="serif", size=14)) +
        theme(legend.text=element_text(family="serif", size = 14)) +
        theme(legend.title=element_text(family="serif", size = 14)) +
        theme(strip.text.x = element_text(family="serif", size = 14)) +
        theme(strip.text.y = element_text(family="serif", size = 14)) +
        theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l =0))) +
        theme(axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l =0))) +
        theme(panel.spacing = unit(1, "lines"))
      
      if(blind.drugs){
        if(slides){
          ggsave(paste0(output_dir,"tmle_", outcome, "_", condition, "_use_SL_",use.SL, "_itt_overlap_",colnames(g_preds)[j],"_estimator_",estimator,"blinded_slides.png"), plot=last_plot()+ggtitle(paste0("Implementation: ", estimator), subtitle=paste0("Patients assigned to ", colnames(g_preds)[j])), scale=1.25)
        } else{
          ggsave(paste0(output_dir,"tmle_", outcome, "_", condition, "_use_SL_",use.SL, "_itt_overlap_",colnames(g_preds)[j],"_estimator_",estimator,"_blinded.png"), plot=last_plot(), scale=1.25)
        }
      }else{
        ggsave(paste0(output_dir,"tmle_", outcome, "_", condition, "_use_SL_",use.SL, "_itt_overlap_",colnames(g_preds)[j],"_estimator_",estimator,".png"), plot=last_plot(), scale=1.25)
      }
    }
  }
  
  ## Print share bounded
  print(paste0("% truncated, multinomial: ", colSums((g_preds > gbound[2]) | g_preds < gbound[1])/n))
  
  ## Print share bounded
  print(paste0("% truncated, multiple binary: ", colSums((g_preds_bin > gbound[2]) | g_preds_bin < gbound[1])/n))
  
  ## Print ESS
  print(paste0("ESS (multinomial, ATE): ", ess_ate_tmle))
  print(paste0("ESS (multiple binary, ATE): ", ess_ate_tmle_bin))
  
  if(est.LMTP==TRUE){
    print(paste0("% truncated, LMTP: ", colSums((do.call(cbind.data.frame, Ahat_lmtp) > gbound[2]) | do.call(cbind.data.frame, Ahat_lmtp) < gbound[1])/n))
    print(paste0("ESS (LMTP, ATE): ", ess_ate_lmtp))
  }
  
  ## Print multinomial preds
  print("Multinomial predictions")
  print(tableContinuous(data.frame(g_preds), cumsum=FALSE, stats= c("n", "min", "mean", "max", "s"), longtable = FALSE,prec =3))
  
  print("Multiple binary predictions")
  print(tableContinuous(data.frame(g_preds_bin), cumsum=FALSE, stats= c("n", "min", "mean", "max", "s"), longtable = FALSE,prec =3))
  
  if(est.LMTP==TRUE){
    print("LMTP predictions")
    print(tableContinuous(data.frame(do.call(cbind.data.frame, Ahat_lmtp)), cumsum=FALSE, stats= c("n", "min", "mean", "max", "s"), longtable = FALSE,prec =3))
  }
}

## Forest plot for ATEs

# Create data for plot

if(blind.drugs){
  comparisons <- blinded.labels[-1]
}else{
  comparisons <- proper(colnames(g_preds))[2:6]
}

ates.dat <- data.frame(x =rep(comparisons,n.estimators),
                       y = c(ATE_tmle,ATE_tmle_bin), 
                       y.lo = c(ATE_CI_tmle_lower,ATE_CI_tmle_bin_lower), 
                       y.hi = c(ATE_CI_tmle_upper,ATE_CI_tmle_bin_upper))

ates.dat$x <- as.factor(ates.dat$x)

ates.dat$Analysis <- c(rep("Multinomial (super learner)",each=length(comparisons)), rep("Binomial (super learner)",each=length(comparisons)))

saveRDS(ates.dat,paste0(output_dir,"tmle_",outcome,"_",outcome.type, "_use_SL_",use.SL,  "_", condition, "_", use.SL, "_","ates_dat.rds"))

if(use.SL==TRUE & condition!="none"){
  
  estimand <- "CATE"
  
  cates.dat <- data.frame(x =rep(comparisons,2),
                          y = c(CATE_tmle,CATE_tmle_bin),
                          y.lo = c(CATE_CI_tmle_lower,CATE_CI_tmle_bin_lower), 
                          y.hi = c(CATE_CI_tmle_upper,CATE_CI_tmle_bin_upper))
  
  cates.dat$x <- as.factor(cates.dat$x)
  
  cates.dat$Analysis <- c(rep("Multinomial (super learner)",each=length(comparisons)),rep("Binomial (super learner)",each=length(comparisons)))
  
  cate.plot <- ForestPlot(cates.dat,
                          xlab=estimand,ylab="Treatment drug") +
    labs(color="TMLE:") +
    scale_x_discrete(limits = rev) +
    theme(legend.position = "bottom") +
    theme(plot.title = element_text(hjust = 0.5, family="serif", size=16), plot.subtitle = element_text(hjust = 0.5, family="serif", size=12)) +
    theme(axis.title=element_text(family="serif", size=16)) +
    theme(axis.text.y=element_text(family="serif", size=14)) +
    theme(axis.text.x=element_text(family="serif", size=14,angle = 45, vjust = 0.5, hjust=1)) +
    theme(strip.text.x = element_text(family="serif", size = 14)) +
    theme(strip.text.y = element_text(family="serif", size = 14)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l =0))) +
    theme(axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l =0))) +
    theme(panel.spacing = unit(1, "lines"))
  
  cate.plot <- cate.plot +
    annotate("segment", x = 0.5, xend = 0.5, y = 0.002, yend = 0.025,
             colour = "blue", size = 1, arrow = arrow(length = unit(0.1, "inches"), ends="last", type="closed")) +
    annotate("segment", x = 0.5, xend = 0.5, y = -0.002, yend = -0.025,
             colour = "blue", size = 1, arrow = arrow(length = unit(0.1, "inches"), ends="last", type="closed")) +
    annotate("text", x = 0.65, y = -0.012, label = "Favors treatment drug")
  
  cate.plot <- cate.plot + annotate("text", x = 0.65, y = 0.012, label = "Favors reference drug") 

  if(slides){
    if(condition=="none"){
      cate.plot <- cate.plot + ggtitle(paste0("Outcome: ", outcome.labels[[outcome]]))
    }else{
      cate.plot <- cate.plot + ggtitle(paste0("Outcome: ", outcome.labels[[outcome]]), subtitle=paste0("Condition: ", condition.labels[[condition]]))
    }
  }
  if(blind.drugs){
    if(slides){
      ggsave(paste0(output_dir,"tmle_", outcome, "_",condition, "_use_SL_",use.SL, "_itt_CATE_blinded_slides.png"), plot=cate.plot, scale=1.25)
    }else{
      ggsave(paste0(output_dir,"tmle_", outcome, "_",condition, "_use_SL_",use.SL,"_itt_CATE_blinded.png"), plot=cate.plot, scale=1.25)
    }
  }else{
    if(slides){
      ggsave(paste0(output_dir,"tmle_", outcome, "_",condition, "_use_SL_",use.SL,"_itt_CATE_slides.png"), plot=cate.plot, scale=1.25)
    }else{
      ggsave(paste0(output_dir,"tmle_", outcome, "_",condition, "_use_SL_",use.SL,"_itt_CATE.png"), plot=cate.plot, scale=1.25)
    }
  }
  
}else{
  estimand <- "ATE"

  ate.plot <- ForestPlot(ates.dat,
                         xlab=estimand,ylab="Treatment drug") +
    labs(color="TMLE:") +
    scale_x_discrete(limits = rev) +
    theme(legend.position = "bottom") +
    theme(plot.title = element_text(hjust = 0.5, family="serif", size=16), plot.subtitle = element_text(hjust = 0.5, family="serif", size=12)) +
    theme(axis.title=element_text(family="serif", size=16)) +
    theme(axis.text.y=element_text(family="serif", size=14)) +
    theme(axis.text.x=element_text(family="serif", size=14,angle = 45, vjust = 0.5, hjust=1)) +
    theme(strip.text.x = element_text(family="serif", size = 14)) +
    theme(strip.text.y = element_text(family="serif", size = 14)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l =0))) +
    theme(axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l =0))) +
    theme(panel.spacing = unit(1, "lines"))
  
  ate.plot <- ate.plot +
    annotate("segment", x = 0.5, xend = 0.5, y = 0.002, yend = 0.025,
             colour = "blue", size = 1, arrow = arrow(length = unit(0.1, "inches"), ends="last", type="closed")) +
    annotate("segment", x = 0.5, xend = 0.5, y = -0.002, yend = -0.025,
             colour = "blue", size = 1, arrow = arrow(length = unit(0.1, "inches"), ends="last", type="closed")) +
    annotate("text", x = 0.65, y = -0.012, label = "Favors treatment drug")
  
  ate.plot <- ate.plot + annotate("text", x = 0.65, y = 0.012, label = "Favors reference drug") 
  
  if(slides){
    if(condition=="none"){
      ate.plot <- ate.plot + ggtitle(paste0("Outcome: ", outcome.labels[[outcome]]))
    }else{
      ate.plot <- ate.plot + ggtitle(paste0("Outcome: ", outcome.labels[[outcome]]), subtitle=paste0("Condition: ", condition.labels[[condition]]))
    }
  }
  
  if(blind.drugs){
    if(slides){
      ggsave(paste0(output_dir,"tmle_", outcome, "_",condition, "_est_binomial_",est.binomial,"_itt_ATE_blinded_slides.png"), plot=ate.plot, scale=1.25)
    }else{
      ggsave(paste0(output_dir,"tmle_", outcome, "_",condition,"_est_binomial_",est.binomial,"_itt_ATE_blinded.png"), plot=ate.plot, scale=1.25)
    }
  }else{
    if(slides){
      ggsave(paste0(output_dir,"tmle_", outcome, "_",condition, "_est_binomial_",est.binomial,"_itt_ATE_slides.png"), plot=ate.plot, scale=1.25)
    }else{
      ggsave(paste0(output_dir,"tmle_", outcome, "_",condition, "_est_binomial_",est.binomial,"_itt_ATE.png"), plot=ate.plot, scale=1.25)
    }
  }
}