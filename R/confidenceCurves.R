
# The equations in this function are derived from
# Marschner, I. "Confidence distributions for treatment effects in Clinical Trials:
# Posteriors without Priors", Statistics in Medicine, 2024;43:1271-1289.

# input function:
# theta.estimator = directly enter the point estimator on the log odds scale
# treat.var = variance of estimator
# standard.error = if variance is not supplied, with sample size, supply standard error or confidence interval
# confidence.upper = upper 95% confidence interval
# confidence.lower = lower 95% confidence interval
# num.ctrl = number of subjects in the control group
# num.trmt = number of subjects in the treatment group
# num.resp.ctrl = number of control subject that responded
# num.resp.trmt = number of treatment subjects that responded
# sample.size = enter if didn't supply num.ctrl and num.trmt
# neutral.effect = what represents a neutral effect (1, or 0)
# dir.benefit: 0 = lower is better, 1 = higher is better
# directory: existing location of where you want to save the image (end with /)
# show: Show 'BENEFIT' (default) or 'EQUIV' (equivalence) on graph


makeConfidenceCurves <- function(theta.estimator=NULL,
                                 estimator.type=NULL,
                                 treat.var=NULL,
                                 standard.error=NULL,
                                 confidence.upper=NULL,
                                 confidence.lower=NULL,
                                 sample.size=NULL,
                                 num.resp.ctrl=NULL,
                                 num.resp.trmt=NULL,
                                 num.ctrl=NULL,
                                 num.trmt=NULL,
                                 directory="",
                                 show='BENEFIT', pval='TWO-SIDED',
                                 min.effect=log(1.05),
                                 neutral.effect=0,
                                 dir.benefit=0,
                                 save.plot=FALSE,
                                 return.plot=FALSE,
                                 tag=""){

  if (is.null(theta.estimator)){
    if (is.null(sample.size) &
        !is.null(num.ctrl) & !is.null(num.trmt)){
      sample.size =  num.ctrl + num.trmt
    } else {
      stop("To get standard error from variance, specify sample size.")
    }
    # ##############################################
    # CALCUATE TREATMENT EFFECT ESTIMATOR
    # ##############################################

    if (!typeof(estimator.type) == "character"){
      stop("Enter estimator type (odds ratio, risk difference)")
    } else if (tolower(estimator.type) == 'odds ratio'|
          grepl("odds",tolower(estimator.type))){
      estimator.type = "Log Odds Ratio"
      a = num.resp.ctrl
      b = num.resp.trmt
      c = num.ctrl - num.resp.ctrl
      d = num.trmt - num.resp.trmt
      theta.estimator = log((a*d)/(b*c))
      # standard error after Marschner (2024)
      standard.error = sqrt(((1/a) + (1/b) + (1/c) + (1/d)))
      # variance
      treat.var = (standard.error**2) * sample.size
      } else if (to.lower(estimator.type) == "risk difference" |
                 (grepl("risk", tolower(estimator.type)) &
                  grepl("diff", tolower(estimator.type)))){
        estimator.type = "Risk Difference"
        crisk = num.resp.ctrl/num.ctrl
        trisk = num.resp.trmt/num.trmt

        theta.estimator = trisk - crisk

        # variance after Marschner et al
        var = (crisk * (1-crisk) +
                 (crisk+risk_diff)*(1 - crisk - risk_diff))/2
      } else {
      stop("Enter risk difference or odds ratio for estimator type.")
      }
  }

  ################
  # STANDARD ERROR
  ################

  # check inputs for standard error

  if (is.null(standard.error)){
    if (!is.null(confidence.upper) & !is.null(confidence.lower)){
      # else check if confidence interval was supplied
      standard.error = (confidence.upper - confidence.lower) / 3.92
    } else if (!is.null(treat.var)){
      if (is.null(sample.size)){
        if (!is.null(num.ctrl) & !is.null(num.trmt)){
          sample.size =  num.ctrl + num.trmt
        } else{
          stop("To get standard error from variance, specify sample size.")
        }
      } else{
        standard.error = sqrt(treat.var / sample.size)
      }
    } else{
      stop("Error estimation not supplied.")
    }
  }
  #########################
  # CONFIDENCE DISTRIBUTION
  #########################

  # now building the curves from (Marschner, 2024) from the theta estimator

  # confidence distribution function
  # from one-sided confidence intervals
  cdf = function(theta){
    temp = (theta.estimator - theta)/ (standard.error)
    return(1 - pnorm(temp)) # 1 - alpha = confidence
  }

  # confidence density function
  cd = function(theta){
    temp = (theta.estimator - theta)/ (standard.error)
    return((1/standard.error) * dnorm(temp))
  }

  # confidence curve
  # from two-sided confidence intervals
  # used to derive the p-value/type 1 error
  # it is also called "pvalue function"
  # and is equivalent to 2 x (1 - cdf())
  # where cdf() is the cumulative or confidence distribution function
  # assuming the distribution of theta under the null is symmetric about 0,

  cc = function(theta){
    temp = (theta.estimator - theta)/ (standard.error)
    temp = 1 - pnorm(temp)
    return(abs(1 - (2*temp)))
  }


  # two-sided p-value
  pval.func.two.tailed = function(theta){
    return(1-cc(theta)) }

  # one-sided p-value
  # sided-ness depends on the location of the estimator
  pval.func.one.tailed = function (theta, theta.hat){
    if (dir.benefit==0){
      if (theta.hat < 0){
        return(1-cdf(theta))
      } else {
        return(cdf(theta))
      }
    } else if (dir.benefit==1){
      if (theta.hat > 0){
        return(cdf(theta))
      } else {
        return(1-cdf(theta))
      }
    }
  }

  # ####################
  # CONSTRUCT THE CURVES
  # ####################

  # let the standard error determine the limits of the x-axis
  # and ensure it is centred on theta estimator

  x.min = theta.estimator - (4 * standard.error)
  x.max = theta.estimator + (4 * standard.error)

  # make sure that x.max goes past zero
  while (x.max < 0){
    x.min = x.min - standard.error
    x.max = x.max + standard.error
  }

  # setting the spacing of the axis ticks
  if (standard.error < 0.1){
    x.ticks = 0.001
  } else {x.ticks = 0.01}
  x = seq(x.min, x.max, x.ticks)  # x-axis

  #################################################
  # COMPUTING INTERCEPT POINTS FOR GRAPH ANNOTATION
  #################################################

  # For annotation of graphs

  # first compute the function for x
  cc.x = lapply(x,function(x) ceiling(cc(x) * 100))
  cdf.x = lapply(x, function(x) ceiling(cdf(x) * 100))

  # CDF: get intercepts for multiple confidence interval
  cdf.int.x = vector(mode="list", 5)
  ints = c(5, 25, 50, 75, 95)
  for (i in 1:length(ints)){
    temp = which.min(abs(as.numeric(cdf.x) - ints[i]))
    cdf.int.x[i] = x[min(temp)]
  }


  # CC: get the 95% CI
  # Closest(DescTools) allows us to approximate 95
  ind = DescTools::Closest(as.numeric(cc.x), 95, which=TRUE)

  # if we only get one intercept:
  if (length(ind)==1){
    # explicitly cut the x.axis in half and ensure one results from each half
    c1 = as.numeric(cc.x)[seq(1, ceiling(length(cc.x)/2))]
    c2 = as.numeric(cc.x)[seq(ceiling(length(cc.x)/2) + 1, length(cc.x))]
    ind.lower = DescTools::Closest(c1, 95, which=TRUE)
    ind.upper = DescTools::Closest(c2, 95, which=TRUE)  + ceiling(length(cc.x)/2)
    ind = c(ind.lower, ind.upper)
  }

  # take the first and last values of ind to get the two intercepts
  cc.min = signif(x[ind[1]], 2)
  cc.max = signif(x[tail(ind, 1)],2)


  # find the x-location nearest to the theta estimator
  loc = which.min(abs(as.numeric(x) - theta.estimator))

  ###########################
  # COMPUTE CONFIDENCE VALUES
  ###########################

  # for dir.benefit  = 0:
  # confidence in benefit: conf(theta<0)
  # confidence in futility: conf(theta > -lmb)

  # functions:

  int_density = function(fn, min, max){
    auc = integrate(fn, min, max)
    conf = auc$value
    return(conf)
  }

  conf.disp = function(conf){
    if (conf>0.9999 | conf<0.0001){
      return(round(conf, digits=6)*100)
    } else{
      return(round(conf, digits=4) * 100)
    }
  }

  # compute confidence
  # depending on direction of benefit and the position of neutral effect
  if (dir.benefit == 0){
    conf.benefit = int_density(cd, -Inf, neutral.effect)
    conf.lmb = int_density(cd, -min.effect, Inf)
  } else {
    conf.benefit = int_density(cd, neutral.effect, Inf)
    conf.lmb = int_density(cd, -Inf, min.effect)}

  # for neat display
  conf.benefit.disp = conf.disp(conf.benefit)
  conf.lmb.disp = conf.disp(conf.lmb)

  # equivalence (never been used)
  # (conf(0<theta<log(min.effect)))
  conf.equiv = int_density(cd, -min.effect, min.effect)
  conf.equiv.disp = conf.disp(conf.equiv)

  # p-values under null
  p.value.one.tailed = pval.func.one.tailed(0, theta.estimator)[[1]]
  p.value.two.tailed = pval.func.two.tailed(0)[[1]]
  if (p.value.one.tailed < 0.001){
    p.value.two.tailed = formatC(p.value.two.tailed, format='e', digits=4)
    p.value.one.tailed = formatC(p.value.one.tailed, format='e', digits=4)
  } else {
    p.value.two.tailed = signif(p.value.two.tailed, digits=2)
    p.value.one.tailed = signif(p.value.one.tailed, digits=2)
  }

  # one- or two-sided
  if(toupper(pval)=='TWO-SIDED'){
    p.value = p.value.two.tailed
  } else if (toupper(pval)=='ONE-SIDED'){
    p.value = p.value.one.tailed
  }

  ##############
  # CREATE PLOTS
  #############

  # as dataframe
  x = data.frame(x=x)

  # ----  ---- ----
  # CDF

  # put label in correct format
  # standard error is used to fix locations

  # estimator
  theta.label = list('paste(hat(theta))')
  # lmb
  delta.label = list('paste(delta)')
  # p-value
  label.p.one = data.frame(
    x = theta.estimator - (2*standard.error),
    y = 0.4,
    label = paste("p = ", p.value.one.tailed, sep='')
  )

  plot1 = ggplot2::ggplot(x, ggplot2::aes(x)) +
    ggplot2::xlab("Treatment effect") +
    ggplot2::ylab("One-sided Confidence Distribution Function") +
    ggplot2::stat_function(fun=function(x) cdf(x), linewidth=1.2) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_text(angle = 90),
          line) +
    ggplot2::annotate('segment', x=x.min, y=0.05, xend=as.numeric(cdf.int.x[1]), yend=0.05,
             arrow=ggplot2::arrow(length=ggplot2::unit(0.2, "cm"), type="closed", ends="last")) +
    ggplot2::annotate('segment', x=x.min, y=0.25, xend=as.numeric(cdf.int.x[2]), yend=0.25,
             arrow=ggplot2::arrow(length=ggplot2::unit(0.2, "cm"), type="closed", ends="last")) +
    ggplot2::annotate('segment', x=x.min, y=0.5, xend=as.numeric(cdf.int.x[3]), yend=0.5,
             arrow=ggplot2::arrow(length=ggplot2::unit(0.2, "cm"), type="closed", ends="last")) +
    ggplot2::annotate('segment', x=x.min, y=0.75, xend=as.numeric(cdf.int.x[4]), yend=0.75,
             arrow=ggplot2::arrow(length=ggplot2::unit(0.2, "cm"), type="closed", ends="last")) +
    ggplot2::annotate('segment', x=x.min, y=0.95, xend=as.numeric(cdf.int.x[5]), yend=0.95,
             arrow=ggplot2::arrow(length=ggplot2::unit(0.2, "cm"), type="closed", ends="last")) +
    ggplot2::geom_hline(yintercept=0.95,  linetype="dashed") +
    ggplot2::geom_hline(yintercept=0.75,  linetype="dashed") +
    ggplot2::geom_hline(yintercept=0.50, linetype="dashed") +
    ggplot2::geom_hline(yintercept=0.25, linetype="dashed") +
    ggplot2::geom_hline(yintercept=0.05, linetype="dashed") +
    ggplot2::scale_x_continuous(n.breaks=5, expand = c(0, 0))+
    ggplot2::scale_y_continuous(n.breaks=6) +
    ggplot2::annotate(
      x = theta.estimator + (2*standard.error), y = 0.05, label = "5% CI", geom = "text",
      color = "black",
      lineheight = .3,
      vjust = 1.2,
      hjust=-0.1,
      size=3,
    ) +
    ggplot2::annotate(
      x = theta.estimator+ (2*standard.error), y = 0.25, label = "25% CI", geom = "text",
      color = "black",
      lineheight = .3,
      vjust = 1.2,
      hjust=-0.1,
      size=3,
    ) +
    ggplot2::annotate(
      x = theta.estimator+ (2*standard.error), y = 0.50, label = "50% CI", geom = "text",
      color = "black",
      lineheight = .3,
      hjust=-0.1,
      vjust = 1.2,
      size=3,
    ) +
    ggplot2::annotate(
      x = theta.estimator+ (2*standard.error), y = 0.75, label = "75% CI", geom = "text",
      color = "black",
      lineheight = .3,
      hjust=-0.1,
      vjust = 1.2,
      size=3,
    ) +
    ggplot2::annotate(
      x = theta.estimator+ (2*standard.error), y = 0.95, label = "95% CI", geom = "text",
      color = "black",
      lineheight = .3,
      hjust=-0.1,
      vjust = 1.2,
      size=3,
    ) +
    ggplot2::annotate(
      x = theta.estimator, y = -Inf, label = theta.label, geom = "text",
      parse = TRUE,
      color = "black",
      vjust = 1.1
    ) +
    ggplot2::geom_vline(xintercept=theta.estimator, linetype="dashed", linewidth=0.8,
               color="aquamarine2") +
    ggplot2::geom_label(data=label.p.one, ggplot2::aes( x=x, y=y, label=label),
               color="black",
               fill = "blanchedalmond",
               lineheight = 0.9,
               size=4, label.padding = ggplot2::unit(0.5, "lines"),
               label.size=0) +
    ggplot2::coord_cartesian(clip = "off")

  # ----  ---- ----
  # PDF

  # for filling under the curve and writing confidence levels
  if (show=='BENEFIT'){
    dnorm_limit = function(x) {
      y = cd(x)
      if (dir.benefit==0){
        y[x > neutral.effect] = NA
      } else {y[x < neutral.effect] = NA}
      return(y)
    }
    label = paste("Conf(BENEFIT)=","\n", conf.benefit.disp,"%", sep='')
  } else if (show=='EQUIV'){
    dnorm_limit =function(x) {
      y = cd(x)
      y[x > min.effect & x< - min.effect] = NA
      return(y)
    }
    label = paste("Conf(EQUIV)=","\n", conf.equiv.disp,"%", sep='')
  } else if (show=='LMB'){
    dnorm_limit = function(x){
      y = cd(x)
      if (dir.benefit==0){
        y[x < -min.effect] = NA
      } else {y[x > min.effect] = NA}
      return(y)
    }
    label = paste("Conf(LMB)=","\n", conf.lmb.disp,"%", sep='')
  }

  # moved dashed line accordingly
  if ((show == "BENEFIT") | (show == "EQUIV")){
    dashed.line.intercept=neutral.effect
  } else if (show == "LMB"){
    dashed.line.intercept=min.effect
  }

  plot2 = ggplot2::ggplot(x, ggplot2::aes(x)) +
    ggplot2::xlab("Treatment effect") +
    ggplot2::ylab('Confidence Density Function') +
    ggplot2::stat_function(fun=function(x) cd(x), linewidth=1.2) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(angle = 90))+
    ggplot2::scale_x_continuous(n.breaks=5, expand = c(0, 0))+
    ggplot2::scale_y_continuous(n.breaks=6) +
    ggplot2::stat_function(fun = dnorm_limit, geom = "area", fill = "blue", alpha = 0.2)+
    ggplot2::geom_vline(xintercept=dashed.line.intercept, linetype="dashed", linewidth=0.6,
               color="darkgrey") +
    ggplot2::annotate(
      x = theta.estimator, y = cd(theta.estimator)/3, label=label, geom = "text",
      color = "black",
      lineheight = .9,
      # vjust = .8,
      size = 4,
    ) +
    ggplot2::annotate(
      x = theta.estimator, y = -Inf, label = theta.label, geom = "text",
      parse = TRUE,
      color = "black",
      vjust = 1.1
    ) +
    ggplot2::coord_cartesian(clip = "off")

  # ----  ---- ----
  # CC / p-val function

  # 95% CI label
  label.loc = data.frame(
    x = theta.estimator + (2*standard.error),
    y = 0.3,
    label = paste("95% CI:","\n","(",as.name(cc.min),",",as.name(cc.max),")", sep='')
  )
  # location of dots/points
  point.loc = data.frame(
    x=c(cc.min, cc.max),
    y=c(0.95, 0.95)
  )

  plot3 = ggplot2::ggplot(x, ggplot2::aes(x)) +
    ggplot2::stat_function(fun=function(x) cc(x), linewidth=1.2) +
    ggplot2::xlab("Treatment effect") +
    ggplot2::ylab('Two-sided Confidence Distribution Function')+
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(angle = 90))+
    ggplot2::scale_x_continuous(n.breaks=5, expand = c(0, 0))+
    ggplot2::scale_y_continuous(n.breaks=6) +
    ggplot2::annotate(
      x = theta.estimator, y = -Inf, label = theta.label, geom = "text",
      parse = TRUE,
      color = "black",
      vjust = 1.1
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::annotate('segment', x=cc.min, y=0.95, xend=cc.max, yend=0.95) +
    ggplot2::geom_point(data=point.loc, ggplot2::aes(x=x, y=y),
                        size=4, color="black", fill="orange", shape=21) +
    ggplot2::geom_label(data=label.loc, ggplot2::aes( x=x, y=y, label=label),
               color="black",
               fill = "orange",
               lineheight = 0.9,
               size=4,
               label.padding = ggplot2::unit(0.5, "lines"),
               label.size=0)

  # ----  ---- ----
  # NULL

  # get the z-statistic
  z.score = qnorm(conf.benefit)

  dnorm_limit_sd = function(x) {
    y = dnorm(x, sd=standard.error)  ## assumed same SD ##
    if (theta.estimator < 0){
      y[x > theta.estimator & x < -theta.estimator] = NA
    } else {
      y[x > -theta.estimator & x < theta.estimator] = NA
    }

    return(y)
  }

  # create new limits
  x.min.shift = x.min - theta.estimator
  x.max.shift = x.max - theta.estimator
  x = data.frame(x=seq(x.min.shift, x.max.shift, x.ticks))

  # format z score label
  label.z = data.frame(
    x = theta.estimator - 4*x.ticks, # set location automatically relative to theta estimator
    y = (dnorm(0, sd=standard.error)/2),
    label = paste("z = ", round(qnorm(conf.benefit), digits=2), sep='')
  )

  # format p-value label
  label.p =  paste("p (two-sided) = ", "\n", p.value.two.tailed, sep='')


  plot4 = ggplot2::ggplot(x, ggplot2::aes(x)) +
    ggplot2::xlab("Treatment effect under null") +
    ggplot2::ylab("Probability Density Function") +
    ggplot2::theme_bw() +
    ggplot2::stat_function(fun=function(x) dnorm(x, sd=standard.error), linewidth=1.2) +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(angle = 90))+
    ggplot2::scale_x_continuous(n.breaks=5, expand = c(0, 0))+
    ggplot2::scale_y_continuous(n.breaks=6) +
    ggplot2::annotate(
      x = theta.estimator, y = -Inf, label = theta.label, geom = "text",
      parse = TRUE,
      color = "black",
      vjust = 1.1
    ) +
    ggplot2::stat_function(fun = dnorm_limit_sd, geom = "area", fill = "blue4", alpha = 0.4)+
    ggplot2::geom_vline(xintercept=theta.estimator, linetype="dashed", linewidth=0.6,
               color="darkgrey") +
    ggplot2::geom_label(data=label.z, ggplot2::aes( x=x, y=y, label=label),
               color="black",
               fill = "orange",
               lineheight = 0.9,
               size=4,
               label.padding = ggplot2::unit(0.5, "lines"),
               angle=90,
               label.size=0)+
    ggplot2::annotate(
      x = 0, y = cd(theta.estimator)/3, label=label.p, geom = "text",
      color = "black",
      lineheight = .9,
      size = 4,
    )+
    ggplot2::coord_cartesian(clip = "off")

  # plot the grid
  cowplot::plot_grid(plot1, plot2, plot3,plot4, labels="AUTO" )

  ############
  # SAVE PLOTS
  ############

  if(save.plot=='TRUE'){
    # set the directory
    dir.create(file.path(directory), showWarnings = FALSE, recursive = TRUE)
    # setwd(file.path(directory))
    if ((! endsWith(directory, "/")) & (! endsWith(directory, "\\"))){
      directory = paste0(directory, "/")
    }
    # save the grid
    ggplot2::ggsave(paste(directory, "confidence_curves_",tag,".png", sep=''),
           dpi=300, device="png", bg="white",
           height=8, width=10, units="in")
  }

  ################################
  # RETURN LIST OF RELEVANT VALUES
  ################################

  p.type = 'two-sided'
  if(pval=='ONE-SIDED'){p.type = 'one-sided'
  }

  return.results = list(min.meaningful.effect=min.effect,
                        mean=round(theta.estimator, digits=4),
                        s.error = round(standard.error, digits=4),
                        conf.benefit = conf.benefit,
                        conf.lack.meaningful.benefit = conf.lmb,
                        p.value=p.value,
                        p.value.test=p.type,
                        ninetyfive.percent.CI.lower=cc.min,
                        ninetyfive.percent.CI.upper=cc.max

  )
  if (return.plot == TRUE){
    return(list(text=return.results, cdf=plot1, pdf=plot2, cc=plot3, null=plot4))
  } else{
    return(return.results)
  }
}

testConfidenceCurves <- function(num.ctrl=50,
                                 num.trmt=50,
                                 vary.ctrl=seq(16,20, by=2),
                                 vary.trmt=seq(26, 30, by=2),
                                 vary.lmb = c(-0.05, -0.1),
                                 estimate.type = 'odds ratio',
                                 dir.benefit = 0,
                                 neutral.effect = 0,
                                 directory='./test'){
  df <- data.frame()
  for (i in vary.lmb){
    for (j in vary.trmt){
      for (k in vary.ctrl){
      list.out <- confidenceCurves::makeConfidenceCurves(num.resp.ctrl = k,
                                       num.resp.trmt = j,
                                       num.ctrl = num.ctrl,
                                       num.trmt = num.trmt,
                                       estimator.type = 'odds ratio',
                                       min.effect = i,
                                       directory = directory,
                                       pval = 'ONE-SIDED',
                                       show='BENEFIT',
                                       save.plot=FALSE,
                                       return.plot=TRUE,
                                       dir.benefit = dir.benefit,
                                       neutral.effect=neutral.effect,
                                       tag=paste("delta",i,"treatresp",j,"ctrlresp", k, sep="" )
      )
      df <- rbind(df, data.frame(list.out$text))
      print(list.out$pdf)
      }
    }
  }
  return(df)
}


