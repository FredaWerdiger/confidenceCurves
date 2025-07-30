#' @title Frequentist confidence analysis for any treatment effect
#'
#' @description This package performs frequentist confidence analysis, given a point estimate and associated error estimation, to answer the question:
#' How much confidence can we have in a particular treatment effect of interest?
#'
#' @details This is a package to perform frequentist confidence analysis on a observed treatment estimate. You may either supply a point estimate and associated precision estimate
#' via standard error, variance and sample size, and 95\% CI interval, or enter outcome data directly (the latter option is only available for binary data). Then, define a neutral effect,
#' and a meaningful clinical effect, and the direction of interest (above or below these) and the function will calculate how much confidence one can have in the associated treatment effect (e.g., beneficial, lacking meaningful benefit).
#' Also returned is the traditional frequentist p-value.
#'
#'
#' @references Equations used in this package are derived from Marschner, I. "Confidence distributions for treatment effects in Clinical Trials: Posteriors without Priors", Statistics in Medicine, 2024;43:1271-1289.

# input function:
#' @param theta.estimator Enter the point estimate, assumed to follow a Normal distribution.
#' @param treat.var Variance associated with the point estimate. Must be supplied with sample size.
#' @param standard.error Standard error associated with the point estimate.
#' @param confidence.lower Lower boundary of 95\% confidence interval.
#' @param confidence.upper Upper boundary of 95\% confidence interval.
#' @param num.ctrl When specifying binary input data; Number of subjects in the control group.
#' @param num.trmt When specifying binary input data; Number of subjects in the treatment group.
#' @param num.resp.ctrl When specifying binary input data; Number of responders in control group (who experienced the outcome).
#' @param num.resp.trmt When specifying binary input data; Number of responders in treatment group (who experienced the outcome).
#' @param sample.size Sample size. Can be calculated from num.ctrl and num.trmt.
#' @param neutral.effect Value corresponding to no effect. Default is 0. If using odds ratio, data should be on the log scale for a neutral effect of 0.
#' @param dir.benefit Direction (0 or 1) around the neutral effect corresponding to benefit. 0: less than no effect value (default); 1: more than no effect value.
#' @param directory Character string expressing directory where you want to save the confidence curve family. Default is "".
#' @param min.effect The minimally clinically interesting effect (meaningful benefit). Default is -0.05.
#' @param dir.min.effect Direction (0 or 1) around the min effect that you're interested in. 0: less than min effect value; 1: more than min effect value. Default assumes LACK of meaningful benefit.
#' @param equiv If interested in expressing confidence in treatment equivalence, you can specify
#'  two numbers as c(a,b) to bound the equivalency region. Default uses min.effect and -min.effect as a and b
#' @param show On the confidence density function, which region to display in shaded blue: BENEFIT' (default), 'LMB' (lack of meaningful benefit),
#' 'MB' (meaningful benefit) or 'EQUIV' (equivalence).
#' @param save.plot save the plot as png to directory, TRUE or FALSE (default).
#' @param return.plot Return the plots from the function, TRUE or FALSE (default).
#' @param tag phrase to append to the image filename as <directory>/confidence_curves_<tag>.png. Default is "".
#' @return List of quantities associated with confidence analysis (under $text) and (if supplied TRUE to return.plot) four confidence curves

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
                                 min.effect=-0.05,
                                 neutral.effect=0,
                                 dir.benefit=0,
                                 dir.min.effect=NULL,
                                 equiv=NULL,
                                 save.plot=FALSE,
                                 return.plot=FALSE,
                                 tag=""){

  if (is.null(theta.estimator)){

    if (!(is.numeric(num.resp.ctrl) &
          is.numeric(num.resp.trmt) &
          is.numeric(num.ctrl) &
          is.numeric(num.trmt))){

      stop("Specify either point estimate or binary outcome data.")

    }
    # ##############################################
    # CALCUATE TREATMENT EFFECT ESTIMATOR
    # ##############################################
    sample.size = num.ctrl + num.trmt

    if (!typeof(estimator.type) == "character"){
      stop("Enter estimator type (odds ratio, risk difference)")
    } else if (
      tolower(estimator.type) == 'odds ratio'|
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
      } else if (
        tolower(estimator.type) == "risk difference" |
        (grepl("risk", tolower(estimator.type)) &
         grepl("diff", tolower(estimator.type)))){
        estimator.type = "Risk Difference"
        crisk = num.resp.ctrl/num.ctrl
        trisk = num.resp.trmt/num.trmt

        risk_diff = trisk - crisk

        # variance after Marschner et al
        treat.var = (crisk * (1-crisk) +
                 (crisk+risk_diff)*(1 - crisk - risk_diff))/2

        standard.error = sqrt(treat.var / sample.size)

        theta.estimator = risk_diff
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
    conf.lmb = int_density(cd, min.effect, Inf)
  } else {
    conf.benefit = int_density(cd, neutral.effect, Inf)
    conf.lmb = int_density(cd, -Inf, min.effect)}

  # for neat display
  conf.benefit.disp = conf.disp(conf.benefit)
  conf.lmb.disp = conf.disp(conf.lmb)

  # equivalence
  # if thresholds are not specified then use min effect.
  if (!is.null(equiv)){
      conf.equiv = int_density(cd, equiv[1], equiv[2])
      conf.equiv.disp = conf.disp(conf.equiv)
  } else {
    conf.equiv = int_density(cd, min(min.effect, -min.effect),
                             max(min.effect, - min.effect))
    conf.equiv.disp = conf.disp(conf.equiv)
  }


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
  # arrows for benefit and lmb
  if (dir.benefit == 0){
    x.end.benefit = -standard.error/2
    if ((show == 'MB') | (dir.min.effect==dir.benefit)){
      x.end.lmb = min.effect - standard.error/2
    } else{
      x.end.lmb = min.effect + standard.error/2
    }
  } else {
    x.end.benefit = standard.error/2
    if ((show == 'MB') | (dir.min.effect==dir.benefit)){
      x.end.lmb = min.effect + standard.error/2
    } else {
      x.end.lmb = min.effect - standard.error/2
    }
  }

  # plot 1, cumulative distribution function
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
    # ggplot2::geom_hline(yintercept=0.95,  linetype="dashed") +
    # ggplot2::geom_hline(yintercept=0.75,  linetype="dashed") +
    # ggplot2::geom_hline(yintercept=0.50, linetype="dashed") +
    # ggplot2::geom_hline(yintercept=0.25, linetype="dashed") +
    # ggplot2::geom_hline(yintercept=0.05, linetype="dashed") +
    ggplot2::scale_x_continuous(n.breaks=5, expand = c(0, 0))+
    ggplot2::scale_y_continuous(n.breaks=6) +
    ggplot2::annotate(
      x = theta.estimator - (3*standard.error), y = 0.05, label = "5% CI", geom = "text",
      color = "black",
      lineheight = .3,
      vjust = -0.2,
      hjust=-0.1,
      size=3,
    ) +
    ggplot2::annotate(
      x = theta.estimator - (3*standard.error), y = 0.25, label = "25% CI", geom = "text",
      color = "black",
      lineheight = .3,
      vjust = 1.2,
      hjust=-0.1,
      size=3,
    ) +
    ggplot2::annotate(
      x = theta.estimator - (3*standard.error), y = 0.50, label = "50% CI", geom = "text",
      color = "black",
      lineheight = .3,
      hjust=-0.1,
      vjust = 1.2,
      size=3,
    ) +
    ggplot2::annotate(
      x = theta.estimator - (3*standard.error), y = 0.75, label = "75% CI", geom = "text",
      color = "black",
      lineheight = .3,
      hjust=-0.1,
      vjust = 1.2,
      size=3,
    ) +
    ggplot2::annotate(
      x = theta.estimator - (3*standard.error), y = 0.95, label = "95% CI", geom = "text",
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
    ggplot2::annotate(
      x = min.effect, y = -Inf, label = delta.label, geom = "text",
      parse = TRUE,
      color = "black",
      vjust = 1.1
    ) +
    ggplot2::geom_vline(xintercept=theta.estimator, linetype="dashed", linewidth=0.8,
               color="aquamarine2") +
    ggplot2::geom_vline(xintercept=neutral.effect, linetype="dashed", linewidth=0.8,
                        color="blue") +
    ggplot2::annotate('segment', x=neutral.effect, y=0.2, xend=x.end.benefit, yend=0.2,
                      color="blue", size=1.5, linejoin="mitre",lineend="butt",
                      arrow=ggplot2::arrow(length=ggplot2::unit(0.4, "cm"))) +
    ggplot2::geom_vline(xintercept=min.effect, linetype="dashed", linewidth=0.8,
                        color="forestgreen") +
    ggplot2::annotate('segment', x=min.effect, y=0.4, xend=x.end.lmb, yend=0.4,
                      color="forestgreen", size=1.5, linejoin="mitre", lineend="butt",
                      arrow=ggplot2::arrow(length=ggplot2::unit(0.4, "cm"))) +
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
    label = paste("Conf(REGION)=","\n", conf.benefit.disp,"%", sep='')
  } else if (show=='EQUIV'){
    dnorm_limit =function(x) {
      y = cd(x)
      y[x > max(equiv)] = NA
      y[x < min(equiv)] = NA
      return(y)
    }
    label = paste("Conf(REGION)=","\n", conf.equiv.disp,"%", sep='')
  } else if (show=='LMB'){
    dnorm_limit = function(x){
      y = cd(x)
      if (dir.benefit==0){
        y[x < min.effect] = NA
      } else {y[x > min.effect] = NA}
      return(y)
    }
    label = paste("Conf(REGION)=","\n", conf.lmb.disp,"%", sep='')
  } else if(show == 'MB'){
    dnorm_limit = function(x){
    y = cd(x)
    if (dir.benefit==0){
      y[x > min.effect] = NA
    } else {y[x < min.effect] = NA}
    return(y)
  }
  label = paste("Conf(REGION)=","\n", 100 - conf.lmb.disp,"%", sep='')
  }

  # moved dashed line accordingly
  if ((show == "BENEFIT") | (show == "EQUIV")){
    dashed.line.intercept=neutral.effect
  } else if ((show == "LMB") | (show == "MB")){
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
                        conf.meaningful.benefit = 1 - conf.lmb,
                        conf.equivalent = conf.equiv,
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
      print(list.out$cdf)
      }
    }
  }
  return(df)
}


