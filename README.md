# Confidence Curves
Constructing the set of Confidence Curves from observed data

## Overview
A function written in R to calculate the frequentist confidence distribution resulting from a point estimate and its associated error. 

## 'Posteriors without Priors'

 `"We are 99% confidence the treatment is beneficial"`

### What is a confidence distribution?
The confidence distribution is constructed by stacking every one-sided confidence interval (from 0-100) vertically to generate a cumulative distribution function. Figure 1 shows an example from real trial data where a point estimate of -0.22 (log odds ratio) was observed. On this scale, values less than zero represented benefit of a treatment to the patient. We see the stack of one-sided intervals, included the 99.89% confidence interval which intersects with zero - the value representing no effect. Any values within the 99.89% confidence interval represent benefit to the patient.

<figure>
  <p align="center">
      <img width="560" height="450" alt="confidence_distribution_INTERACT_annotated" src="https://github.com/user-attachments/assets/41d0d671-7832-439d-8cfd-b6d89d571bea" title="Example Confidence Distribution Function"/>
    <figcaption>Figure 1. Confidence distribution function from observed data that has a point estimate of -0.22 (a log odds ratio). The red dashed line represents the observed effect, and the blue line represents the line of no effect. Any treatment effect values less than zero represent benefit of treatment to the patient, and are in the <i>region of benefit</i>.</figcaption>
  </p>
</figure>

### Why would we do this?
By constructing this distribution, we have the flexibility to reason about a variety of treatment effects. For example, our treatment effect of interest in this example is "treatment benefit". Given that the 99.98% confidence interval contains all values that represent treatment benefit, we may say: ***We are 99.89% confident that the treatment has benefit***.

The traditional frequentist p-value does not allow us to make such a useful statement. Frequentist confidence statements are interpretable and intuitive and especially useful in the context of an adapative trial. 
