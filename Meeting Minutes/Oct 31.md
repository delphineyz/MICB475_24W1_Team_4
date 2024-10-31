## Agenda
- go over proposal feedback and plan next steps 
  - comment from marker:  "I am a bit concerned about such a high sample loss after rarefaction. It may be a good thing to discuss with your mentoring TA."
  - potentially think more about how we are doing the filtering; currently, are we only doing taxanomic filtering? Video in module 6 talks more about it and how, for example, if sample collection is different we could maybe filter differently; even though it isn't in our case, it'd be something to keep in mind.
  - also mentioned that we should avoid using ANCOM, even though we had Hans' support
  - TA commented that a T-test is wrong for AIM1 because our data is probably non-parametric, thus should we should consider Wilcoxon/Mann-Whiteney test right? (won't use Kruskall-wallis because we're comparing two groups, constant and then varying RH)
  - should we do our alpha/beta diversity tests on QIME, or on R?
  - for box 3F, how should we change it: "Generate bubble plot and stacked bar plot to show microbial taxa abundance" 
  - random notes
    - Jagroop, did you do feature-based filtering where you got rid of rare ASVs, maybe that is why we lost a lot of samples?   
    - **![image](https://github.com/user-attachments/assets/1b566b25-ab0e-4117-9358-b73a62785ebe)
**


## Meeting Minutes

- find some higher level taxa that are under or over-represented to throw into our hypothesis
- find a couple of metabolic pathways that can be over or under-represented with high humidity "Humidty has an effect on specific microbes which can lead to certain pathways being over or under-represented due to losing or gaining those specific microbes mentioned before."
- add the comments for the research question
- for aim 1, bring more evidence from other papers as well
- comment 14: use a taxa bar plot
- comment 15: we will repeat the same methodological section as aim 1 but then make sure to clarify how aim 1 and 2 differ again
- comment 16: could do indicator species analysis and then differential abundance or just differential abundance only
- comment 18: use ancom only, not ancom and deseq2
- comment 19: link how functional potential links to microbial composition
- rephrase wording to make more explicit about what constant and relative are and in aim 2 mention that we will look at either constant only or constant and relative
  - can give acronym: constant relative humidity (CRH) and varying relative humidity (VRH)
- change the "this will inform us on which samples" is too general instead write about how we will keep either constant only or constant and varying will be collapsed into the same range to increase our sample size
- write that we will either collapse into one or filter out varying RH that depends on aim 1 in box 2d and remove the taking out non microbial taxa
- talk about how we created low medium and high categories and list when we did it in the appraoch and we can have a data wrangling section for the appraoch where we talk about modifiy adding a new column for the RH and then 2a would be to filter out samples with varying humdity if needed if not continue (use the box names to our advantage)
- only mention as a taxa bar plot
- differential abundance is a pairwise comparision
- aldex2 use log ratio
- bar plot of log fold 2 change across the condiitons with a bar plot from ancom
- we should have a reference sample (such as 50 percent RH)
- determine what statistical anaylsis we would use for picrust
- qiime, picrust, dada2, ancom, any tool we use we need to cite them all
