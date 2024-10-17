## Agenda
- comments and edits on our research proposal rough draft
- docs: https://docs.google.com/document/d/1AUAA0RJr9fe_u5RRp7WyGOJfFhPB-yHMRZGim7aQSOw/edit?usp=sharing

### Questions 
- would love some assitance in making the QIME code that we used in class specific to paired-end sequences
- denoising step: we can't seem to separate the forward and reverse reads for the denoise step
- does PICRUSt2 require the use of OTUs?
- questions about the README.md in our github: what should that look like?

## Meeting Minutes
- reviewed proposal
  - good title
  - introduction
    - could mention the functional potential of microbes that produce acidic or corrosive products affect instruments (with PICRUSt)
      - could be a sentence or a figure, depending on the results
    - The fungal part of the introduction may be out of the scope of our intro (1-2 sentences max, not an entire paragraph)
    - humidity affecting microbial diversity in general is fine
    - start very general and narrow down to more specific (eg. microbes, is, and then get more narrow)
    - the last paragraph should be an overview of what we have done and what we have found for the manuscript
  - hypothesis
    - biofilm formation is part of the acidic corrosion but may not need to be mentioned in the hypothesis
    - shift it to composition instead of diversity for the hypothesis and the research question (diversity is an aspect of microbial composition)
    - make sure the biofilm/microbial-induced corrosion part of the hypothesis is backed up by a biological basis
  - research aims
    - need to fit the hypothesis and build a story (should have reference papers that point you in that direction)\
    - we should talk about the biological basis to explain what we are doing (eg. "found that increased relative humidity increased   microbial-induced corrosive microbes, so we performed an indicator species analysis")
    - aim 1
      - acts as a control for confounding variable, if varying humidity shows varying results then perhaps we need to remove otherwise we can include, expecting there to be no changes between the two so that we can pool all of the samples together
    - aim 2
      -  can use Bray Curtis or Weighted Unifrac for beta diversity and Shannon for alpha diversity
        - either including phylogenetic relationship with the microbes or not
      -  should we should include a tax bar plot (not very quantitative) but a good idea to include
      - make sure to explain how this is different from aim 1 because we will look at alpha and beta diversity in both
    - aim 3
      - high, medium and low can be categorized based on the humidity percentages we have
      - could show bubble plots to show change
      - could also use indicator species and deseq2 which can show the relative change (fold change) in abundance
      - quick way: look at abundance and say it is different (20% vs 50%)
      - better way: reference frames paper, removes bias from sequencing, rob knight lab has a lot of papers on this, uses log10 to look at relative abundance changes
      - could do an indicator taxa analysis and could find some microbes of interest
        - first find certain microbes you are interested in and then perform deseq2
    - aim 4
      - look at things mentioned specifically in microbial-induced corrosion or biofilm formation or acidic metabolites
      - PICRUt predicts metabolic potential but it doesnt mean that it actually uses that metabolic pathway
      - limited by resolution of 16S sequencing
      - only predicting potential
  - data takes time to run so use a detached screen to do our analysis
