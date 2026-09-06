Yes, we have some edits to make because the figures shown in the paper are based on running this code and reporting uncorrected P values. The other issue was that an FDR correction was too strict.

We had some questions for changes to implement and the Wood lab provided the following answers below. What we need is the code re run (can be altered to be Matlab or Python if easier) and we need the following adjustments made based on the Wood lab's feedback/ notes to our questions:
For the Luminex 32 cytokine panel:
What is the appropriate threshold for deciding which of the analytes to show on the heatmap and use in further comparisons? Should all 32 be shown?
Prune out targeted based on linear range criteria
Any focus on particular analytes? - ones that are significant vs pulling out myeloid related ones
PLSDA - cannot do stats on linked variable - show bar or violin plots instead for individual analytes - remove p values on the violin plots by heat map
Permutation analysis is statistically sound following PLSDA (permutation approach from Kristie's papers)
For 5xFAD/ WT cytokines - batch adjustment for sex after the linear range
Are there sex-specific differences? Can show batched heat map if needed
when examining individual analytes for significant differences between the 3 flicker groups (No Ficker, 20Hz, 40Hz), how should these be corrected for multiple comparisons? (N= 8 animals per group).
ANOVA and multiple comparisons correction (Tukeys) and if the data is not not parametric  use KW and Dunns Post Hoc Corrections
Alternatively: Multiple effects model for bar plots and for latent variables as well to control for sex ?
For Spearman correlations between an analyte and amyloid beta load (For example amyloid I-40 correlated with M-CSF). The statistical analysis is Spearman, the correlation line is Pearson (linear regression trend line)
What is the appropriate top and bottom correlation (Rho value) to show a relationship between the two? (we also have p values)
Use the p values
How should this be corrected for multiple comparisons? (for example: FDR based on the number of analytes, so 32 comparisons?)
Correlation matrices FDR across everything compared - is probably too much - limit to pre-selected cytokines (from the initial heat map possibly)
Alternatively, FDR adjusted by row (so 32 corrections)
Alternatively, picked out based on Myeloid and control for those comparisons - controlling for specific comparisons limits to just those selected (e.g. the 5 myeloid related cytokines not all 32)
For Spearman correlations between an analyte and a morphology parameter, a morphology parameter (For example MIG correlated with Number of Microglia)
What is the appropriate top and bottom correlation (Rho value) to show a relationship between the two?
In biology .46 is okay, but the p value is the driver for what is included
Statistical analysis is Spearman, the correlation line is Pearson (linear regression trend line)
Labeled Males and Females on scatter plots




Ashley

-----------------------------------
Ashley Prichard, PhD
Research Health Scientist (RR&D)
Durham VA Medical Center: Research Division
US Department of Veterans Affairs
Singer Lab Affiliate, Georgia Institute of Technology| Emory University

Email: aprichard3@gatech.edu | ashley.prichard@va.gov
Website: www.ashleyprichard.com

From: Xu, William O <williamxu@gatech.edu>
Sent: Saturday, July 11, 2026 1:13 AM
To: Prichard, Ashley M <aprichard3@gatech.edu>
Subject: Re: Interest in processing data using R?
 
Did a review of the code in the box you sent me and was able to get an understanding. I think I can make it work. Besides the outline you sent me in the original email, are there any more specifics you guys are trying to look at with the updated analysis? 

William

Get Outlook for Mac
From: Xu, William O <williamxu@gatech.edu>
Date: Tuesday, July 7, 2026 at 3:47 PM
To: Prichard, Ashley M <aprichard3@gatech.edu>
Subject: Re: Interest in processing data using R?

Thanks for sending, I should be able to make sense of it, will just take a bit to overview the code + manuscript. Will get back to you. 

William

Get Outlook for Mac
From: Prichard, Ashley M <aprichard3@gatech.edu>
Date: Tuesday, July 7, 2026 at 3:33 PM
To: Xu, William O <williamxu@gatech.edu>
Subject: Re: Interest in processing data using R?

Hi William,
Attached is a draft of the manuscript. The R code and prior analyses of the cytokine data is in the Sara Temp folder https://gatech.box.com/s/jumxuptgrfujpi6aakqq7x2429vwd7qw
Do you want to take a look at it and see if you are still interested/ can make it work? The code is not super commented, and as I mentioned, I am not a super fan of "R".

Ashley

-----------------------------------
Ashley Prichard, PhD
Research Health Scientist (RR&D)
Durham VA Medical Center: Research Division
US Department of Veterans Affairs
Singer Lab Affiliate, Georgia Institute of Technology| Emory University

Email: aprichard3@gatech.edu | ashley.prichard@va.gov
Website: www.ashleyprichard.com

From: Xu, William O <williamxu@gatech.edu>
Sent: Friday, July 3, 2026 4:56 PM
To: Prichard, Ashley M <aprichard3@gatech.edu>
Subject: Re: Interest in processing data using R?
 
Yes, I would be interested, whatever you can send for context would be greatly appreciated! Manuscript + Code would be an amazing place to start off! 

William 

Get Outlook for Mac
From: Prichard, Ashley M <aprichard3@gatech.edu>
Date: Thursday, July 2, 2026 at 4:34 PM
To: Xu, William O <williamxu@gatech.edu>
Cc: Bales, Katie Leigh <klbales@emory.edu>
Subject: Interest in processing data using R?

Hi William,
I wanted to gauge your interest in helping with a manuscript to process cytokine data that is time sensitive? The manuscript is in preparation (written draft), we have morphology data post flicker, and we want to revise how the cytokine data was processed for revised figures before submitting the paper at the end of the summer. 

For background, the data is from a study using hippocampus samples to analyze post flicker effects on 32 different cytokines. Mice were in 12 groups (the same as in the RNA data) across Sex, Genotype, and Flicker Frequencies (12 groups altogether).

Code already exists from the Wood lab in R to process the data, however Sara has completed her PhD and is working on transitioning to a new job etc. and does not have the time to re-run the analyses with the modifications we need.  I do not like R and have little experience with it myself. The analysis does not have to be in R per-se, but the statistical comparisons between the groups need to be within Genotype and possibly separated by sex, plus FDR corrected for multiple comparisons. We then wanted to correlate the relevant cytokines identified with the changes in cell morphology within subject. 

Let me know if you are interested or not, and we can chat further. I can also send more information, example of the manuscript and code etc. if you want to see what you'd be working with.

Ashley

-----------------------------------
Ashley Prichard, PhD
Research Health Scientist (RR&D)
Durham VA Medical Center: Research Division
US Department of Veterans Affairs
Singer Lab Affiliate, Georgia Institute of Technology| Emory University

Email: aprichard3@gatech.edu | ashley.prichard@va.gov
Website: www.ashleyprichard.com
