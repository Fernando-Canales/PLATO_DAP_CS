# Changes Checklist
## Purpose of this document
To have an easy-to-look checklist for the main corrections of the paper.

### So far:

- All editorials have been corrected
- All minor comments have been corrected
- Major Point 7 needs to be corrected because
	- I need to re do Section 8, now with the results from the Variable transit scenario
### What should be inside the Response Letter?
Excerpt from David Moulliet's email about the referee report:
>	In your cover letter, please list all the changes made (including new or removed coauthors), and mark the changes in your manuscript using boldface or colored text. However, please note that your source files should only include the "clean" revised version.  A separate pdf with the marked changes (bold or diff file) should be uploaded as an "optional file."

The Response letter is the cover letter. There I will put a reply to all the Major, Minor and Editorial points made by the referee. In each case I would also indicate the respective change in the manuscript.
### Here I talk about code stuff

- [x] Major 1: Add EB occurrence rate of 1%
- [x] Major 2: Add distribution of transit depths and durations

## Here I talk about paper/latex stuff

### Week of the 14th of July 2025

- [x] [[Major_Comments_Response#First Major point|Major 1]]: State in the latex document that is not our aim to estimate a realistic EB rate for PLATO. Such thing has be done already (see Bray et al. 2023). I can also mention the Kepler and TESS papers from Prsa. 
Change the following plots and the corresponding figure captions:
- [x] Nominal mask
- [x] Secondary mask
- [x] Extended mask
- [x] Efficiency plot for the double-aperture photometry (baseline scenario)
- [x] Efficiency plot for the centroid shift (baseline scenario)
- [x] Efficiency plot for the double-aperture photometry (distribution of transit depth and durations)
- [x] Efficiency plot for the centroid shift (distribution of transit depths and durations)
### Week of the 21st of July 2025

Fix [[Major_Comments_Response#Second Major point| Major 2]] but in terms of:
- [x] Changing Section 7 of the paper (Results)  to include the variable EB transit parameters
	- [x] Write section 7.3: Scenario B: variable transit parameters
		- [x] Show the results about the value distribution of transit depths and duration
			- [x] Efficiency plot
			- [x] Efficiency table
	- [x] Write section 7.4: Comparison between scenarios
		- [x] Putting a table comparing the efficiency of both scenarios
- [ ] Changing Section 8 of the paper  (Discussion) to include the variable EB transit parameters

Fix  [[Major_Comments_Response#Fifth Major point|Major 5c]] and  [[Minor_Comments_Response#Tenth minor point|Minor 10]]:
- [x] Text explaining Eq. (39)
- [x] Caption of Fig. 7 related to Eq. (39)
- [x] Modify the ABSTRACT, more concretely the Results part of it. 

Fix [[Minor_Comments_Response#First minor point|Minor 1]]:
- [x] Delete the line that says that Kepler extracted photometry on-board (Kepler did not do that!)

Fix [[Major_Comments_Response#Ninth Major point| Major 9]]:
- [x] Refer to Borucki et al. for the definition of 'Super-Earth'
- [x] Recall that a planet with a significance higher than 7.1 could be, in  principle, detectable  for PLATO

### Week of the 28th of July 2025

Fix [[Major_Comments_Response#Eight Major point|Major 8]]:
- [x] Create the matrix/table with all the percentages of all the metrics
- [x] Reword the comparison: 'detected by X and not by Y' instead of 'detected only by X and not by Y'
- [ ] Write conclusions from the new matrix/table

Fix [[Minor_Comments_Response#Second minor point|Minor 2]]:
- [x] Explain the L3 data products -> not only the catalog from the P1 sample.

Fix [[Minor_Comments_Response#Third minor point|Minor 3]]:
- [x] Explain the required two quantities for COB and one for FX

Fix [[Minor_Comments_Response#Fifth minor point|Minor 5]]:
- [x] Cite my thesis, Fig. 7.1 to show the study of 2-pixel-ring extended mask efficiency
- [x] State that the current strategy is not the optimal way to build an extended mask

Fix [[Minor_Comments_Response#Eight minor point|Minor 8]]:
- [x] Explain the 10^-6  in the transit depth

Fix [[Minor_Comments_Response#Eleventh minor point|Minor 11]]:
- [x] Explain that the cutoff in magnitude comes from GAIA

Fix [[Minor_Comments_Response#Fourteenth minor point|Minor 14]]:
- [x] State the importance of having a mask library within budget

Fix [[Editorials_Response#First editorial|Editorial 1]]:
- [x] Change 'seismic techniques' and write 'asteroseismology'

Fix [[Editorials_Response#Second editorial|Editorial 2]]:
- [x] Cite the Nasa Exoplanet Archive for the number of exoplanets

Fix [[Editorials_Response#Third editorial|Editorial 3]]:
- [x] Mention clearly that PLATO will increase the detection number of Earth-like planets

Fix [[Editorials_Response#Fourth editorial|Editorial 4]]:
- [x] Write 'space-based and ground-based transit surveys'

Fix [[Editorials_Response#Fifth editorial|Editorial 5]]:
- [x] Assess the transit signal for validation

Fix [[Editorials_Response#Sixth editorial|Editorial 6]]:
- [x] Mention the two pointings for the P5 figure

Fix [[Editorials_Response#Seventh editorial|Editorial 7]]:
- [x] Cite Rauer et al. paper

Fix [[Editorials_Response#Eighth editorial|Editorial 8]]:
- [x] Transit depth in terms of shadowed stellar disc

Fix [[Editorials_Response#Ninenth editorial|Editorial 9]]:
- [x] Explain the 'weighted by the aperture' for centroid shifts

Fix [[Editorials_Response#Tenth editorial|Editorial 10]]:
- [x] Cite Samadi et al. (2019) for the microscanning

### Week of the 4th of August 2025

From [[Meeting_July_31_2025]]:
- [x] Remove any mention of the scenarios in the Abstract 
- [ ] Trim the Abstract. The Abstract is already too long. In total, it has more than 400 words. The word limit for the Abstract set by A&A is of 300 words
- [x] Remove scenario A from the results
- [x] Improve Appendix B
- [x] Diagnostic check: compute the difference between the values and the quadratic sum of the errors
- [x] Fix [[Minor_Comments_Response#Fourth minor point|Minor point 4]]
- [x] Remove discussion about Gamma, since we are only showing Scenario B in the paper
- [x] Check the code to be sure Gamma = 1 when using distribution of transit parameters
- [x] Set a deadline to send the paper to co-authors for revisions (14th of August)
- [ ] Re-do Section 8 now with the results from Scenario B

Fix [[Major_Comments_Response#Sixth Major point|Major point 6]]:
- [x] Explain the magnitude bins in the efficiency plots

Fix [[Major_Comments_Response#Eight Major point|Major Point 8]]:
- [x] Write the conclusions from the new added Matrix

Fix [[Minor_Comments_Response#Sixth minor point|Minor Point 6]]:
- [x] Write the explanation about the extended mask not being able to detect a true positive

Fix [[Minor_Comments_Response#Thirteenth minor point|Minor point 13]]:
- [x] Remove secondary masks from Figure B.2.
- [x] Mention in the Efficiency plots why we include the secondary masks even if the plot involves the target magnitude2

Fix  [[Major_Comments_Response#Third Major point|Major point 3]]:
- [ ] Contact Juan
- [ ] Investigate and read documentation for the PSF business
- [x] Section 2.3

### Week of the 11th of August

To-do:
- [x] Remove results from fixed transit parameters from Table 2
- [x] Remove any mention to fixed transit parameters result from the paper that is not the mention of my thesis (including gamma and <> factors)
- [x] Incorporating changes in Section 8
- [x] Trim the abstract
- [ ] Sending the paper to coauthors on Friday
- [x] Putting triangles as the contaminants in Fig. 1

Start to investigate what is going on at P=10.5 and P=11.5 for the:
- [x] Secondary flux measurements
- [x] Secondary centroid shifts
- [x] Extended flux measurements

### Week of 18th of August
  - [ ] P = 10.70 instead of P = 11 in the efficiency plots
  - [x] Table 3 for Secondary mask methods only
  - [ ]  Improve the answers from the referee using Réza's comments
  - [x]  Introduce the noise condition for the extended and nominal cobs as well in the code
  - [x]  Introduce the noise condition for the extended and nominal cobs in the paper efficiency equations 
  - [x]  Table 3 and new table 3 conclusion
  - [x]  Mention the 1 hour average explanation in Sect. 2.4 
  - [ ]  Don't mention real camera PSF in the referee response and in the paper
  - [x]  Check which grows faster, if the nominal significance or the extended significance or the secondary significance at P = 10.5 and P = 11.5
### Week of 25th of August
- [ ] Send an email to A&A about an extension
- [ ] Implement the new scaled over transit duration and number of transits efficiency condition for the COBs