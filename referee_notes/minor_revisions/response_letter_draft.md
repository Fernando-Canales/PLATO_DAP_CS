# Response Letter Draft 

Dear referee,  

Thank you very much for the constructive revision of the manuscript. We greatly appreciated this second round of revisions.  We have carefully considered all your comments and have revised the manuscript accordingly. This letter contains our detailed point-by-point responses.  

We hope that our responses make the paper suitable now for publication in A&A.  

Best regards,  
Dr. Fernando Gutiérrez Canales  

## Minor Points

**Comment:** I also concede that assigning an EB to all background stars, while unrealistic, gives a conservative upper bound to the results. However, I don’t see this mentioned explicitly in the revised paper. **Please explicitly state that the assumption that all contaminants are EBs means that your results are a conservative upper limit.** The obvious place to do so is when you state this assumption, in the final paragraph of section 5.2.  

**Response:** This is a valid point. However, it is a bit subtle. Our aim is not to set a realistic estimation for the expected false positives for PLATO. Our aim is to compare different metrics in order to assess which is the most effective at detecting false positives. To do this, we need to consider all situations in which we can potentially have an EB-generated transit. As consequence, we consider that every contaminant is an EB, and then evaluate for each metric its ability to detect a FP. We do, however, concede that this can be seen as a worst-case scenario approach. Having said that, at the end of section 5.2 we have added the following text: "We note that assuming all contaminants are EBs represents a worst-case scenario for contamination. This allow us to assess the maximum capability of each detection method under challenging conditions."

> RS:  I don't agree with the referee's point of view. How can our estimates be considered an upper limit? 
Once again, our aim is to compare different metrics in order to assess which is the most effective at detecting false positives.
To do this, we need to consider all situations in which we can potentially have an EB-generated transit.  Consequently, we consider 
 that every contaminant is an EB, and then evaluate for each metric its ability to detect a false positive.     The effectiveness (or efficiency) of a metric in detecting false positives is then the ratio between the number of cases where a false positive is detectable and the number of transits detectable in the nominal flux.  <


**Comment:** Referencing the lead authors’ thesis for the detailed demonstration that the metrics’ efficiency does not significantly change for different PSFs is a good step. However, **is there a reference that could be added for the Monte Carlo simulations that generated the camera-to-camera variations, or were those simulations performed as part of this work?**  

**Response:**  This is a valid point. These Monte Carlos simulations were carried out by a member of the systems team using an optical model of the camera (built with the optical software Zemax). The corresponding PSF data set was distributed to the PLATO consortium. A technical note describing this data set exists but since it is an internal PLATO document is not suitable for citation. We however added the following text to Section 8.4 in the paper: "Such PSF data sets were provided by Dr. Thierry Appourchaux from IAS."

>RS: These Monte Carlo simulations were carried out by a member of the systems team, using an optical model of the camera (built with Zemax software).
The corresponding PSF data set was distributed to the PLATO consortium.  A technical note describing this data set has been included in the list of references. <

**Comment:** Thank you for these extensive changes. However, I'm still not entirely convinced by Figure 5a. While the extended mask centroid is an improvement over the nominal centroid case, it is not clearly so. For the brightest stars (P=10), the results are consistent between these two methods, as they are at P=11. But at other magnitudes there is a greater difference.  

In particular, the results for P=12.5 catch my eye. Here the extended mask centroid efficiency increases compared to slightly brighter (P=12) stars. Why is that? Is it a numerical artefact?  

>RS: I think you can show that the increase of the efficiency in the extended centroid between 12 and 12.5 remains within 3 sigma, so it is not statically significant.<


Furthermore, the relative efficiency of the various extended mask centroids and nominal mask centroids varies considerably, with the ordering of most to least efficiency changing from magnitude bin to magnitude bin. **Could you comment on why this is the case?** This variation means that the statement on page 13 (left-hand column), the paragraph referring to Figure 5a, that “the extended centroids are the most efficient for either 24 or 6 cameras” is an incorrect generalisation; counterpoints can be seen in the figure. **Please justify this general statement or explain it further in the text.**  

>RS: I think the statement that “the extended centroids are the most efficient for either 24 or 6 cameras” is not    relevant when referring to the figure 5. since indeed the   relative efficiency of the various centroid metric varies with the magnitude;
The fact that overall the extended centroid is the more efficient than the nominal centroid is related to the averaged  efficiency reported in table 2. I then propose to remove this sentence. 
<
 
Regarding the manuscript changes listed in the response, it is not clear to me which changes the third bullet point is referring to, as I can find no discussion of magnitude-dependent camera sensitivity.  

**Response:**  We thank the referee for this comment. However, the increase of efficiency in the extended centroid between 12 and 12.5 remains within 3 sigma. Therefore, it is not a significant change. Furthermore, we removed the sentence "the extended centroids are the most efficient for either 24 or 6 cameras" since it is in fact not relevant when referring to Fig. 5 since indeed the relative efficiency of the various centroid metric varies with the magnitude.  The fact that, overall, the extended centroid is the more efficient than the nominal centroid is related to the averaged efficiency reported in Table 2. Furthermore, we found an additional typo next to the sentence mentioned by the referee. We found the sentence that reads "The nominal centroids are very close in efficiency to the nominal ones." when it should say "The nominal centroids are very close in efficiency to the extended ones." This has been corrected as well in the paper.   

Furthermore, the third bullet point in our response to Major Point 5 refers to the strange peaks in efficiency at P=10.5 and P=11 that our first results used to have. Such strange efficiency peaks vanished when implementing the new results (EBs sampled transit duration and depth and the other changes detailed in the response to the major revisions). We concede that the third bullet point of the section "Manuscript changes" of our Response to Major point 5 was misleading as the referee spot it and we appreciate it.   


**Comment:** However, the new text discussing the results in Table 2 could be improved. It states that “_Table 2 shows that nominal and extended centroids could be prioritized where computationally feasible given their high efficiency and their capability to detect FP in several contaminants_”. However, why would you prioritise these methods when secondary flux (SFX) is both significantly more efficient and computationally cheaper, as is argued in the next sentence of the text! **Please could you clarify this suggestion.**  

**Response:**  This is a valid point. We are glad the referee spot it. The sentence "Table 2 shows that nominal and extended centroids could be prioritized where computationally feasible given their high efficiency and their capability to detect FP in several contaminants" has been removed of the paragraph. Instead, we have added the following text after Table 2: "Secondary flux is particularly effective when only one contaminant is capable of generating a FP, as secondary masks are specifica lly designed for single-contaminant scenarios. However, when multiple contaminants can generate FPs, nominal and extended centroids become preferable due to their ability to simultaneously monitor several contaminants and their higher efficiency compared to extended flux, though only where computational resources allow it". And then we refer to Sect. 8.2 where we discuss these considerations when presenting the overall strategy for selecting the best metrics.

>RS:   SFX is the most effective at detecting FP when there is only one contaminant capable of generating a FP.
However, the secondary mask can only be defined for a single contaminant. In cases where several contaminants are likely to generate a PF, another metric should be preferred. In this case, the nominal and extended centroids should be preferred to the extended flow because of their greater efficiency, but only if resources allow. This point is discussed in details in section 8.2<

**Comment:** However, with the deletion of that paragraph, there is now no explanation of the penultimate row of Table 2, which was included to allow a more direct comparison between the secondary flux and extended flux cases. **Is the penultimate row of Table 2 still needed?** Similarly, I question the use of the final row of Table 2 given the updated results and discussion. **Is the final row of Table 2 still useful?**  

>RS: to be check but it seams indeed that the penultimate row is not commented. Then it should be removed. <

**Response:**  This is a valid point. The penultimate row of Table 2 has been removed. However, the last row of Table 2 it is still useful since constitutes an example for  the discussion at the beginning of Section 8.1 about the differential analysis. 

**Comment:** However, the final listed manuscript change, the removal of the footnote, as not been implemented. The footnote is still present in the revised paper but is now numbered as footnote 8. **Please remove this footnote, as you stated you had done.**  

**Response:** This is a valid point. The footnote is now removed. We apologize for not having done it in the previous round of revisions.  


**Comment:** The use of the Borucki et al. (1996) definitions to derive your Super-Earth values is a good one, and I appreciate the inclusion of a reference to Marchiori et al. (2019) for the depth and duration values for the Jovian and Earth-like planets. However, the inconsistency between the figure and the caption still exists! The Caption has not been corrected – it still states that “_the green circles refer to Earth-like planets using 24 cameras as well_”, while the figure legend refers to 6 cameras (consistent with your response). **Please make sure that the caption is correct.**   

**Response:** This is a valid point. The caption now says: "the green circles refer to Earth-like planets using 6 cameras as well". 

**Comment:** The updated caption for figure 6 does not include the clarification that it is showing the averaged uncertainty. **Please make sure that this change is implemented**.  

**Response:** This is a valid point. The caption for figure 6 now says: "Distribution of  the centroid shift averaged uncertainty (Eq.(37)) along target P magnitude for the extended (blue triangles) and nominal (black circles) masks."

## Editorials

Here we list the editorials spotted from the referee and our responses.  

1.- In the reference list, PLATO Science Requierement Document (2021)  -> PLATO Science **Requirement** Document (2021)  

**Response:** We took care of this typo, now the reference appear as: "PLATO Science Requirements Document (2021)"

2.- What is the difference between the following documents in the reference list?
    
    - ESA, 2021, PLATO Science Requirement Document (2021)
        
    - PLATO SciRD, 2019, PLATO Science Requirements Document, PLATO SciRD (2019)  
  
**Response:** This is a valid point. We are glad the reference spotted this inconsistency. The PLATO Science  Requirements Document is a valuable source for our work. However, due to a legacy mistake both the 2019 and 2021 versions appear in the paper. The correct thing is to show only the 2021 in the references.  Now the only reference in our paper to the PLATO Science Requirements Document is to the one from 2021. The one from 2019 does not longer appear in the paper.  

3.- Is it "on board" or "on-board"? Please be consistent.  

**Response:** This is also a valid point. We have decided for the easiest solution and implemented the form "on-board" to every case.  

4.- Page 8, right-hand column, section 5.2, penultimate paragraph. The text reads "at low Galactic latitudes for TESS for PLATO". Should this just be "for TESS" or just "for PLATO"?  

**Response:** This is also a valid point. The correct sentence should say "at low Galactic latitudes for TESS". It has been implemented in the paper accordingly.  

5.- Please rasterise Figure 6, if possible. There are clearly a very large number of points plotted in this figure, making it slow to load and affecting the performance of the PDF.  

**Response:** This is a valid point. Figure 6 has been rasterized as requested to improve the PDF performance.  

6.- In Table 2, Seconday flux -> **Secondary** flux  

**Response:** This is a valid point. It has been corrected accordingly.  

7.- In Table 2, Extented flux -> **Extended** flux  

**Response:** This is a valid point. It has been corrected accordingly.  

8.- Page 11, right-hand column, Section 6.2.2, last paragraph: "For the same reason than the nominal..." -> "For the same reason **as for** the nominal..."  

**Response:** This is a valid point. It has been corrected accordingly.


