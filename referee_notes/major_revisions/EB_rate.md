
# Assumptions about EBs

### Papers on EB occurrence rates in space missions

I found several papers about EB rates for space missions. So far:
- [Prsa et al. 2011][https://iopscience.iop.org/article/10.1088/0004-6256/141/3/83]  paper on Kepler. It provides an EB occurrence rate of 1.2%. Furthermore, it mentions that:

> Hipparcos found 917 Eclipsing Binaries in a sample of 118, 218 stars, which correspond to 0.82% occurrence rate. Based on our conservative estimates, we see a ~1.2% occurrence rate. This ~50% increase in detection attests for Kepler photometric superiority. Even this rate might be slightly underestimated because of the dominant selection effect: short timescale.


Fig. 12 of the paper shows an EB occurrence rate of ~1.2% for Galactic latitudes above 9 degrees.
-![[images/Pasted image 20250704152748.png]]


- [Prsa et al 2022] [https://iopscience.iop.org/article/10.3847/1538-4365/ac324a] paper on TESS. It provides an EB occurrence rate  in terms of several parameters. For instance, their Fig. 15 shows  a distribution of Detected EB rate[%] as a function of Galactic latitude. The Fig. shows an EB rate of around 3% at a Galactic latitude of 0°.  However,  PLATO will observe at higher Galactic latitudes:  [Bray et al. 2023][https://academic.oup.com/mnras/article/518/3/3637/6825462] stated that:

	> the Northern hemisphere LOP (LOPN), centered on Galactic coordinates l = 65° and b = 30° and the Southern hemisphere LOP (LOPS), centered on l = 253° and b = −30°

 
![[images/Pasted image 20250704150554.png]]

If we however consider Fig. 14 we can take the ratio of EBs (orange bars) over all 2-min targets (blue bars) and obtain a ratio of ~1% which agrees with the Prsa paper on Kepler mentioned above. 
![[images/Pasted image 20250704152326.png]]

- [Mowali et al. 2023][https://ui.adsabs.harvard.edu/abs/2023A&A...674A..16M] paper on GAIA. This is a catalog paper more than an EB rate paper, just like [Söderhjelm and Dischler 2005][http://www.aanda.org/10.1051/0004-6361:20042541] . But it is useful as well. It says that they found 2 million EB candidates in GAIA data. If we consider that GAIA has observed around 2 billion stars. Then we can make a very rough estimate of the EB occurrence rate of 0.1%.  This paper is currently mentioned in the manuscript.

### Changes in the code already implemented

I have changed dap_metrics.py in the following way to implement the EB occurrence rate in the calculations:
- More or less at line 
![[images/Pasted image 20250709104430.png]]
![[images/Pasted image 20250704171240.png]]
  ![[images/Pasted image 20250704171257.png]]



To Do:
- [ ] Maybe something about EB occurrence rate on CoRoT. Shouldn't be that different.
- [ ] Maybe try a 1.5% EB rate on my code or another number different than 1%. Might not change things a lot. 
- [ ] Add a reference to my PhD thesis