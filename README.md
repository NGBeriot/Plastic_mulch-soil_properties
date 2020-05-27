# Plastic_mulch-soil_properties
  Raw data and script for generating the graphs published in : citation

# Content description
  Raw_Data_Beriot.xlsx<br />
	  contains the measurement for all the parameters. Ring_R contains the measurements done on the ring samples, namely ks, n, pb, FC, WDPT. 
	  Cup_R contains all the measurements done on the air dried loose samples, namely pH, EC, ASI.
  
  SAS_results_Beriot.xlsx<br />
    contains the outcomes of the linear mixed effect model implemented in SAS® 9.4.
	  The eight first sheets contain the estimates,	standard Error and associated p-value for the combination of the factors Content*Type*Size (n, pb, ks, FC, WDPT, pH, ASI).
	  The eight next sheets contain the p-value obtained after the pairwise comparison (p.n, p.pb, p.ks, p.FC, p.WDPT, p.pH, p.ASI)

  Plastic_mulch-soil_properties_Analysis.R <br />
	contains the script to perform to: <br />
	  - calculate the estimated bulk density with the measure of porosity and the plastic content<br />
	  - perform the linear regressions between the measured parameters<br />
	  - perform the principal component analysis<br />
	  - plot the estimates for each each parameter and each combination of the factors Content*Type*Size<br />

# Abbreviations
	n, Porosity
	pb, Dry bulk density 
	ks, Saturated hydraulic conductivity
	FC, Field capacity (gravimetric water content at pF2) 
	WDPT, Water drop penetration time
	pH, Potential for hydrogen
	EC, Electrical conductivity 
	ASI, Aggregates stability index
