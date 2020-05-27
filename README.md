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
	  - calculate the estimated bulk density with the measure of porosity and the plastic content<br />
	  - perform the linear regressions between the measured parameters<br />
	  - perform the principal component analysis<br />
	  - plot the estimates for each each parameter and each combination of the factors Content*Type*Size<br />

# Abstract
Plastic residues in soil have attracted growing attention in recent years. The use of plastic mulch films used in agriculture has beenare strongly considered as to be a major source of the plasticthis residues found in soil. Mulching with low-density polyethylene (LDPE) is widely practiced and the resulting macro- and microscopic plastic residues in agricultural soil have aroused concerns for years. Over theIn the lastpast decades, a variety of biodegradable (Bio) plastics have been introduced developed in the hopes of reducing the plastic contamination ofin the terrestrial ecosystem. However, the impact of these biodegradable Bio plastics in agroecosystems have not been sufficiently studied. Therefore, we investigated the impact of macro (Ma, around 5 mm) and micro (Mi, < 1 mm) sized plastic mulch film debris from LDPE and one commercially available type of starch-based biodegradable (Bio) mulch film on soil physicochemical and hydrological properties. We used environmentally relevant concentrations of plastics, ranging from 0 to 2% (w/w), identified by field studies and literature review. We studied the effects of the plastic residue on a sandy soil for one month in a laboratory experiment. The bulk density, porosity, saturated hydraulic conductivity, field capacity and soil water repellency were altered significantly in the presence of the four kinds of plastic debris, while pH, electrical conductivity and aggregate stability were not substantially affected in this study. Overall, our research provides clear experimental evidence that microplastics affect soil properties. The type, size and, content of plastic debris and theiras well as the interactions between these three factors played complex roles in the variations of the measured soil parameters. Living in a plastic eraage, it is crucial to conduct further interdisciplinary studies in order to have a comprehensive understanding of plastic debris in soil and agroecosystems.

# Abbreviations
	n, Porosity
	pb, Dry bulk density 
	ks, Saturated hydraulic conductivity
	FC, Field capacity (gravimetric water content at pF2) 
	WDPT, Water drop penetration time
	pH, Potential for hydrogen
	EC, Electrical conductivity 
	ASI, Aggregates stability index
