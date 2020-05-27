# Process and analyze DATA collected for the paper : 
# Impact of plastic mulch film debris on soil physicochemical and hydrological properties
# PLEASE CHANGE THE WORKING DIRECTORY (WD) BEFORE RUNNING THE SCRIPT.
# For any issue or question, please contact Nicolas Beriot at nicolas.beriot@wur.nl

# Packages
  if(!require(readxl)){install.packages("readxl")}
  if(!require(ggplot2)){install.packages("ggplot2")}
  if(!require(multcompView)){install.packages("multcompView")}
  if(!require(tidyverse)){install.packages("tidyverse")}
  if(!require(factoextra)){install.packages("factoextra")} # Plot the PCA 
rm(list=ls()) # cleaning console
graphics.off() # cleaning plots

# color palette, 
  color_13=c("grey30","steelblue2", "orange2", "royalblue4", "red3" ,"steelblue2", "orange2", "royalblue4", "red3","steelblue2", "orange2", "royalblue4", "red3")

#//////////////////////////////////////////////////////////////////////////////////////////////////////   
# DATA loading ####
  wd="C:/Users/berio001/OneDrive - WageningenUR/PhD/Paper/Paper.0_Soil physics" # Write the directory with SAS_results.xlsx and Raw_data_Beriot.xlsx
  setwd(wd)
  
  # *** Raw data ####
    ring_Samples <-read_excel("Raw_Data_Beriot.xlsx", "Ring_R")
    ring_Samples$Set<-factor(ring_Samples$Set)
    cup_Samples <-read_excel("Raw_Data_Beriot.xlsx", "Cup_R")
  
    Parameters=c("n", "pb", "ks", "FC", "WDPT") # abreviation for the soil parameters
    
    Parameters_labels=c("Porosity [-]",                              # n
                       "Dry bulk density [kg/m³]",                  # pb
                       "Saturated hydraulic conductivity (log10)",  # ks
                       "Field Capacity [-]",                        # FC "Gravimetric water content at pF2 [-]",
                       "Water Drop Penetration Time [s]" )          # WDPT
    
    x_lab=c("Control", "LDPE-Mi", "Bio-Mi", "LDPE-Ma", "Bio-Ma","LDPE-Mi", "Bio-Mi", "LDPE-Ma", "Bio-Ma","LDPE-Mi", "Bio-Mi", "LDPE-Ma", "Bio-Ma")
    
    # Treatment label
      ring_Samples$Treat_conc[ring_Samples$Treat==1]<-"Control"
      ring_Samples$Treat_conc[ring_Samples$Treat==3 & ring_Samples$Set==1]<-"LDPE-Mi_0.5"
      ring_Samples$Treat_conc[ring_Samples$Treat==5 & ring_Samples$Set==1]<-"Bio-Mi_0.5"
      ring_Samples$Treat_conc[ring_Samples$Treat==7 & ring_Samples$Set==1]<-"LDPE-Ma_0.5"
      ring_Samples$Treat_conc[ring_Samples$Treat==8 & ring_Samples$Set==1]<-"Bio-Ma_0.5"
      ring_Samples$Treat_conc[ring_Samples$Treat==3 & ring_Samples$Set==2]<-"LDPE-Mi_1"
      ring_Samples$Treat_conc[ring_Samples$Treat==5 & ring_Samples$Set==2]<-"Bio-Mi_1"
      ring_Samples$Treat_conc[ring_Samples$Treat==7 & ring_Samples$Set==2]<-"LDPE-Ma_1"
      ring_Samples$Treat_conc[ring_Samples$Treat==8 & ring_Samples$Set==2]<-"Bio-Ma_1"
      ring_Samples$Treat_conc[ring_Samples$Treat==3 & ring_Samples$Set==3]<-"LDPE-Mi_2"
      ring_Samples$Treat_conc[ring_Samples$Treat==5 & ring_Samples$Set==3]<-"Bio-Mi_2"
      ring_Samples$Treat_conc[ring_Samples$Treat==7 & ring_Samples$Set==3]<-"LDPE-Ma_2"
      ring_Samples$Treat_conc[ring_Samples$Treat==8 & ring_Samples$Set==3]<-"Bio-Ma_2"
      ring_Samples$Treat_conc<-factor( ring_Samples$Treat_conc, 
                         levels = c("Control", "LDPE-Mi_0.5", "Bio-Mi_0.5", "LDPE-Ma_0.5", "Bio-Ma_0.5", 
                                    "LDPE-Mi_1", "Bio-Mi_1", "LDPE-Ma_1", "Bio-Ma_1", 
                                    "LDPE-Mi_2", "Bio-Mi_2", "LDPE-Ma_2", "Bio-Ma_2"))
      ring_Samples$ks=log10(ring_Samples$ks)

      
  # *** Plot Raw data ####
      # Choose a parameter to proceed: 
      # Parameters : 
      # "n"    "pb"   "ks"   "FC"   "WDPT"
      #  1      2       3     4       5
      p=1
      
      expand_y_letter=c(0.502,1.49,-1.95,0.2,11)
      pd <- position_dodge(0)
      
      # plot / treatment
      ggplot(ring_Samples, aes(x=Treat_conc, y=ring_Samples[,Parameters[p]][[1]], color=Treat_conc)) + 
        geom_point(position=pd, size=1)+   
        geom_jitter( width=0.11, size=1, alpha=0.9)+
        scale_color_manual(values=color_13)+
        expand_limits(y = c(NA, expand_y_letter[p]))+
        ylab(Parameters_labels[p])+
        theme(
          plot.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=16),
          axis.text.x  = element_text(angle=90, vjust=0.5, hjust=1, size=14), #hjust : align the labels 
          axis.text.y  = element_text(size=14),
          legend.position="none") #
      
      
      # plot / time replicat
      ggplot(ring_Samples, aes(x=Treat_conc, y=ring_Samples[,Parameters[p]][[1]], color=Set)) + 
        geom_point(position=pd, size=1)+   
        geom_jitter( width=0.11, size=1, alpha=0.9)+
        expand_limits(y = c(NA, expand_y_letter[p]))+
        ylab(Parameters_labels[p])+
        theme(
          plot.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=16),
          axis.text.x  = element_text(angle=90, vjust=0.5, hjust=1, size=14), #hjust : align the labels 
          axis.text.y  = element_text(size=14)) 
      
#//////////////////////////////////////////////////////////////////////////////////////////////////////   
# Estimated bulk density claculation ####    

  #plastic content: 
    ring_Samples$content[ring_Samples$Set==1]=0.005 # 0.5%
    ring_Samples$content[ring_Samples$Set==2]=0.01  # 1%
    ring_Samples$content[ring_Samples$Set==3]=0.02  # 2%
    ring_Samples$content[ring_Samples$Treat==1]=0   # control : 0%
  
  #plastic particle density pp:   
    ring_Samples$pp[ring_Samples$Treat==3 | ring_Samples$Treat==7] = 0.91
    ring_Samples$pp[ring_Samples$Treat==5 | ring_Samples$Treat==8] = 1.45
    ring_Samples$pp[ring_Samples$Treat==1 ] = 2.48
  
  #soil particle density ps:    
    ps=mean(  ring_Samples$pb[ring_Samples$Treat==1 ]/(1-ring_Samples$n[ring_Samples$Treat==1 ])  )  #ps~=2.48

  # Calculations of pb: 
    # Eq.S6: Relationship between dry bulk density and porosity considering soil particles only
    # Eq.S6: ??_b= (1-porosity) ??_S
      ring_Samples$pb_calc_0= (1-ring_Samples$n) *ps 
    
    # Eq.S7: Relationship between dry bulk density and porosity considering soil particles and plastic  
    # Eq.S7: ??_b= (1-porosty) * (??_P+1) *  (??_P * ??_S )/(??_P * ??_S  + ??_P )
      ring_Samples$pb_calc_1= (1 -ring_Samples$n) *
                              (1 +ring_Samples$content) * 
                              (ring_Samples$pp *ps) /
                              (ring_Samples$pp +ring_Samples$content *ps)
      
    # Linear regression equation: 
    # lm: ??_b= a1 + porosity *a2
      A=lm(pb~n, ring_Samples)
      ring_Samples$pb_calc_2= A[[1]][1] +ring_Samples$n * A[[1]][2]   

  
  # Calculation of the difference with the measured pb
    ring_Samples$pb_error_0=abs(ring_Samples$pb-ring_Samples$pb_calc_0)
    ring_Samples$pb_error_1=abs(ring_Samples$pb-ring_Samples$pb_calc_1)
    ring_Samples$pb_error_2=abs(ring_Samples$pb-ring_Samples$pb_calc_2)
    # mean difference 
      mean(ring_Samples$pb_error_0) 
      mean(ring_Samples$pb_error_1)
      mean(ring_Samples$pb_error_1)
      
  # Root-mean-square deviation measured/ observed
    RMSE = function(m, o){ 
      sqrt(mean((m - o)^2))
    }
 
    RMSE(ring_Samples$pb, ring_Samples$pb_calc_0)/mean(ring_Samples$pb)*100
    RMSE(ring_Samples$pb, ring_Samples$pb_calc_1)/mean(ring_Samples$pb)*100
    RMSE(ring_Samples$pb, ring_Samples$pb_calc_2)/mean(ring_Samples$pb)*100
    
  # Coefficient of determination
    Rsquare=function(y,z){
      m=mean(y)
      R=1-( sum((y-z)^2) / sum((y-m)^2) )
      return(R)
    }

    Rsquare(ring_Samples$pb,ring_Samples$pb_calc_0)
    Rsquare(ring_Samples$pb,ring_Samples$pb_calc_1)
    Rsquare(ring_Samples$pb,ring_Samples$pb_calc_2)
    
    Rsquare(ring_Samples$pb, (-1.77*ring_Samples$n+2.12))
    Rsquare(ring_Samples$n,(0.0454*ring_Samples$ks+0.576))
    Rsquare(ring_Samples$pb,(-0.145*ring_Samples$ks+0.895))
    
#//////////////////////////////////////////////////////////////////////////////////////////////////////   
# Linear regression analysis ####
  # Subset data
    P=subset(ring_Samples, select = c(n, pb, pb_calc_0, pb_calc_1, ks,Treat_conc ) )
    P$Treat_conc<-factor(P$Treat_conc, 
                         levels = c("Control", "LDPE-Mi_0.5", "Bio-Mi_0.5", "LDPE-Ma_0.5", "Bio-Ma_0.5", 
                                    "LDPE-Mi_1", "Bio-Mi_1", "LDPE-Ma_1", "Bio-Ma_1", 
                                    "LDPE-Mi_2", "Bio-Mi_2", "LDPE-Ma_2", "Bio-Ma_2"))

    
  # Linear regression y=ax+b // function to display the equation on graph
      lm_eqn <- function(df,x,y){
        m <- lm(y~x,df);
        eq <- substitute(italic(y) == a *" "*italic(x)*" + "* b *","~~italic(r)^2~"="~r2, 
                         list(a = format(unname(coef(m)[2]), digits = 3),
                              b = format(unname(coef(m)[1]), digits = 3),
                              r2 = format(summary(m)$r.squared, digits = 3)))
        as.character(as.expression(eq));
      }
 
  #*** Porosity and pb ####
    A=lm(pb~n, P) # Linear regression    
    
    # display the equation Eq.S6  on graph
     eq_S6 <- as.character(as.expression( substitute(
                        italic(y) == a *" "*italic(x)*" + "* b *","~~italic(r)^2~"="~r2, 
                          list(a = format(unname(ps), digits = 3),
                               b = format(unname(ps), digits = 3),
                              r2 = format( Rsquare(P$pb, (-ps*P$n + ps ) ), digits = 3) )
              )))
    
    # PLOT
      ggplot(P, aes(x=n, y=pb, color=Treat_conc, shape=Treat_conc, fill = Treat_conc)) + 
        geom_point(size=6, alpha=0.8,stroke=2)+                            # Alpha set the transparancy, stroke change the width of the shapes 
        scale_shape_manual(values=c(3,2,2,17,17,0,0,15,15,1,1,16,16))+    
        scale_color_manual(values=color_13)+
        scale_fill_manual(values=color_13)+
        geom_abline(intercept =  A[[1]][1], slope =  A[[1]][2], size=1  ) +
        geom_abline(intercept =  ps, slope =  -ps, size=1.5, color="green3", linetype="dashed"  ) + # Add a green dashed line for Eq.S6 
        geom_text(x =0.38, y =1.2,              
                  label = lm_eqn(P,x=P$n, y=P$pb),
                  parse = TRUE, colour="black", 
                  show.legend = FALSE, size=6 ) +  
        geom_text(x =0.38, y =1.25,
                  label = eq_S6,
                  parse = TRUE, colour="green3",
                  show.legend = FALSE, size=6 ) +
        xlab(Parameters_labels[1])+
        ylab(Parameters_labels[2])+
        theme(
          plot.title = element_blank(),
          axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          axis.text.x  = element_text(size=14), 
          axis.text.y  = element_text(size=14),
          legend.title = element_blank(),
          legend.text = element_text(size=14),
          legend.spacing.x = unit(0.5, unit = "cm"))
    # Export as png  width= 1500, height = 750
  
  
  #*** ks and Porosity #### 

    A=lm(n~(ks), P) # Linear regression  
    
    # PLOT
      ggplot(P, aes(x=(ks), y=n, color=Treat_conc, shape=Treat_conc, fill = Treat_conc)) + 
        geom_point(size=6, alpha=0.8,stroke=2)+                            # Alpha set the transparancy, stroke change the width of the shapes 
        scale_shape_manual(values=c(3,2,2,17,17,0,0,15,15,1,1,16,16))+    
        scale_color_manual(values=color_13)+
        scale_fill_manual(values=color_13)+
        geom_abline(intercept =  A[[1]][1], slope =  A[[1]][2], size=1  ) +
            geom_text(x =-2.5, y =0.37,               
                  label = lm_eqn(P,x=(P$ks), y=P$n),
                  parse = TRUE, colour="black", 
                  show.legend = FALSE, size=6 ) +   
        xlab(Parameters_labels[3])+
        ylab(Parameters_labels[1])+
        theme(
          plot.title = element_blank(),
          axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          axis.text.x  = element_text(size=14), 
          axis.text.y  = element_text(size=14),
          legend.title = element_blank(),
          legend.text = element_text(size=14),
          legend.spacing.x = unit(0.5, unit = "cm"))
    # Export as png  width= 1500, height = 750
  
  #*** ks and pb #### 
      A=lm(pb~(ks), P) # Linear regression 
      
      # PLOT
        ggplot(P, aes(x=(ks), y=pb, color=Treat_conc, shape=Treat_conc, fill = Treat_conc)) + 
          geom_point(size=6, alpha=0.8,stroke=2)+                            # Alpha set the transparancy, stroke change the width of the shapes 
          scale_shape_manual(values=c(3,2,2,17,17,0,0,15,15,1,1,16,16))+    
          scale_color_manual(values=color_13)+
          scale_fill_manual(values=color_13)+
          geom_abline(intercept =  A[[1]][1], slope =  A[[1]][2], size=1  ) +
          geom_text(x =-3.5, y =1.2,                 
                    label = lm_eqn(P,x=(P$ks), y=P$pb),
                    parse = TRUE, colour="black", 
                    show.legend = FALSE, size=6 ) +   
          xlab(Parameters_labels[3])+
          ylab(Parameters_labels[2])+
          theme(
            plot.title = element_blank(),
            axis.title.x = element_text(size=16),
            axis.title.y = element_text(size=16),
            axis.text.x  = element_text(size=14), 
            axis.text.y  = element_text(size=14),
            legend.title = element_blank(),
            legend.text = element_text(size=14),
            legend.spacing.x = unit(0.5, unit = "cm"))
      # Export as png  width= 1500, height = 750
    
#//////////////////////////////////////////////////////////////////////////////////////////////////////    
# Principal Component Analysis ####
  # Subset data        
    P=subset(ring_Samples, select = c(n, pb, ks, FC, WDPT ,Treat_conc ) )
    colnames(P)[1]="porosity"
    P$Treat_conc<-factor(P$Treat_conc, 
                           levels = c("Control", "LDPE-Mi_0.5", "Bio-Mi_0.5", "LDPE-Ma_0.5", "Bio-Ma_0.5", 
                                      "LDPE-Mi_1", "Bio-Mi_1", "LDPE-Ma_1", "Bio-Ma_1", 
                                      "LDPE-Mi_2", "Bio-Mi_2", "LDPE-Ma_2", "Bio-Ma_2"))
  # PCA
    ir.pca<-prcomp(subset(P, select=-Treat_conc),center = TRUE, scale = TRUE )
    summary(ir.pca)
    ir.pca[2]
  
  # PLOT
    fviz_pca_biplot(ir.pca, axes = c(1, 2), habillage=P$Treat_conc, palette = color_13, addEllipses = TRUE,     #change axes = to plot for other axes.
                    label = "var", col.var = "black", repel = TRUE, legend.title = "Treatment")+ #, title = tit) + #tit="titel PCA"
      scale_shape_manual(values=c(3,2,2,17,17,0,0,15,15,1,1,16,16))    
  # Export as png, width= 1000, height = 750  
  
#//////////////////////////////////////////////////////////////////////////////////////////////////////  
# Estimates and pair wise comparison #### 
    # Choose a parameter to proceed: 
    # Parameters : 
    # "n"    "pb"   "ks"   "FC"   "WDPT"
    #  1      2       3     4       5
    p=3 
    
    # *** Load the p-values ####
    
    Asterisk_control=data.frame( matrix("",13,length(Parameters)) )
    colnames(Asterisk_control)=Parameters
    levels(Asterisk_control$n)=c("", "*")
    p_val<-data.frame(read_excel("SAS_results_Beriot.xlsx", paste("p.",Parameters[p],sep = "") ))
    rownames(Asterisk_control)<-p_val$Name_txt
    
    for (pj in 1:length(Parameters)){
      levels(Asterisk_control[,Parameters[pj]])=c("", "*")
      p_val<-data.frame(read_excel("SAS_results_Beriot.xlsx", paste("p.",Parameters[p],sep = "") ))
      p_val[p_val=="<.0001"]=0.00001
      rownames(p_val)<-p_val$Name_txt
      PT<-as.matrix(subset(p_val, select=-Name_txt))
      PT<- as.matrix(apply(PT, 2, as.numeric))
      rownames(PT)<-p_val$Name_txt
      
      # asterisk
      Asterisk_control[PT[,"Control"]<0.001,Parameters[pj]]<-"*"
    }    
    
    rownames(Asterisk_control)<-p_val$Name_txt
    
    # *** Load estimates ####
      estimates <-read_excel("SAS_results_Beriot.xlsx",Parameters[p])
      #re-ordering, adding letters and re-labeling the factors
        Order_letter=as.vector(as.matrix.data.frame(estimates[,Parameters[p]]))
        estimates$Treat<-factor(estimates$Treat, levels = estimates$Treat[order(Order_letter)] )
        estimates<- estimates[order(Order_letter, decreasing=FALSE), ]
        estimates$letter=Asterisk_control[,Parameters[p]]
      
      levels( estimates$Treat)[levels( estimates$Treat)=="Ck"]="Control"
      levels( estimates$Treat)[levels( estimates$Treat)== "µPE.1"]="LDPE-Mi_0.5"
      levels( estimates$Treat)[levels( estimates$Treat)== "µBio.1"]="Bio-Mi_0.5"
      levels( estimates$Treat)[levels( estimates$Treat)== "MPE.1"]="LDPE-Ma_0.5"
      levels( estimates$Treat)[levels( estimates$Treat)== "MBio.1"]="Bio-Ma_0.5"
      levels( estimates$Treat)[levels( estimates$Treat)== "µPE.2"]="LDPE-Mi_1"
      levels( estimates$Treat)[levels( estimates$Treat)== "µBio.2"]="Bio-Mi_1"
      levels( estimates$Treat)[levels( estimates$Treat)== "MPE.2"]="LDPE-Ma_1"
      levels( estimates$Treat)[levels( estimates$Treat)== "MBio.2"]="Bio-Ma_1"
      levels( estimates$Treat)[levels( estimates$Treat)== "µPE.3"]="LDPE-Mi_2"
      levels( estimates$Treat)[levels( estimates$Treat)== "µBio.3"]="Bio-Mi_2"
      levels( estimates$Treat)[levels( estimates$Treat)== "MPE.3"]="LDPE-Ma_2"
      levels( estimates$Treat)[levels( estimates$Treat)== "Mbio.3"]="Bio-Ma_2"
    
    
    # *** Plot ####
      V_letter=c(-2.8, -1.8, -2.4,-1.1, -1.2 ) #verticale offset for the asterix
      expand_y_letter=c(0.502,1.49,-1.95,0.2,11)
      pd <- position_dodge(0)
    
      ggplot(estimates, aes(x=estimates$Treat, y=estimates$Estimate)) + 
        geom_errorbar((aes(ymin=estimates$Estimate-estimates$StE, 
                           ymax=estimates$Estimate+estimates$StE)),
                      colour=color_13, width=.3, size=2, position=pd) + # create the bar, size=4
        geom_point(position=pd, colour=color_13, size=4)+   # create a point in the midel of the bar, size=4
        geom_text(aes(label = estimates$letter),fontface="bold", position = position_dodge(0.9), size=7, vjust =V_letter[p] )+  # Add the letters, vjust: vertical offset   -4.5, -3.5, -4, -2 V_letter[p]=c(-4, 2.8, 3.3,-1.5 ) #, hjust=-0.1
        expand_limits(y = c(NA, expand_y_letter[p]))+
        ylab(Parameters_labels[p])+
        #ylab(expression('Saturated hydraulic conductivity (log  '[10]*')'))+    #special for ks
        #expand_limits(y = c(NA,max(estimates$Estimate)))+
        theme(
          plot.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=16),
          axis.text.x  = element_text(angle=90, vjust=0.5, hjust=1, size=14), #hjust : align the labels 
          axis.text.y  = element_text(size=14) )
      #Export as png,  width= 500, height = 500
    