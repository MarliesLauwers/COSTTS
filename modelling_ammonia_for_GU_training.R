rm(list=ls())
library(ggplot2)
library(deSolve)
library(kimisc)
library(dplyr)
library(stringr)
#library(mosaicCalc)
#library(mosaic)
library(chron)
setwd("/Users/kreozotica/Desktop/Manure_NH3/") # open file with ambient parameters 
ambient_param = read.table(paste0("manure_val.txt"),header=T)
source("concentrat.R")
t_interval = c(0, 86000);
#init_cond = c(120,  0);
Sm=2.83   #[Mol/m3*sec] max conversin rate at high urea concentration
Km=2000   #[Mol/m3] Michaelis const
#ti = (0:23)*3600;
A=1       #[m^2] emitted surface area
V=0.01    # [m^3] volume of layer
##########################################################################
##########################################################################
Kd=0.81*10^(-10)  #adapted acid dissociation constant for NH3
##########################################################################
##########################################################################
##########################################################################
################## playing around with pH values #########################
for ( i in 1:12){ 
  v=ambient_param$Vel[i] #[m/sec] air velocity 
  T=ambient_param$T_air[i]+273
  ti=ambient_param$T_air[i]
  UC=ambient_param$UC[i];
  
  m=1.07^(T-293);
  k = 48.439*(v^0.8)*(T^-1.4) # mass transfer coefficent [m/s]
  Henry_k = 1384*(1.053^(293-T)) #g/cm3*atm Henry's const  
  
  emi_diff_pH<-NULL
  ############ trying different pH values #############
  for ( ii in seq(6.8, 8.8, by=0.4)){  
    pH=ii
    #print(pH)
    f = 1/(1+(10^(-pH)/(Kd*m)));
    ############### INTEGRATION ##################
    initial_cond <- c(U = UC, C = 0)
    times <- seq(from = 0, to = 86400, by = 60)
    
    out   <- ode(y = initial_cond, times = times, func = concentrat, parms = NULL)
    
    HMS = seconds.to.hms(out[,1] )  
    
    integration_output<- data.frame(HMS,out[,2],out[,3])
    colnames(integration_output)<-c("time_only", "U", "Emi_cntr")
    
    ############### calculating Emi out of Integration Param ##################
    E=(k*A*integration_output$Emi_cntr*f)/Henry_k # amammonia emisions
    emi_ammonia_mg_sec=E*17.03*1000*60; #17.03[g/mol] - molar mass of Ammonia|
    integration_output<- cbind(integration_output,emi_ammonia_mg_sec)
    
    emi_diff_pH<-rbind(emi_diff_pH,data.frame(ii,integration_output))# combining all integration outputs into one data frame (with diff.pH)
    
  }
  colnames(emi_diff_pH)[1]<-c("pH") # rename the first column in the data frame
  time=as.POSIXct(emi_diff_pH$time_only,format="%H:%M:%S",tz="Europe/Berlin") # transfering the time data into POSIX
  emi_diff_pH<-cbind(emi_diff_pH,time) # create a data frame with the new column
  emi_diff_pH$pH_double = emi_diff_pH$pH #doublicating PH column
  
  emi_diff_pH$pH_double <- paste(emi_diff_pH$pH_double, c("pH"), sep = " ") # add the word PH for nice ggploting 

  ###################################################################
  ########### ADD calculations for measured pH ##############
  ###################################################################
  
  source("concentrat_msrd.R")
  pH_msrd=ambient_param$pH[i]
  T_msrd=ambient_param$T_air[i]+273
  f_msrd = 1/(1+(10^(-pH_msrd)/(Kd*m)));
  k_msrd = 48.439*(v^0.8)*(T_msrd^-1.4) # mass transfer coefficent [m/s]
  Henry_k_msrd = 1384*(1.053^(293-T_msrd)) #g/cm3*atm Henry's const  
  ############### INTEGRATION ##################
  initial_cond <- c(U = UC, C = 0)
  times <- seq(from = 0, to = 86400, by = 60)
  out_msrd   <- ode(y = initial_cond, times = times, func = concentrat_msrd, parms = NULL)
  HMS_msrd = seconds.to.hms(out_msrd[,1] )  
  integration_output_msrd<- data.frame(HMS_msrd,out_msrd[,2],out_msrd[,3])
  colnames(integration_output_msrd)<-c("time_only", "U", "Emi_cntr")
  ############### calculating Emi out of Integration Param ##################
  E_msrd=(k_msrd*A*integration_output_msrd$Emi_cntr*f_msrd)/Henry_k_msrd # amammonia emisions
  emi_ammonia_mg_sec_msrd=E_msrd*17.03*1000*60; #17.03[g/mol] - molar mass of Ammonia|
  integration_output_msrd<- cbind(pH_msrd,integration_output_msrd,emi_ammonia_mg_sec_msrd)
  
  colnames(integration_output_msrd)[1]<-c("pH") # rename the first column in the data frame
  integration_output_msrd<-cbind(integration_output_msrd,time) # create a data frame with the new column
  integration_output_msrd$pH_double = integration_output_msrd$pH #doublicating PH column
  
  integration_output_msrd$pH_double <- paste(integration_output_msrd$pH_double, c("pH_msrd"), sep = " ") # add the word PH for nice ggploting
  
  
    #################### PLOTTING ####################
  plot<-ggplot(emi_diff_pH) + geom_line(aes(x=time, y=emi_ammonia_mg_sec, color=pH_double), size = 1.2)  +
    
    geom_line(aes(x=integration_output_msrd$time,y=integration_output_msrd$emi_ammonia_mg_sec,color=integration_output_msrd$pH_double),linetype="dotted",  size = 1.2) + 
    
    
    
    scale_x_datetime(date_breaks = "4 hour", date_labels = "%I:%M %p") + 
    labs( x = "time [hours]", y = "Emissions [mg/min]") + ylim(0, 50) +
    theme(text=element_text(size=14))  #family="Times"))
  
  
  plot<-plot+ggtitle(paste0("Model 2: AE, different pH levels for SL #",i ))
  
  
  
  print(plot)
  
  #### working way how to save in a good resoultion with converting 
  ggsave(plot, filename = paste("dff_pH_SL",i,".tiff", sep=""), dpi=600,   # ES&T. 300-600 at PLOS One,
         compression = "lzw",  # PLOS One 44MB w/out compression, 221KB w/compression
         family="Arial",  # ggplot default.  Could set to others depending on journal
         type="cairo")  # PLOS One 44MB w/out compression, 221KB w/compression)
  dev.off() 
  
}

dev.off() 


####################################################################################
####################################################################################
####################################################################################
################## playing around with  Temperature values #########################
####################################################################################
####################################################################################
####################################################################################

rm(list=ls())
library(ggplot2)
library(deSolve)
library(kimisc)
library(dplyr)
library(stringr)
library(ggrepel)
library(patchwork)
#library(mosaicCalc)
#library(mosaic)
library(chron)
setwd("/Users/kreozotica/Desktop/Manure_NH3/") # open file with ambient parameters 
ambient_param = read.table(paste0("manure_val.txt"),header=T)
source("concentrat.R")
 
t_interval = c(0, 86000);
#init_cond = c(120,  0);
Sm=2.83   #[Mol/m3*sec] max conversin rate at high urea concentration
Km=2000   #[Mol/m3] Michaelis const
#ti = (0:23)*3600;
A=1       #[m^2] emitted surface area
V=0.01    # [m^3] volume of layer
##########################################################################
##########################################################################
Kd=0.81*10^(-10)  #adapted acid dissociation constant for NH3
##########################################################################
##########################################################################
##########################################################################
################## playing around with pH values #########################
for ( i in 1:12){ 
  v=ambient_param$Vel[i] #[m/sec] air velocity 
  #T=ambient_param$T_air[i]+273
  #ti=ambient_param$T_air[i]
  UC=ambient_param$UC[i];
  pH=ambient_param$pH[i];
  emi_diff_T<-NULL
  ############ trying different pH values #############
  for ( ii in seq(6.0, 12.0, by=1)){  
    T=ii + 273 # need to present the temperature in Kelvin
    m=1.07^(T-293);
    f = 1/(1+(10^(-pH)/(Kd*m)));
    k = 48.439*(v^0.8)*(T^-1.4) # mass transfer coefficent [m/s]
    Henry_k = 1384*(1.053^(293-T)) #g/cm3*atm Henry's const  
    ############### INTEGRATION ##################
    initial_cond <- c(U = UC, C = 0)
    times <- seq(from = 0, to = 86400, by = 60)
    
    out   <- ode(y = initial_cond, times = times, func = concentrat, parms = NULL)
    
    HMS = seconds.to.hms(out[,1] )  
    
    integration_output<- data.frame(HMS,out[,2],out[,3])
    colnames(integration_output)<-c("time_only", "U", "Emi_cntr")
    

    ############### calculating Emi out of Integration Param ##################
    E=(k*A*integration_output$Emi_cntr*f)/Henry_k # amammonia emisions
    emi_ammonia_mg_sec=E*17.03*1000*60; #17.03[g/mol] - molar mass of Ammonia|
    integration_output<- cbind(integration_output,emi_ammonia_mg_sec)
    
    emi_diff_T<-rbind(emi_diff_T,data.frame(ii,integration_output))# combining all integration outputs into one data frame (with diff.pH)
    
  }
  
  ###################################################################
  ################### preparation for plotting  #####################
  ###################################################################
  
  colnames(emi_diff_T)[1]<-c("T") # rename the first column in the data frame
  time=as.POSIXct(emi_diff_T$time_only,format="%H:%M:%S",tz="Europe/Berlin") # transfering the time data into POSIX
  emi_diff_T<-cbind(emi_diff_T,time) # create a data frame with the new column
  emi_diff_T$T_double = emi_diff_T$T #doublicating PH column
  
  emi_diff_T$T_double <- paste(emi_diff_T$T_double, c("T"), sep = " ") # add the word PH for nice ggploting 
   
  ###################################################################
  ########### ADD calculations for measured temerature ##############
  ###################################################################
  source("concentrat_msrd.R")
  T_msrd=ambient_param$T_air[i]+273
  T_msrd_Cel<-ambient_param$T_air[i]
  m_msrd=1.07^(T-293);
  f_msrd = 1/(1+(10^(-pH)/(Kd*m_msrd)));
  k_msrd = 48.439*(v^0.8)*(T_msrd^-1.4) # mass transfer coefficent [m/s]
  Henry_k_msrd = 1384*(1.053^(293-T_msrd)) #g/cm3*atm Henry's const  
  ############### INTEGRATION ##################
  initial_cond <- c(U = UC, C = 0)
  times <- seq(from = 0, to = 86400, by = 60)
  out_msrd   <- ode(y = initial_cond, times = times, func = concentrat_msrd, parms = NULL)
  HMS_msrd = seconds.to.hms(out_msrd[,1] )  
  integration_output_msrd<- data.frame(HMS_msrd,out_msrd[,2],out_msrd[,3])
  colnames(integration_output_msrd)<-c("time_only", "U", "Emi_cntr")
  ############### calculating Emi out of Integration Param ##################
  E_msrd=(k_msrd*A*integration_output_msrd$Emi_cntr*f_msrd)/Henry_k_msrd # amammonia emisions
  emi_ammonia_mg_sec_msrd=E_msrd*17.03*1000*60; #17.03[g/mol] - molar mass of Ammonia|
  integration_output_msrd<- cbind(T_msrd_Cel,integration_output_msrd,emi_ammonia_mg_sec_msrd)
  
  colnames(integration_output_msrd)[1]<-c("T") # rename the first column in the data frame
  integration_output_msrd<-cbind(integration_output_msrd,time) # create a data frame with the new column
  integration_output_msrd$T_double = integration_output_msrd$T #doublicating PH column
  
  integration_output_msrd$T_double <- paste(integration_output_msrd$T_double, c("Tm"), sep = " ") # add the word PH for nice ggploting 
  ############################################################
  ######################## PLOTTING ##########################
  ############################################################
  plot<-ggplot(emi_diff_T) + geom_line(aes(x=time, y=emi_ammonia_mg_sec, color=T_double), size = 1.2) +
    
   
   
     geom_line(aes(x=integration_output_msrd$time,y=integration_output_msrd$emi_ammonia_mg_sec,color=integration_output_msrd$T_double),linetype="dotted", size = 1.2) + 
    
     
    
    scale_x_datetime(date_breaks = "4 hour", date_labels = "%I:%M %p") + 
    labs( x = "time [hours]", y = "Emissions [mg/min]") + ylim(0, 40) +
    
    theme(text=element_text(size=14))
  
  
  plot<-plot+ggtitle(paste0("Model 2: AE, different temperature levels for SL #",i ))
  
   
  
  print(plot)
  
  #### working way how to save in a good resoultion with converting 
  ggsave(plot, filename = paste("dff_temp_SL",i,".tiff", sep=""), dpi=600,   # ES&T. 300-600 at PLOS One,
         compression = "lzw",  # PLOS One 44MB w/out compression, 221KB w/compression
         family="Arial",  # ggplot default.  Could set to others depending on journal
         type="cairo")  # PLOS One 44MB w/out compression, 221KB w/compression)
  dev.off() 
  
}

 

