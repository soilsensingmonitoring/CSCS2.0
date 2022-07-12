## CSCS2.0
  
## Create a function for allocating an unknown soil to existing soils in the Comprehensive Soil Classification System (CSCS2.0)
  
# Input: the full names of the input excel spreadsheet;
# Newname: "old" will allocate the unknown soil profiles to CSCS using the old names (e.g. ST_HAPLUSTEPTS); 
           # "new" will allocate the unknown soil profiles using the new names (e.g. BAJESOZEM);
# Number: the number of closest soil taxa to display;
# Plot: if the resulting plots will be saved as png files to the folder;
# Label: if the new soil taxa will be labelled on the PC plot.


CSCS_allocation <- function (Input, Newname, Number, Plot, Label) 
{
    ## Turn off warnings globally as some packages give warnings that are not necessary
    options(warn=-1)
    
    ## But note that turning off warning messages globally might not be a good idea.
    ## To turn warnings back on, use
    ## options(warn=0)
    
    
    ## Create a function to install required packages
    pkgTest <- function(x)
    {
      if (!require(x,character.only = TRUE))
      {
        install.packages(x,dep=TRUE)
        if(!require(x,character.only = TRUE)) stop("Package not found")
      }
    }
    
    ## names of required packages 
    lib <- c("matrixStats", "readxl", "devtools", "aqp", "plyr", "dendextend", "ggplot2", "Cairo")
    
    ## Install packages
    for (i in 1:length(lib))
    {
      pkgTest(lib[i])
    }
    
    
    ## Install "ithir"
    if (!require(ithir))
    {
      install_bitbucket("brendo1001/ithir/pkg")
    }
    
    ## Load all required packages
    
    library(matrixStats)
    library(readxl)
    library(devtools)
    library(ithir)
    library(plyr)
    library(dendextend)
    library(ggplot2)
    library(aqp)
    library(Cairo)
    
    ## Clear console
    cat("\014")
    
    
    ## Load existing data
    load(file = "data_cscs_2.0.RData")
    
    ## Note: The existing Universal (Comprehensive) Soil Classification System (CSCS) contains the 414 soil properties of selected 493 soil taxa from ST, ASC, NZ and WRB
    ## We used the existing soil dataset as a reference to generate a mean and SD of the 414 soil properties
    ## mean and SD for ice content (145:162 elements) is set to 0 and 1, respectively.
    
    CSCS_mean
    CSCS_SD
    
    colnames(cscs_old.pc) ## 414 PCs of the existing 493 soil profiles
    str(Base.data.combined.comps) ## PC formula of the existing 493 soil profiles, which can be used to get the PCs of any new soil data
    
    
    ############# Input unknown data
    
    df_orginal <- read_excel(Input)
    var.name <- colnames(df_orginal)[-c(1:3)]
    head(df_orginal)
    
    
    
    
    
    
    
    
    
    
    
    
    ############# Spline fitting input data to the same depth
    profileID <- unique(df_orginal$Soil.ID)
    
    ## Create a spline fitting output dataframe, first column is the profile ID, next 414 columns are the fitted soil properties at different depths
    df_fitted <- matrix(NA, nrow = length(profileID), ncol = 415)
    df_fitted <- as.data.frame(df_fitted)
    colnames(df_fitted) <- var_names
    df_fitted$oldnames <- profileID
    df_fitted[1,]
    
    ## We will use a uniform lambda for fitting all the soil profiles
    lambda <- 0.01
    
    ## For loop for spline fitting
    for (i in 1:length(profileID))
    {  
      ## Subset the profile in turn
      df_subset<-df_orginal[df_orginal$Soil.ID==profileID[i],]
      df_subset
      
      ## Get the ID of the profile
      temp <- df_fitted[df_fitted$oldnames == profileID[i],]$oldnames
      temp
      
      
      
      ## For each of the variable, use ea_spline function to fit the soil properties to the standard depth
      ## As we have used i in the outer loop, we cannot use the same index. We use j for the inner loop
      for (j in 1:length(var.name))
      { 
        ## Convert the subset data into a format for spline fitting and Subset non-NAs
        df_subset2 <- df_subset[!is.na(df_subset[,var.name[j]]),]
        df_subset2
        
        
        ## If the profile has more than or equal to three non-NAs, then spline-fitting 
        if (nrow(df_subset2[,var.name[j]]) >= 3)   
        {  
          depths(df_subset2) <- Soil.ID ~ Upper.Boundary.cm + Lower.Boundary.cm
          
          df_Fit <- ea_spline(df_subset2, var.name = var.name[j], d= t(c(0,5,10,15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80 , 90, 100, 110, 130, 150)), lam = lambda, vlow = 0, show.progress=FALSE)
          df_Fit_df <- as.data.frame(df_Fit$harmonised)[, -c(1, 20)] ## extract the fitted soil properties at all the standard depths
          df_Fit_df
        }
        
        
        ## if all depths are NAs for this soil property, then gap-filling with mean values from CSCS mean values  
        else
          
        {   
          index_first <- (j-1)*18+1
          index_last  <- (j-1)*18+18
          df_Fit_df <- as.data.frame(t(CSCS_mean[index_first:index_last]))
          df_Fit_df
        }
        
        temp <- cbind(temp,  df_Fit_df)
      }
      
      df_fitted[df_fitted$oldnames == profileID[i],] <-  temp
      
    }
    
    
    
    
    
    
    
    
    
    
    
    
    ## Fill bedrock values (below the measured depths) for all the soil profiles
    ## Get the soil values of already filled profiles (removing ID column)
    temp_bedrock <- df_fitted[,-1]
    
    ## For each row (soil profile data), find the NAs and then fill them with the mean values of the reference data
    for (i in 1:nrow(temp_bedrock))
    { 
      temp_bedrock[i, is.na.data.frame(temp_bedrock[i,])] <- CSCS_mean[is.na.data.frame(temp_bedrock[i,])]
      
    }
    
    ## Combine the column (profile) ID with the filled values
    df_fitted <- cbind(df_fitted[,1], temp_bedrock)
    colnames(df_fitted) <- var_names
    write.csv(df_fitted, "Spline_fitted_profiles.csv", row.names=F)
    
    ## Clear console
    cat("\014")
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    #### Project new data onto the CSCS
    ## Get the splined-fitted values of all the unknown profiles
    newdata <- df_fitted[,-1]
    newdata2 <- newdata
    
    
    ## Scale the splined-fitted values of all the unknown profiles using the Mean and SD of the reference data
    for (i in 1: length(CSCS_mean)) 
      for (j in 1:nrow(newdata))
      {
        newdata2[j,i]<-(newdata[j,i]-CSCS_mean[i,1])/CSCS_SD[i,1];
      }
    
    
    ## Predict the scaled data to the PC spaces of the existing USCS
    Proj.newdata.x <- predict(Base.data.combined.comps, as.data.frame(newdata2))
    row.names(Proj.newdata.x) <- df_fitted[,1]
    row.names(Proj.newdata.x)
    
    ## Clear console
    cat("\014")
    
    if (Plot == T)
    {
      ## Plot the unknown soil profiles in the USCS in PC space
      Cairo::Cairo(
        30, #length
        30, #width
        file = "PC_CSCS_New_Soil.png",
        type = "png", #tiff
        bg = "transparent", #white or transparent depending on your requirement 
        dpi = 600,
        units = "cm" #you can change to pixels etc 
      )
      plot(cscs_old.pc[1:100,1:2], col="black", pch=19, cex=1.1, xlim=c(-100,60),ylim=c(-100,120))  ## ST 
      points(x=cscs_old.pc[101:104, 1], y=cscs_old.pc[101:104, 2], col="blue", pch=19, cex=1.1)     ## WRB
      points(x=cscs_old.pc[105:231, 1], y=cscs_old.pc[105:231, 2], col="orange", pch=19, cex=1.1)   ## ASC
      points(x=cscs_old.pc[232:274, 1], y=cscs_old.pc[232:274, 2], col="cyan", pch=19, cex=1.1)     ## NZ
      points(x=cscs_old.pc[275:320, 1], y=cscs_old.pc[275:320, 2], col="grey", pch=19, cex=1.1)     ## Fr
      points(x=cscs_old.pc[321:375, 1], y=cscs_old.pc[321:375, 2], col="yellow", pch=19, cex=1.1)     ## Ru
      points(x=cscs_old.pc[376:383, 1], y=cscs_old.pc[376:383, 2], col="pink", pch=19, cex=1.1)     ## Br
      points(x=cscs_old.pc[384:398, 1], y=cscs_old.pc[384:398, 2], col="brown", pch=19, cex=1.1)     ## Kr
      points(x=Proj.newdata.x[,1], y=Proj.newdata.x[,2],col="red", pch=19, cex=1.1)           ## new unknown
      
      if (Label == T)
      {
        text(x=Proj.newdata.x[,1], y=Proj.newdata.x[,2], labels = df_fitted[,1], data= Proj.newdata.x, cex=1, font=2, pos=1, col ='red')
      }
      dev.off()
      
      
    }
    
    ## Get the names of the USCS profiles and the unknown profiles
    row.names(Proj.newdata.x) <- df_fitted[,1]
    
    if (Newname == "old")
    {row.names(cscs_old.pc) <- old_names } 
    if (Newname == "new")
    {row.names(cscs_old.pc) <- new_names }
    
    
    ## Output the PCs of the USCS and the unknown profiles
    fr.pcs <- rbind(cscs_old.pc, Proj.newdata.x)
    write.csv(fr.pcs,"pcs_CSCS_new_profiles.csv",row.names = T)
    write.csv(Proj.newdata.x,"pcs_new_profiles.csv",row.names = T)
    
    
    
    ############### Hard classification   
    ## Calculate distance matrix
    hard.distance.unknownx <- matrix(NA, nrow = nrow(Proj.newdata.x), ncol = nrow(cscs_old.pc) ) 
    
    for (i in 1:nrow(Proj.newdata.x))
      for (j in 1:nrow(cscs_old.pc))
      {  
        hard.distance.unknownx[i,j]<-as.numeric(dist(rbind(Proj.newdata.x[i,],cscs_old.pc[j,]),method = "euclidean"))##get distance index, not zero
      }
    
    

    
    
    ## For each of the new profiles, get the names of all the USCS profiles
    hard.distance.unknown.sorted.x<-matrix(NA, nrow = nrow(Proj.newdata.x) , ncol = nrow(cscs_old.pc))
    hard.distance.unknown.sorted.dis<-matrix(NA, nrow = nrow(Proj.newdata.x) , ncol = nrow(cscs_old.pc))
    
    for (i in 1: nrow(Proj.newdata.x))
    {
      hard.distance.unknown.sorted.x[i,] <- row.names(cscs_old.pc)
    }
    
    
    ## For each of the new profiles, sort the names based on the closet N taxa in the USCS
    for (i in 1: nrow(Proj.newdata.x))
    {
      hard.distance.unknown.sorted.x[i,] <- hard.distance.unknown.sorted.x[i, order(hard.distance.unknownx[i,])]
      hard.distance.unknown.sorted.dis[i,] <- hard.distance.unknownx[i, order(hard.distance.unknownx[i,])]
    }
    
    ## Get the closet N taxa names in the USCS
    hard.distance.unknown.sorted.x[,1:Number]
    
    allocated.names <- cbind(row.names(Proj.newdata.x), hard.distance.unknown.sorted.x[,1:Number])
    colnames(allocated.names) <- c("ProfileID", paste("Closest_taxon_", 1:Number, sep = ""))
    
    ## Get the closet N taxa distances in the USCS
    hard.distance.unknown.sorted.dis[,1:Number]
    
    allocated.distance<-cbind(row.names(Proj.newdata.x), hard.distance.unknown.sorted.dis[,1:Number])
    colnames(allocated.distance) <- c("ProfileID", paste("Closest_distance_", 1:Number, sep = ""))
    
    
    write.csv(allocated.names, "allocated.names.csv",row.names = F)
    write.csv(allocated.distance, "allocated.distance.csv",row.names = F)
    
    
    
    
    #### Dendrogram
    tree.un <- hclust(as.dist(as.matrix(dist(cscs_old.pc))), method="ward.D2")
    
    if (Plot == T)
    {
      
      Cairo::Cairo(
        120, #length
        40, #width
        file = "Tree_CSCS.png",
        type = "png", #tiff
        bg = "transparent", #white or transparent depending on your requirement 
        dpi = 600,
        units = "cm" #you can change to pixels etc 
      )
      plot(tree.un, hang = -1, cex = 0.6)
      dev.off()
      
      # Rectangle dendrogram using ggplot2
      
      dend <- tree.un %>%
        as.dendrogram %>%
        set("branches_k_color", k=8) %>% set("branches_lwd", 1) %>%
        set("labels_colors") %>% set("labels_cex", c(.3, 0.5))
      
      
      
      
      Cairo::Cairo(
        120, #length
        40, #width
        file = "dendrogram_CSCS.png",
        type = "png", #tiff
        bg = "transparent", #white or transparent depending on your requirement 
        dpi = 600,
        units = "cm" #you can change to pixels etc 
      )
      
      plot(dend)
      
      dev.off()
      
    }  
    
    ## Print the allocated names and distances
    cat("\014")
    hard.distance.unknown.sorted.x[,1:Number]
    hard.distance.unknown.sorted.dis[,1:Number]
    
    return(allocated.names)
  }
  
  
  
## Input a new soil dataset and allocate the soil profiles to existing soils and plot them
  
output <- CSCS_allocation(Input = "Input_UK.xlsx", Newname = "new", Number = 5, Plot = T, Label = F)
  
output
