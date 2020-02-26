# Plotting Event Rate vs. Time for FCS Files (Flow or Mass / CyTOF)


#########################################################
### Installing and loading required packages
#########################################################

if (!require("svDialogs")) {
  install.packages("svDialogs", dependencies = TRUE)
  library(svDialogs)
}

if (!require("flowCore")) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("flowCore")
}

if (!require("tidyverse")) {
  install.packages("tidyverse", dependencies = TRUE)
  library(tidyverse)
}

if (!require("reshape2")) {
  install.packages("reshape2", dependencies = TRUE)
  library(reshape2)
}

if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}

if (!require("tcltk")) {
  install.packages("tcltk", dependencies = TRUE)
  library(tcltk)
}

#########################################################
### Script starts here
#########################################################

# Hack to make plots nicer in RStudio Plots window (but not zoom)
trace(grDevices:::png, quote({
  if (missing(type) && missing(antialias)) {
    type <- "cairo-png"
    antialias <- "subpixel"
  }
}), print = FALSE)

# Clear environment
rm(list = ls(all = TRUE))

# Data Import from file chosen by user

#library(svDialogs) # Moved to top 
# Get user input for file
testfile <- dlg_open()
# Convert to string value
testfile <- capture.output(testfile)[7]
{
  
  if ((testfile)=="character(0)")
    stop("File input cancelled")
  
  #Remove invalid characters from file input location
  testfile <- gsub("[\"]","",testfile)
  testfile<-substring (testfile,5)
  
  #Set file and directory
  filename <- basename (testfile)
  dir <- dirname (testfile)
  
  # Set working directory accoding to file chosen
  setwd(dir)
  
  #library(flowCore) #Moved to top
  
  # this read.FCS() function imports the flow data:
  raw_fcs<-read.FCS(filename, alter.names = TRUE)
  
  # Preparation work for arcsinh transform (columns is also used later for naming changes)
  # Create list of parameters
  columns<-colnames(raw_fcs)
  # Remove "Time" column to avoid it being transformed
  columns<-setdiff(columns,"Time")
  # Remove "Cell_Length" and Gaussians column to avoid it being transformed
  columns<-setdiff(columns,"Event_length")
  columns<-setdiff(columns,"Cell_length")
  columns<-setdiff(columns,"Center")
  columns<-setdiff(columns,"Offset")
  columns<-setdiff(columns,"Width")
  columns<-setdiff(columns,"Residual")
  ## Remove FSC and SSC
  removefscssc<-grep("FSC|SSC",columns,value=TRUE)
  columns<-columns[! columns %in% removefscssc]
  
  # Read data into a data frame
  FCSDATA <- as.data.frame(exprs(raw_fcs))

  #Remove unnecessary parameter text
  names(FCSDATA)[-1] <- sub("Di", "", names(FCSDATA)[-1])
  names(FCSDATA)[-1] <- sub("Dd", "", names(FCSDATA)[-1])
  # Create list of channel / parameter descriptions 
  params<-parameters(raw_fcs)[["desc"]]
  # Replace parameters with descriptions, keeping things like Time, Event Length unchanged
  colnames(FCSDATA)[!is.na(params)] <- na.omit(params)
  
  # Determine whether data is CyTOF or Flow by presence of FSC
  # isflow will be 0 for a CyTOF or greater than 1 if flow
  isflow <-sum(grep("FSC",colnames(FCSDATA)))
  # Determine whether data is pre CyTOF 3 (Helios) by presence of "Cell_length", rather than "Event_length"
  isCyTOF2 <-sum(grep("Cell_length",colnames(FCSDATA)))
  
  ## Remove Time, Event_Length & Gaussian Parameters
  removecolumns <- c("Event_length", "Center", "Offset", "Width", "Residual", "Cell_length")
  FCSDATA <- FCSDATA[,!(names(FCSDATA) %in% removecolumns)]
  
  
  ## Remove FSC and SSC
  # library(tidyverse) # Moved to top
  FCSDATA <- FCSDATA %>% select(-contains("FSC"))
  FCSDATA <- FCSDATA %>% select(-contains("SSC"))
  
  # Get number of cell events (based on "191 or 193" - i.e. Iridium)
  if(isflow==0){
    # Note that this only works correctly because "Time" has been removed by a previous step - otherwise the position would be wrong.
    irpos <- grep("191|193",colnames(FCSDATA))
    # Divide by two because we're looking for two columns
    cellevents <- sum(FCSDATA[,irpos]!=0)/2
    kcellevents <- round(cellevents/1000,0)
  }else{
    cellevents <- nrow(FCSDATA)
    kcellevents <- round(cellevents/1000,0)
  }
  
  #For converting FCS time to mins - flow uses 10ms units, CyTOF uses ms
  if(isflow>0){
    div = (60*100)
  }else{
    div = (60*1000)
  }
  
  # Get position of Time in dataset (to handle Flow and CyTOF)
  TimePos <- which(colnames(FCSDATA)=="Time")
  
  # Find total acquisition time
  maxtime<-round(max(FCSDATA[,TimePos])/div,1)
  
  # Now that we have the total time, we can calculate the number of cell events/sec
  eventspersec <- round(cellevents/maxtime/60,0)
  
  
  # Would be nice to know cell events, as well as all events
  if(isflow==0){
    FCSDATA <- as.data.frame(FCSDATA[,c(TimePos,irpos)])
    # Create cellevents column
    FCSDATA$CellEvents <- FCSDATA$`191Ir`> 0 | FCSDATA$`193Ir`> 0
    # Remove other columns
    FCSDATA <- as.data.frame(FCSDATA[,c(TimePos,4)])
  }

  if(isflow!=0){
    # Keep only the time column - we don't care about the others now
    FCSDATA <- as.data.frame(FCSDATA[,TimePos])
    colnames(FCSDATA) <- "Time"
  }
  
  # Sort by Time
  FCSDATA <- as.data.frame(FCSDATA[order(FCSDATA[,1]),])
  
  # Create data for plot, accounting for the downsampling
  # Set resolution based on length of dataset - longer than 1/2 hour = 10 sec resolution
  if (isflow==0){
    Resolution <- ifelse(maxtime>30,10,1)
  }else{
    # Flow tends to be shorter acquisition times, so use >1 mins = 1 sec resolution
    Resolution <- ifelse(maxtime>1,1,0.1)
  }
  
  EventsPerSecSeries <- NULL
  for (i in seq(div/60*Resolution,max(FCSDATA[,1]),div/60*Resolution)){
    EventsPerSecSeries[length(EventsPerSecSeries)+1] <- sum(FCSDATA[,1] < i & FCSDATA[,1] > (i-div/60*Resolution)) 
  }
  # Convert to seconds
  EventsPerSecSeries <- EventsPerSecSeries/Resolution
  
  
  # And for cell events if CyTOF
  if (isflow==0){
    CellsPerSecSeries <- NULL
    for (i in seq(div/60*Resolution,max(FCSDATA[,1]),div/60*Resolution)){
      CellsPerSecSeries[length(CellsPerSecSeries)+1] <- sum(FCSDATA[,1] < i & FCSDATA[,1] > (i-div/60*Resolution)) * sum(FCSDATA[,1][FCSDATA[,2]==TRUE] < i & FCSDATA[,1][FCSDATA[,2]==TRUE] > (i-div/60*Resolution)) / sum(FCSDATA[,1] < i & FCSDATA[,1] > (i-div/60*Resolution))
    
    }
    # Convert to seconds
    CellsPerSecSeries <- CellsPerSecSeries/Resolution
  }
  
  
  EventsForPlot <- as.data.frame(cbind("Events" = EventsPerSecSeries,"Time" = c(1:length(EventsPerSecSeries))*Resolution))
  if (isflow==0){
    EventsForPlot <- cbind(EventsForPlot,CellsPerSecSeries)
  }

  ## Plot x as Time and Y as Events/sec
  EventsPlot <- ggplot(EventsForPlot,aes(x = Time, y = Events)) +
    geom_line(colour="blue", alpha=0.5)+
    stat_smooth(geom="line",method="loess", colour="blue", alpha=0.5, span=0.25)+
    # Hide legend
    theme(legend.position = "none") +
    # Change Y axis label
    ylab("Events/sec") +
    # Change X axis label
    xlab("Time (sec)")+
    ggtitle(filename) +
    # Add mean events/sec label
    geom_label(colour="black",
               fontface="bold",
               size=2.5,
               alpha=0.5,
               mapping=aes(max(EventsForPlot$Time)/2,(max(EventsForPlot$Events))*1.05,
                           label=paste("Cell events/sec =",eventspersec, ", Total = ",kcellevents, " k")))+
    # Average
    geom_hline(yintercept = mean(EventsForPlot$Events), linetype="dashed", colour="blue")+
    # +/- 2SD
    geom_hline(yintercept = mean(EventsForPlot$Events) + 2 * sd(EventsForPlot$Events), linetype="dashed", colour="blue")+
    geom_hline(yintercept = mean(EventsForPlot$Events) - 2 * sd(EventsForPlot$Events), linetype="dashed", colour="blue")+
    # Axes
    coord_cartesian(expand=FALSE)+
    ylim(0,max(EventsForPlot$Events*1.1))
  
  if (isflow==0){
    # For Cell Events
    EventsPlot +
      geom_line(data=EventsForPlot,aes(x=Time, y=CellsPerSecSeries), colour="red", alpha=0.5)+
      stat_smooth(geom="line",method="loess", colour="red", alpha=0.5, se=FALSE, span=0.25)+
      # Average
      geom_hline(yintercept = mean(na.omit(EventsForPlot$CellsPerSecSeries)), linetype="dashed", colour="red", alpha=0.5)+
      # +/- 2SD
      geom_hline(yintercept = mean(na.omit(EventsForPlot$CellsPerSecSeries)) + 2 * sd(na.omit(EventsForPlot$CellsPerSecSeries)), linetype="dashed", colour="red", alpha=0.5)+
      geom_hline(yintercept = mean(na.omit(EventsForPlot$CellsPerSecSeries)) - 2 * sd(na.omit(EventsForPlot$CellsPerSecSeries)), linetype="dashed", colour="red", alpha=0.5)
  }else{
      # Only plot events, not events + cellevents
      EventsPlot
     }
  
} # End of file cancel loop