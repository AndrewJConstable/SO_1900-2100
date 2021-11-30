library(ncdf4)
library(tidyverse)
library(raster)
library(here)
library(lubridate) # for date functions
library(scales)  # generating scales for plotting etc.
library(ggplot2)
library(dplyr)
library(latex2exp) # for writing the degree Celsius
library(sf)
library(magick)
library(ggforce)


########################
# Files ################

pathSST<-"/Users/andreworca/Desktop/_wAAD/_r/MEASO/ACCESS/SO_1900-2100/Temperature/Means/"
filesSST<-c(
  "MonthlyMeans_1900-1909.nc"
  ,"MonthlyMeans_1910-1919.nc"
  ,"MonthlyMeans_1920-1929.nc"
  ,"MonthlyMeans_1930-1939.nc"
  ,"MonthlyMeans_1940-1949.nc"
  ,"MonthlyMeans_1950-1959.nc"
  ,"MonthlyMeans_1960-1969.nc"
  ,"MonthlyMeans_1970-1979.nc"
  ,"MonthlyMeans_1980-1989.nc"
  ,"MonthlyMeans_1990-1999.nc"
  ,"MonthlyMeans_2000-2009.nc"
  ,"MonthlyMeans_2010-2014.nc"
  ,"MonthlyMeans_2015-2024.nc"
  ,"MonthlyMeans_2025-2034.nc"
  ,"MonthlyMeans_2035-2044.nc"
  ,"MonthlyMeans_2045-2054.nc"
  ,"MonthlyMeans_2055-2064.nc"
  ,"MonthlyMeans_2065-2074.nc"
  ,"MonthlyMeans_2075-2084.nc"
  ,"MonthlyMeans_2085-2094.nc"
  ,"MonthlyMeans_2095-2100.nc"
) # end file list
SSTfiles<-as.data.frame(list(
  path = sapply(filesSST,function(p,d) paste(d,p,sep=""), pathSST)
  ,year0 = as.numeric(sapply(filesSST,function(p) substr(p,14,17)))
  ,year1 = as.numeric(sapply(filesSST,function(p) substr(p,19,22)))
))
SSTfilesYrs<-do.call(rbind,lapply(seq(1,nrow(SSTfiles),1)
                                  ,function(i,f){
                                    cbind(seq(f[i,2],f[i,3],1)
                                          ,rep(i,(f[i,3]-f[i,2]+1)))},SSTfiles))

pathSIC<-"/Users/andreworca/Desktop/_wAAD/_r/MEASO/ACCESS/SO_1900-2100/SIC/Means/"
filesSIC<-c(
  "MonthlyMeans_SIC_1850-2014.nc"
  ,"MonthlyMeans_SIC_2015-2100.nc"
) # end file list

SICfiles<-as.data.frame(list(
  path = sapply(filesSIC,function(p,d) paste(d,p,sep=""), pathSIC)
  ,year0 = as.numeric(sapply(filesSIC,function(p) substr(p,18,21)))
  ,year1 = as.numeric(sapply(filesSIC,function(p) substr(p,23,26)))
))
SICfilesYrs<-do.call(rbind,lapply(seq(1,nrow(SICfiles),1)
                                  ,function(i,f){
                                    cbind(seq(f[i,2],f[i,3],1)
                                          ,rep(i,(f[i,3]-f[i,2]+1)))},SICfiles))

pathGWL<-"/Users/andreworca/Desktop/_wAAD/_r/MEASO/ACCESS/SO_1900-2100/GWL/"
filesGWL<-c(
  "tas_Amon_ACCESS-ESM1-5_r1i1p1f1_gn_185001-230012_gwl20.nc"
  ,"tas_Amon_ACCESS-ESM1-5_r1i1p1f1_gn_185001-230012_gwl30.nc"
  ,"tas_Amon_ACCESS-ESM1-5_r1i1p1f1_gn_185001-230012_ymonmean.nc"
) # end file list 

############################################
# input Global Warming Level (GWL)
ncfname <- paste(pathGWL,filesGWL[2],sep="")
ncin <- nc_open(ncfname)
print(ncin)

# read in dimensions

lon<-ncvar_get(ncin,"lon")
nlon<-dim(lon)
lat<-ncvar_get(ncin,"lat")
nlat<-dim(lat)
# depth<-ncvar_get(ncin,"lev")
# ndepth<-dim(depth)

t<-ncvar_get(ncin,"time")
nt<-dim(t)
tunits<-ncatt_get(ncin,"time","units") # days since 1850-01-01
Year<-year(as.Date(t,origin="1850-01-01"))

# read in variable array
ncVar<-"tas"
ncData<-ncvar_get(ncin,ncVar)
NAvalue<-ncatt_get(ncin,ncVar,"_FillValue")[[2]]
ncVarUnits<-ncatt_get(ncin,ncVar,"units")
nc_close(ncin)

rGWL <- brick(ncfname, stopIfNotEqualSpaced = F)
#reset minima and maxima of dimensions of each raster in order to enable use of 'area'
xmin(rGWL)<-0
xmax(rGWL)<-360
ymin(rGWL)<-(-90)
ymax(rGWL)<-90

# calculate weighted means
CellArea<-area(rGWL[[1]])
rGWL_wt_vals<-rGWL*CellArea
GWL<-cellStats(rGWL_wt_vals,sum)/cellStats(CellArea,sum)
GWL<-GWL-GWL[1]
plot(Year,GWL,type="l") # plot GWL as a check

#Reproject to South Polar Stereographic (EPSG:3976)
#Define projection
south_stereo <- "+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# latitude circle for plotting
pLat<- -60
pLonInc<-0.1
pCircle<-st_linestring(cbind(seq(-180,180,pLonInc),rep(pLat,(360/pLonInc+1))))
pCircle<-st_sfc(pCircle,crs=4326)
pCircle<-st_transform(pCircle,south_stereo)
pCircle<-as.data.frame(st_coordinates(pCircle))
#scale circle (compared to output - hard wire fix of scaling issue)
  # xMin<-min(aData[,"x"])
  # xCmin<-min(pCircle[,"x"])
  # scaleC<- #(xMin*30/60/xCmin-1)*2.5+1
pCircle<-1.293*pCircle

############################################
# Routine
############################################

SIC_file<-0 # SIC file not loaded yet
SST_file<-0 # SST file not loaded yet

doYears<-seq(2010,2014,1)
doMonths<-seq(1,12,1)

printPlot<-TRUE
savePlot<-TRUE

outputDir<-"/Users/andreworca/Desktop/_wAAD/_r/MEASO/ACCESS/Output Animation"
################################
#do loop

for (y in doYears){
  for (m in doMonths){
    
    # determine which file the year is in
    whichSSTfile<-SSTfilesYrs[which(SSTfilesYrs[,1]==y),2]
    if(SST_file!=whichSSTfile){
      SST_file<-whichSSTfile
      sst<-brick(SSTfiles[SST_file,1], stopIfNotEqualSpaced = F)
      sstDate<-as.vector(sapply(names(sst),function(d) paste(substr(d,2,5),substr(d,7,8),substr(d,10,11),sep="-")))
      sstYears<-year(sstDate)
      sstMonths<-month(sstDate)
      sst_reproj <- projectRaster(sst, crs =south_stereo) # reproject to polar
    }
    
    whichSICfile<-SICfilesYrs[which(SICfilesYrs[,1]==y),2]
    if(SIC_file!=whichSICfile){
      SIC_file<-whichSICfile
      sic<-brick(SICfiles[SIC_file,1], stopIfNotEqualSpaced = F)
      sicDate<-as.vector(sapply(names(sic),function(d) paste(substr(d,2,5),substr(d,7,8),substr(d,10,11),sep="-")))
      sicYears<-year(sicDate)
      sicMonths<-month(sicDate)
      sic_reproj <- projectRaster(sic, crs =south_stereo) # reproject to polar
    }
    
    # subset raster
    p<-which(sstYears==y & sstMonths==m)
    pName<-names(sst_reproj)[p]
    
    sst_df_reproj <- as.data.frame(sst_reproj[[pName]], xy = T)
    
    aData<-sst_df_reproj
    maxSST<-10
    
    aData[(!is.na(aData[,3]) & aData[,3]>maxSST),3] <- maxSST
    
    # colour ramp for just sst
    pColours <- c("#2c7bb6","#f7f7f7","#ca0020","black","black")
    pValues  <- rescale(c(-2,2,4,5,maxSST),to = c(0,1), from=c(-2,maxSST))
    
    minSIC<-25
    sic_df_reproj <- as.data.frame(sic_reproj[[pName]], xy = T)
    mData<-merge(aData,sic_df_reproj,by=c("x","y"))
    useSIC <- (!is.na(mData[,4]) & mData[,4]>=minSIC & mData[,4]<=100)
    aData<-mData[,c(1:3)]
    aData[useSIC,3]<-mData[useSIC,4]
    names(aData)<-c("x","y",pName)
    pColours <- c("#2c7bb6","#f7f7f7","#ca0020","black","black","blue","darkblue","darkblue")
    pValues  <- rescale(c(-2,2,4,5,(minSIC-0.1),minSIC,80,100),to = c(0,1), from=c(-2,100))
    
    
    pYear<-year(sstDate[p])
    pGWL<-as.vector(GWL[which(Year==pYear)])
    pGWLcolour<-ifelse(pGWL<=1.5,"blue",ifelse(pGWL>1.5 & pGWL <=2,"red","darkred")) 
    pGWL_txt<-paste(format(round(pGWL, 1), nsmall = 1),"\u00B0C",sep="")
    
    # location of annotations
    names(pCircle)<-names(aData)
    
    
    aPlot<-  ggplot(aData,aes_string("x", "y", fill = pName)) +
      geom_raster()+
      scale_fill_gradientn(
        colours=pColours
        ,values=pValues
        ,na.value = "white")+
      theme_minimal()+
      labs(fill = "Scale") +
      theme(legend.position="none",aspect.ratio=1
            ,axis.text.x = element_blank()
            ,axis.text.y = element_blank()
            ,axis.ticks = element_blank()
            ,axis.title.x = element_blank()
            ,axis.title.y = element_blank()
            ,plot.title = element_text(size = 18, face = "bold",colour="#044c73")
      ) + # end theme
      ggtitle(paste(pYear,month(sstDate[p],label=TRUE,abbr=FALSE),sep=" : ")) +
      annotate("text", x=max(aData[,"x"])/10,y=0, 
               label=pGWL_txt,size = 6,color = pGWLcolour) +
      annotate("text", x=max(aData[,"x"])/10,y=max(aData[,"y"])/7,
               label="GWL",size = 3,color = pGWLcolour)

    names(pCircle)<-names(aData)
    aPlot<-aPlot+geom_polygon(data=pCircle,colour="white",fill=NA)
    
    if(printPlot) print(aPlot)
    if(savePlot) ggsave(paste(outputDir,"/",pName,".png",sep=""))
    
    #####################
  } # end do month
} # end do year
###############################################################################
###############################################################################
# create animation
outputDir<-"/Users/andreworca/Desktop/_wAAD/_r/MEASO/ACCESS/Output animation"

imgs <- list.files(outputDir, full.names = TRUE)
img_list <- lapply(imgs, image_read)

## join the images together
img_joined <- image_join(img_list)

## animate at 2 frames per second
img_animated <- image_animate(img_joined, fps = 20)

## view animated image
img_animated

## save to disk
image_write(image = img_animated,
            path = "COP26.gif")


###############################################################################
# generating mean raster for a given month for a range of GWL
###############################################################################
vGWL<-c(0,1.5,2,4)
GWLmonth  <- 2

for(g in 1:length(vGWL)){

  pGWL<-vGWL[g]  # value to be plotted
  GWLrange  <- c((pGWL-0.3),pGWL)
  pGWLcolour<-ifelse(pGWL<=1.5,"blue",ifelse(pGWL>1.5 & pGWL <=2,"red","darkred")) 

GWLyears  <- Year[which(GWL>=GWLrange[1] & GWL<=GWLrange[2])] # do as sequence of years from 1900 to 1950 for first one
useSSTyear<- SSTfilesYrs[,1] %in% GWLyears

SIC_file<-0 # SIC file not loaded yet
SST_file<-0 # SST file not loaded yet

doYears<-SSTfilesYrs[useSSTyear,1]
if(pGWL==0) doYears<-c(1900:1950)
m<-GWLmonth
printPlot<-TRUE
savePlot<-TRUE

outputDir<-"/Users/andreworca/Desktop/_wAAD/_r/MEASO/ACCESS"

################################
#do loop to combine rasters into a stack

rStack_sst<-NULL
rStack_sic<-NULL

for (y in doYears){
  
  # determine which file the year is in
  whichSSTfile<-SSTfilesYrs[which(SSTfilesYrs[,1]==y),2]
  if(SST_file!=whichSSTfile){
    SST_file<-whichSSTfile
    sst<-brick(SSTfiles[SST_file,1], stopIfNotEqualSpaced = F)
    sstDate<-as.vector(sapply(names(sst),function(d) paste(substr(d,2,5),substr(d,7,8),substr(d,10,11),sep="-")))
    sstYears<-year(sstDate)
    sstMonths<-month(sstDate)
  }
  
  whichSICfile<-SICfilesYrs[which(SICfilesYrs[,1]==y),2]
  if(SIC_file!=whichSICfile){
    SIC_file<-whichSICfile
    sic<-brick(SICfiles[SIC_file,1], stopIfNotEqualSpaced = F)
    sicDate<-as.vector(sapply(names(sic),function(d) paste(substr(d,2,5),substr(d,7,8),substr(d,10,11),sep="-")))
    sicYears<-year(sicDate)
    sicMonths<-month(sicDate)
  }
  
  # subset raster
  p<-which(sstYears==y & sstMonths==m)
  pName<-names(sst)[p]
  
  if(is.null(rStack_sst)){
    rStack_sst<-sst[[pName]]
    rStack_sic<-sic[[pName]]
  } else {
    rStack_sst<-stack(rStack_sst,sst[[pName]])
    rStack_sic<-stack(rStack_sic,sic[[pName]])
  }
} # end doYears

# mean of cells
rMean_sst<-mean(rStack_sst)
rMean_sic<-mean(rStack_sic)


# plotting        

sst_reproj <- projectRaster(rMean_sst, crs =south_stereo) # reproject to polar
sic_reproj <- projectRaster(rMean_sic, crs =south_stereo) # reproject to polar


sst_df_reproj <- as.data.frame(sst_reproj, xy = T)

aData<-sst_df_reproj

aData[(!is.na(aData[,3]) & aData[,3]>maxSST),3] <- maxSST

# colour ramp for just sst
pColours <- c("#2c7bb6","#f7f7f7","#ca0020","black","black")
pValues  <- rescale(c(-2,2,4,5,maxSST),to = c(0,1), from=c(-2,maxSST))

sic_df_reproj <- as.data.frame(sic_reproj, xy = T)
mData<-merge(aData,sic_df_reproj,by=c("x","y"))
useSIC <- (!is.na(mData[,4]) & mData[,4]>=minSIC & mData[,4]<=100)
aData<-mData[,c(1:3)]
aData[useSIC,3]<-mData[useSIC,4]
names(aData)<-c("x","y",pName)
pColours <- c("#2c7bb6","#f7f7f7","#ca0020","black","black","blue","darkblue","darkblue")
pValues  <- rescale(c(-2,2,4,5,(minSIC-0.1),minSIC,80,100),to = c(0,1), from=c(-2,100))

MonthName<-unique(sstMonths)[m]

aPlot<-  ggplot(aData,aes_string("x", "y", fill = pName)) +
  geom_raster()+
  scale_fill_gradientn(
    colours=pColours
    ,values=pValues
    ,na.value = "white")+
  theme_minimal()+
  labs(fill = "Scale") +
  theme(legend.position="none",aspect.ratio=1
        ,axis.text.x = element_blank()
        ,axis.text.y = element_blank()
        ,axis.ticks = element_blank()
        ,axis.title.x = element_blank()
        ,axis.title.y = element_blank()
        ,plot.title = element_text(size = 18, face = "bold",colour="#044c73")
  ) + # end theme
  ggtitle(paste("Mean : ",month(sstDate[p],label=TRUE,abbr=FALSE),sep="")) +
  annotate("text", x=max(aData[,"x"])/10,y=0, 
           label=paste(format(round(pGWL, 1), nsmall = 1),"\u00B0C",sep=""),size = 7,color = pGWLcolour) +
  annotate("text", x=max(aData[,"x"])/10,y=max(aData[,"y"])/7,
           label="GWL",size = 4.5,color = pGWLcolour)

names(pCircle)<-names(aData)
aPlot<-aPlot+geom_polygon(data=pCircle,colour="white",fill=NA)

if(printPlot) print(aPlot)
if(savePlot) ggsave(paste(outputDir,"/",paste("Mean_",month(sstDate[p],label=TRUE,abbr=TRUE),"_",pGWL,sep=""),".png",sep=""))
} # end gwl vector

