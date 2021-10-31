## Import dataset

unzip(FB018B_Surface.zip,exdir=getwd())
library(rgdal) 
FB018A<-readOGR("FB018A_Surface.shp")
FB018B<-readOGR("FB018B_Surface.shp")
RR_FB018A<-readOGR("RR_FB018A.shp")
RR_FB018B<-readOGR("RR_FB018B.shp")

library(ggplot2)
library(gridExtra)

## Assess size and preservation proprtion of debris

df_FB018A<-data.frame(RR_FB018A$Dimension)
df_FB018B<-data.frame(RR_FB018B$Dimension)
colnames(df_FB018A)<-"Size";colnames(df_FB018B)<-"Size"

tiff("Hist_Size_Debris.tif",width=7000,height=3500,units="px",res=800)
histA<-ggplotGrob(ggplot(df_FB018A,aes(x=Size))+
                    geom_histogram(binwidth=4,color="black",fill="grey")+
                    labs(title="FB018A",x="debris size"))
histB<-ggplotGrob(ggplot(df_FB018B,aes(x=Size))+
  geom_histogram(binwidth=4,color="black",fill="grey")+
  labs(title="FB018B",x="debris size"))
grid.arrange(histA,histB,ncol=2)
dev.off()

tiff("Hist_Log-Size_Debris.tif",width=7000,height=3500,units="px",res=800)
histlA<-ggplotGrob(ggplot(log(df_FB018A),aes(x=Size))+
                     geom_histogram(binwidth=0.3,color="black",fill="grey")+
                     geom_vline(aes(xintercept=mean(Size)),color="red",
                                linetype="dashed",size=0.8)+
                     labs(title="FB018A",x="log debris size"))
histlB<-ggplotGrob(ggplot(log(df_FB018B),aes(x=Size))+
                     geom_histogram(binwidth=0.3,color="black",fill="grey")+
                     geom_vline(aes(xintercept=mean(Size)),color="red",
                                linetype="dashed",size=0.8)+
                     labs(title="FB018B",x="log debris size"))
grid.arrange(histlA,histlB,ncol=2)
dev.off()

## Plot geographical location of objects + their size

coord_FB018A<-subset(coordinates(RR_FB018A),select=-c(coords.x3))
coord_FB018B<-subset(coordinates(RR_FB018B),select=-c(coords.x3))
# for some reason (Autocad?) my shp contained a third coordinate...

pA<-fortify(FB018A)
siteA_df<-data.frame(pA[,1:2]);colnames(siteA_df)<-c("x","y")
pointsA_df<-data.frame(coord_FB018A,df_FB018A)
colnames(pointsA_df)<-c("x","y","Size")
pB<-fortify(FB018B)
siteB_df<-data.frame(pB[,1:2]);colnames(siteB_df)<-c("x","y")
pointsB_df<-data.frame(coord_FB018B,df_FB018B)
colnames(pointsB_df)<-c("x","y","Size")

tiff("Huts-Plots.tif",width=7000,height=3500,units="px",res=800)
hutA<-ggplotGrob(ggplot(siteA_df,aes(x,y))+
  geom_polygon(color="black",fill="grey")+
  geom_point(data=pointsA_df,shape=21,colour="#483D8B",aes(size=Size))+
    labs(title="FB018A")+
    scale_size_continuous(breaks=c(1,10,100))+
    coord_fixed(ratio=1))
hutB<-ggplotGrob(ggplot(siteB_df,aes(x,y))+
                   geom_polygon(color="black",fill="grey")+
                   geom_point(data=pointsB_df,shape=21,colour="#483D8B",
                              aes(size=Size))+
                   labs(title="FB018B")+
                   scale_size_continuous(breaks=c(3,30,300))+
                   coord_fixed(ratio=1))
grid.arrange(hutA,hutB,ncol=2)
dev.off()
  
### SPATIAL DEPENDENCY

library(spdep)
library(pgirmess)

# Moran's I correlogram (distance classes)

corr_FB018A<-correlog(coord_FB018A,RR_FB018A$Dimension,method="Moran")
corr_FB018B<-correlog(coord_FB018B,RR_FB018B$Dimension,method="Moran")
pvalA<-(corr_FB018A[,3]<0.05)+0; pvalB<-(corr_FB018B[,3]<0.05)+0

tiff("Correlogram.tif",width=7000,height=2500,units="px",res=800)
corrA<-as.data.frame(cbind(corr_FB018A,pvalA))
corrA<-ggplotGrob(ggplot(corrA,aes(x=dist.class,y=coef))+
  geom_point(data=corrA,aes(colour=p.value),size=4)+
    theme(legend.position = "none")+
    coord_fixed(ratio=200)+
    expand_limits(x=c(0,550),y=c(-1,1.25))+
    labs(title="FB018A")+
  geom_line(size=0.4))
corrB<-as.data.frame(cbind(corr_FB018B,pvalB))
corrB<-ggplotGrob(ggplot(corrB,aes(x=dist.class,y=coef))+
                    geom_point(data=corrB,aes(colour=p.value),size=4)+
                    coord_fixed(ratio=200)+
                    labs(title="FB018B")+
                    expand_limits(x=c(0,550),y=c(-1,1.25))+
                    geom_line(size=0.4))
grid.arrange(corrA,corrB,ncol=2)
dev.off()


# Local Moran's I test of RR

nbA<-dnearneigh(RR_FB018A,0,150,longlat=F)
RR_FB018A_lMoran<-localmoran(log(RR_FB018A$Dimension),nb2listw(nbA,style="W"))
nbB<-dnearneigh(RR_FB018B,0,100,longlat=F)
RR_FB018B_lMoran<-localmoran(log(RR_FB018B$Dimension),nb2listw(nbB,style="W"))

# Plot Local Moran

IiA<-RR_FB018A_lMoran[,1]
PrA<-RR_FB018A_lMoran[,5]
pointsA_dfM<-data.frame(coord_FB018A,IiA,PrA)
colnames(pointsA_dfM)<-c("x","y","Ii","Pr")

IiB<-RR_FB018B_lMoran[,1]
PrB<-RR_FB018B_lMoran[,5]
pointsB_dfM<-data.frame(coord_FB018B,IiB,PrB)
colnames(pointsB_dfM)<-c("x","y","Ii","Pr")

tiff("Huts-Moran.tif",width=7000,height=3500,units="px",res=800)
MoranA<-ggplotGrob(ggplot(siteA_df,aes(x,y))+
                   geom_polygon(color="black",fill="grey")+
                   geom_point(data=pointsA_dfM,aes(size=Pr,color=Ii))+
                   labs(title="FB018A")+
                     scale_size_continuous(breaks=c(0.2,0.5,0.8))+
                     coord_fixed(ratio=1))
MoranB<-ggplotGrob(ggplot(siteB_df,aes(x,y))+
                     geom_polygon(color="black",fill="grey")+
                     geom_point(data=pointsB_dfM,aes(size=Pr,color=Ii))+
                     labs(title="FB018B")+
                     scale_size_continuous(breaks=c(0.001,0.05,0.1,0.5,0.8))+
                     coord_fixed(ratio=1))
grid.arrange(MoranA,MoranB,ncol=2)
dev.off()


## SPATIALLY INDEPENDENT RR

valA<-(RR_FB018A_lMoran[,1]<0)+0
points_valA<-data.frame(coord_FB018A,valA,df_FB018A)
colnames(points_valA)<-c("x","y","spat.indep.","Size")
plotIiA<-subset(points_valA,spat.indep.>0)

valB<-(RR_FB018B_lMoran[,1]<0)+0
points_valB<-data.frame(coord_FB018B,valB,df_FB018B)
colnames(points_valB)<-c("x","y","spat.indep.","Size")
plotIiB<-subset(points_valB,spat.indep.>0)

tiff("Huts-Spat-Indep.tif",width=7000,height=3500,units="px",res=800)
SpatIndepA<-ggplotGrob(ggplot(siteA_df,aes(x,y))+
                     geom_polygon(color="black",fill="grey")+
                     geom_point(data=plotIiA,aes(size=Size),shape=21,
                                colour="#483D8B")+
                       scale_size_continuous(breaks=c(1,10,100))+
                       labs(title="FB018A")+
                     coord_fixed(ratio=1))
SpatIndepB<-ggplotGrob(ggplot(siteB_df,aes(x,y))+
                     geom_polygon(color="black",fill="grey")+
                     geom_point(data=plotIiB,aes(size=Size),shape=21,
                                colour="#483D8B")+
                       scale_size_continuous(breaks=c(3,30,300))+
                       labs(title="FB018B")+
                     coord_fixed(ratio=1))
grid.arrange(SpatIndepA,SpatIndepB,ncol=2)
dev.off()
