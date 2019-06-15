require(automap)
require(ggplot2)
require(RColorBrewer)

#####Geo: The geographic coordinates of each location 
#####value: The observed value of each location
#####df_grid: The prediction grid of the interpolation regions

df<-cbind(Geo,z=value)
data<-SpatialPointsDataFrame(Geo,df)
kriging_result = autoKrige(z~1,data, df_grid)
prediction_spdf = kriging_result$krige_output
a<-coordinates(prediction_spdf)
df<-as.data.frame(cbind(a,value=prediction_spdf$var1.pred))
color<-brewer.pal(11, "RdYlBu")[11:1]
ggplot()+
   geom_polygon(data=china_map1,aes(x=long,y=lat,group=group)
                ,fill="white",colour="black")+
   geom_point(data=df,aes(x1,x2,color=value))+
   scale_color_gradientn(colours = color) +
   coord_map()+
   theme_bw()+
   theme(panel.grid=element_blank())

######Automatic cross-validation
df<-cbind(Geo,z=value)
data<-SpatialPointsDataFrame(Geo,df)
set.seed(1)
aa<-autoKrige.cv(z~1,data,nfold=10)
aa<-as.data.frame(aa$krige.cv_output)
cor.test(aa[,1],aa[,3],method="pear")

