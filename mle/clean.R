library(openxlsx)

df<-read.xlsx("./Iowa 2013 Pooling Data.xlsx")

# remove NA gender and filter the female population
df <- df[!is.na(df$Gender),]
df_female <- df[df$Gender=='F',]

# remove NA age and negative age
df_female <- df_female[!is.na(df_female$Age) & df_female$Age>0, ]

# sort by Pool.ID from small to large
#df_female$Pool.ID <- ifelse(is.na(df_female$Pool.ID),999999,df_female$Pool.ID) 
df_female <- df_female[order(df_female$Pool.ID),]

df_female$Specimen.Type <- ifelse(df_female$Specimen.Type=="Cervix" 
                                  | df_female$Specimen.Type=="Vaginal","Swab","Urine")
# race variable
df_female <- df_female[!is.na(df_female$Race),]
df_female$Race <- ifelse(df_female$Race=="W",1,0)

#unique(df_female$Risk.New.Partner)
df_female$Risk.New.Partner <- ifelse(df_female$Risk.New.Partner=="Y",1,0)

#unique(df_female$Risk.Multiple.Partners)
df_female$Risk.Multiple.Partners <- ifelse(df_female$Risk.Multiple.Partners=="Y",1,0)

#unique(df_female$Risk.Contact)
df_female$Risk.Contact <- ifelse(df_female$Risk.Contact=="Y",1,0)

#unique(df_female$Symptoms)
df_female <- df_female[!is.na(df_female$Symptoms),]
df_female$Symptoms <- ifelse(df_female$Symptoms=="Y",1,0)

# test response
df_female <- df_female[!is.na(df_female$CT.Result),]
df_female <- df_female[!(df_female$CT.Result=="E"),]
df_female$CT.Result <- ifelse(df_female$CT.Result=="P",1,0)
df_female$GC.Result <- ifelse(df_female$GC.Result=="P",1,0)

# remove three special individuals: 103715, 107651, 108239 
# since Specimen Type is different within one group
df_female <- df_female[!(df_female$Pool.ID%in%c(103715, 107651, 108239) & !is.na(df_female$Pool.ID)),]

# group size
grp<-as.numeric(table(df_female$Pool.ID))
cj<-c(rep(grp,grp), rep(1,length(df_female$Pool.ID[is.na(df_female$Pool.ID)])))
df_female<-cbind(cj,df_female)

# group response
grp_id<-unique(df_female$Pool.ID)[-length(unique(df_female$Pool.ID))]
GT.Result<-NULL
for (id in grp_id){
  GT.Result <- rbind(GT.Result,apply(df_female[c("CT.Result","GC.Result")][df_female$Pool.ID==id & !is.na(df_female$Pool.ID),],2,max))
}
GT.Result <- cbind(rep(GT.Result[,1],grp), rep(GT.Result[,2],grp));colnames(GT.Result)<-c("CT.Result","GC.Result")
GT.Result <- rbind(GT.Result, df_female[c("CT.Result","GC.Result")][is.na(df_female$Pool.ID),])
colnames(GT.Result)<-c("CT.Group","GC.Group")
df_female <- cbind(GT.Result,df_female)

keeps<-c("Pool.ID","Specimen.Type",
         "cj","CT.Result","GC.Result","CT.Group","GC.Group",
         "Age","Race","Risk.New.Partner","Risk.Multiple.Partners","Risk.Contact","Symptoms")
df_female <- df_female[keeps]

df_female_swab <- df_female[df_female$Specimen.Type=="Swab",]
df_female_urine <- df_female[df_female$Specimen.Type=="Urine",]

#View(df_female)
#View(df_female_swab)
#View(df_female_urine)

write.csv(df_female_swab,file = "./female_swab.csv")
write.csv(df_female_urine,file = "./female_urine.csv")
