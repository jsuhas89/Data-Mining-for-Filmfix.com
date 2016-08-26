data1<-read.csv("data_knn.csv")
data1<-data1[2:11]
for(i in 1:10){
  data1[,i]<-as.factor(data1[,i])
}
rules<-apriori(data1,parameter = list(support = 0.3, confidence = 0.2))
inspect(rules)