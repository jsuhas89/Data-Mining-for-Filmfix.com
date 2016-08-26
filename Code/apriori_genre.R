data1<-read.csv("data1.csv")
data1<-data1[2:9]
for(i in 1:8){
  data1[,i]<-as.factor(data1[,i])
}
rules<-apriori(data1,parameter = list(support = 0.7, confidence = 1))
inspect(rules)