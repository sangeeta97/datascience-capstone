library(ClustOfVar)
new<- read.csv("C:/Users/sangeeta kumari/Desktop/for poster/Book44.csv", header = TRUE, row.names = "Label")
tree <- hclustvar(new)

> data(decathlon)
> tree <- hclustvar(decathlon[,1:10])
> plot(tree)
> part <- cutreevar(tree,5)
> print(part)

Call:
  cutreevar(obj = tree, k = 4)
read.csv(row.names = Label)

capture.output(summary(mylist), file = "My New File.txt")

part$cluster



name       description                        
"$var"     "list of variables in each cluster"
"$sim"     "similarity matrix in each cluster"
"$cluster" "cluster memberships"              
"$wss"     "within-cluster sum of squares"    
"$E"       "gain in cohesion (in %)"          
"$size"    "size of each cluster"             
"$scores"  "score of each cluster"            
> summary(part)

Call:
  cutreevar(obj = tree, k = 4)



Data: 
  number of observations:  41
number of variables:  10
number of clusters:  4

Cluster  1 : 
  squared loading
100m                   0.68
Long.jump              0.69
400m                   0.67
110m.hurdle            0.64


Cluster  2 : 
  squared loading
Shot.put             0.77
High.jump            0.55
Discus               0.67


Cluster  3 : 
  squared loading
Pole.vault            0.62
1500m                 0.62


Cluster  4 : 
  squared loading
Javeline               1


Gain in cohesion (in %):  54.14
> View(decathlon)
> stability(tree,B=40) 

Call:
  stability(tree = tree, B = 40)


name     description                            
"$matCR" "matrix of corrected Rand indices"     
"$meanR" "vector of mean corrected Rand indices"
> 