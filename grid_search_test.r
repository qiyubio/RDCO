## Author: Qi Yu
## Contact: dryuqi@gmail.com

library(tidyverse)
library(rpart)

d <- mtcars %>% 
  # Convert `am` to factor and select relevant variables
  mutate(am = factor(am, labels = c("Automatic", "Manual"))) %>% 
  select(am, mpg, hp)

ggplot(d, aes(mpg, hp, color = am)) +
  geom_point()

set.seed(245)
n <- nrow(d)
train_rows <- sample(seq(n), size = .8 * n)
train <- d[ train_rows, ]
test  <- d[-train_rows, ]

gs <- list(minsplit = c(2, 5, 10),
           maxdepth = c(1, 3, 8)) %>% 
  cross_df()
gs
mod <- function(...) {
  rpart(am ~ hp + mpg, data = train, control = rpart.control(...))
}

#pmap and pwalk allow you to provide any number of arguments

gs <- gs %>% mutate(fit = pmap(gs, mod))
gs

compute_accuracy <- function(fit, test_features, test_labels) {
  predicted <- predict(fit, test_features, type = "class")
  mean(predicted == test_labels)
}

test_features <- test %>% select(-am)
test_labels   <- test$am

#map_dbl()  return a double vector
gs <- gs %>%
  mutate(test_accuracy = map_dbl(fit, compute_accuracy,
                                 test_features, test_labels))


gs <- gs %>% arrange(desc(test_accuracy), desc(minsplit), maxdepth)


parat<-list( nI=c(10000,15000,20000,25000,30000,35000,40000,45000,50000),
            tel_weight=c(0,1),
            tss_weight=c(0,1),
            h3k4_weight=c(0,1),
            chrs='chr38',sam=1) %>%
  cross_df()
            
source('simulatehotspots_biowulf_onechr.R')
fit = pmap(parat, simulatehotspots_biowulf)



library(doParallel)
registerDoParallel(cores = 24)
library(foreach)


para_Rfile <- file("para.Rfile", "w")
cv.sse <- foreach(nI = seq(10000, 50000, 5000), .combine = rbind) %dopar% {
  foreach(tel_weight = 0:1, .combine = rbind) %dopar% {
  foreach(tss_weight = 0:1, .combine = rbind) %dopar% {
  foreach(h3k4_weight = 0:1, .combine = rbind) %dopar% {
    foutn<-paste(nI,tel_weight,tss_weight,h3k4_weight,"program.R",sep="_")
    fout <- file(foutn, "w")
    
    # TRAIN A PROJECTION PURSUIT REGRESSION WITH VARIOUS SETTINGS AND TRAINING DATA
    #ppreg <- ppr(Y ~ X, data = set1, nterms = n, sm.method = "supsmu", bass = b)
    # CALCULATE SSE WITH VALIDATION DATA
    #test.sse <- sum((set2$Y - predict(ppreg, set2))^2)
    #data.frame(bass = b, nterms = n, sse = test.sse)
    #simulatehotspots_biowulf(nI=nI,tel_weight=tel_weight,tss_weight=tss_weight,h3k4_weight=h3k4_weight,chrs="chr38",sam=1)
    cat(sprintf("source('simulatehotspots_biowulf_wholechr.R')\nsimulatehotspots_biowulf(nI=%s,tel_weight=%s,tss_weight=%s,h3k4_weight=%s,sam=1)\n",nI,tel_weight,tss_weight,h3k4_weight),file=fout)
    cat(sprintf("cd /data/yuq3/simulation/grid_search_whole_genome; module load bedtools macs deeptools R ucsc bwa; R --vanilla < %s >%s.out\n",foutn,foutn),file=para_Rfile)
    
    }
  }
  }
}



