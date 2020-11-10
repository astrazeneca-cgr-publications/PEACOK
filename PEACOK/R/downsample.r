downsample.by.gender <- function(samples, phenotypes, alpha = 0.05, seed = 1234) {
  set.seed(seed)
  # hack for visible global binding
  female.case <- NULL
  male.case <- NULL
  female.control <- NULL
  male.control <- NULL
  genotype <- NULL
  case <- NULL
  control <- NULL
  userID <- NULL
  
  
  phenotypes <- samples %>% inner_join(phenotypes) # the firsr two columns is now userID and is.female
  summary <- data.frame(genotype = character(),  female.case = integer(), male.case = integer(),
                        female.control = integer(), male.control = integer()) 
  for (i in 3:dim(phenotypes)[2]) {
    summary <- summary %>% add_row(genotype = names(phenotypes)[i], 
                                   female.case = length(which(phenotypes$is.female & phenotypes[,i] ==1)),
                                   male.case  = length(which(!phenotypes$is.female & phenotypes[,i] ==1)),
                                   female.control = length(which(phenotypes$is.female & phenotypes[,i] ==0)),
                                   male.control = length(which(!phenotypes$is.female & phenotypes[,i] ==0)))
  }
  
  summary <- summary %>% mutate(case = female.case  + male.case) %>% 
    mutate(control = female.control  + male.control) %>%
    select(genotype, case,control, female.case, female.control, male.case, male.control)
  
  downsample <- phenotypes
  
  for (i in 3:dim(downsample)[2]) {
    current <- summary[i - 2 ,] # shift by 2 for summary to match phenotype columns
    if (current$male.case == 0) { # all female case so set all controls to female
      index <- which(!downsample$is.female & downsample[,i] == 0)
      downsample[index  ,i] <- NA
    } else if (current$female.case == 0) { #all male, so set all controls to male
      index = which(downsample$is.female & downsample[,i] == 0)
      downsample[index  ,i] <- NA
    } else {
      if (current$case > 0 && current$control> 0) {
        x <- matrix(c(current$female.case,current$male.case,current$female.control,current$male.control),nrow = 2)
        fet <- fisher.test(x)
        if (fet$p.value < alpha) { #fet$estimate < 0.5 ||  fet$estimate > 2
          case.ratio <- current$female.case / current$case
          contro.ratio <- current$female.control / current$control
          # there is a need for subsampling
          if (contro.ratio > case.ratio) { #too many females in controls
            contro.size.for.female = round(1.0 * current$male.control * current$female.case / current$male.case)
            temp <- rep(NA,current$female.control)
            temp[sample(current$female.control,contro.size.for.female)] <- 0
            index <- which(downsample$is.female & downsample[,i] == 0)
            downsample[index  ,i] <- temp
          } else {   # too many males in controls
            contro.size.for.male = round(1.0 * current$female.control  * current$male.case / current$female.case)
            temp <- rep(NA,current$male.control)
            temp[sample(current$male.control,contro.size.for.male)] <- 0
            index <- which(!downsample$is.female & downsample[,i] == 0)
            downsample[index  ,i] <- temp
          }
        }
      }
    }
  }
  
  # output
  downsample <- downsample[,-2]
  downsample
}