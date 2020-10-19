setwd("~/GitHub/untargeted_processing")
add <- read.csv("data/adducts.csv")
add <- add[add$mzdiff != 0,]
library(Rdisop)
library(CompoundDb)

load("data/RData/FT_NEG.RData")
features$polarity <- "NEG"
rownames(features) <- paste0("N_", rownames(features))
features$FG <- paste0("N_", features$FG)
features$mzneut <- features$mzmed + 1.007276
features_neg <- features

load("data/RData/FT_POS.RData")
features$polarity <- "POS"
rownames(features) <- paste0("P_", rownames(features))
features$FG <- paste0("P_", features$FG)
features$mzneut <- features$mzmed - 1.007276

features <- rbind(features_neg, features)
rm(features_neg)
features <- features[order(features$FG), ]

features$ann <- NA
features$annotation <- NA
features$FGx <- NA
FG <- unique(features$FG)


for(k in 1:length(FG)
    ){
  ft <- features[features$FG == FG[k], ]
  if(sum(is.na(ft$annotation)) == nrow(ft)){
    # deduce the relationsheeps between mzvalues. Parent ion maybe is the one with more associations
    for(i in 1:nrow(ft)){
      for(j in 1:nrow(ft)){
        mzrange <- (ft$mzmed[j] - ft$mzmed[i]) + 0.01 * c(-1, 1)
        idx <- which(add$mzdiff > mzrange[1] & add$mzdiff < mzrange[2])
        if(any(idx)){
          if(is.na(ft$ann[i])){
            ft$ann[i] <- paste(round(ft$mzmed[j], 4), add$assignation[idx])
          } else {
            ft$ann[i] <- paste(
              ft$ann[i], paste(round(ft$mzmed[j], 4), add$assignation[idx]), 
              sep = "; ")
          }
        } else { # is it a dimer?
          mzrange <- (ft$mzneut[i]*2) + 0.01 * c(-1, 1)
          if(ft$mzneut[j] > mzrange[1] & ft$mzneut[j] < mzrange[2]){
            if(is.na(ft$ann[i])){
              ft$ann[i] <- paste(round(ft$mzmed[j], 4), "[2M*H]*")
            } else {
              ft$ann[i] <- paste(
                ft$ann[i], paste(round(ft$mzmed[j], 4), "[2M*H]*"), 
                sep = "; ")
            }
          }}
         # is it an unassigned MS2 fragment?
          tmp <- unlist(strsplit(ft$MS2[i], "; "))
          tmp <- tmp[grep("frag", tmp)]
          for(l in 1:length(tmp)){
            idx <- unlist(matchWithPpm(as.numeric(gsub("\\[.*", "", tmp[l])), ft$mzmed, ppm = 10))
            if(any(idx)){
              if(is.na(ft$ann[i])){
                ft$ann[i] <- paste(round(ft$mzmed[idx], 4), "[frag]*")
              } else {
                ft$ann[i] <- paste(ft$ann[i], 
                                   paste(round(ft$mzmed[idx], 4), "[frag]*"), 
                                   sep = "; ") 
              } 
          }
        }
        
      } # close "j"
    } # close "i"
    
    for(i in 1:nrow(ft)){
      ft$ann[i] <- paste(unique(unlist(strsplit(ft$ann[i], "; "))), collapse = "; ")
    }
    
    ft$FGx <- FG[k]
    i <- which.max(unlist(lapply(strsplit(ft$ann, "; "), length)))
    ft$annotation[i] <- "[M*H]*"
    for(j in 1:nrow(ft)){
      mzrange <- (ft$mzmed[j] - ft$mzmed[i]) + 0.01 * c(-1, 1)
      idx <- which(add$mzdiff > mzrange[1] & add$mzdiff < mzrange[2])
      if(any(idx)){
        ft$annotation[j] <- add$assignation[idx]
      } else { # is it a dimer?
        mzrange <- (ft$mzneut[i]*2) + 0.01 * c(-1, 1)
        if(ft$mzneut[j] > mzrange[1] & ft$mzneut[j] < mzrange[2]){
          ft$annotation[j] <- "[2M*H]*"
        }
      } # close "is it a dimer?"
    }
    
    # is there any other feature in the other polarity with the same rt-mzneut?
    ft_opp <- features[features$polarity != unique(features$polarity[features$FG == FG[k]]), ]
    rtrange <- mean(ft$rtmed) + 10 * c(-1, 1)
    ft_opp <- ft_opp[ft_opp$rtmed > rtrange[1] & ft_opp$rtmed < rtrange[2], ]
    for(i in 1:nrow(ft_opp)){
      mzrange <- (ft_opp$mzneut[i]) + 0.01 * c(-1, 1)
      idx <- which(ft$mzneut > mzrange[1] & ft$mzneut < mzrange[2])
      if(any(idx)){
        ft_opp$ann[i] <- ft$FG[idx]
        ft_opp$annotation[i] <- ft$annotation[idx]
        ft_opp$FGx[i] <- FG[k]
      }
    }
    ft_opp <- ft_opp[!is.na(ft_opp$annotation), ]
    
    ft <- rbind(ft, ft_opp)
    for(i in 1:nrow(ft)){
      features$ann[rownames(features)==rownames(ft)[i]] <- ft$ann[i]
      features$annotation[rownames(features)==rownames(ft)[i]] <- ft$annotation[i]
      features$FGx[rownames(features)==rownames(ft)[i]] <- ft$FGx[i]
    }
    rm(i, j, l, ft, ft_opp, mzrange, rtrange, idx)
  }
}

for(i in 1:length(features)){
  features$ann[i] <- unique(unlist(strsplit(features$ann[i], "; ")))
}

rm(k, FG, add)
features$annotation[features$polarity == "NEG"] <- gsub("\\*", "-", features$annotation[features$polarity == "NEG"])
features$annotation[features$polarity == "POS"] <- gsub("\\*", "+", features$annotation[features$polarity == "POS"])

# match [M*H]* with internal database
cmps <- read.csv("data/database.csv")
cmps$mzneut <- NA
for(i in 1:nrow(cmps)){
  cmps$mzneut[i] <- getMolecule(as.character(cmps$formula[i]))$exactmass
}
features$cmp <- NA
for(i in 1:nrow(features)){
  tmp <- cmps[unlist(matchWithPpm(features$mzneut[i], cmps$mzneut, ppm = 10)), ]
  rtrange <- features$rtmed[i] + 10 * c(-1, 1)
  tmp <- tmp[tmp$RT > rtrange[1] & tmp$RT < rtrange[2], ]
  features$cmp[i] <- paste(tmp$cmp, collapse = "; ")
}
rm(i, tmp, rtrange, cmps)

features <- features[order(features$FGx), ]

# check why there are missing "FGs" from the previous code
# annotate 13C
# check MS2 relationships splitted in different FG



unique(features$FGx)
sum(features$FGx == "N_FG008")
unique(features$FG[features$FGx == "N_FG008"])

sum(features$FG == "N_FG008")
sum(features$FG == "P_FG007")
