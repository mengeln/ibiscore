NorCal_IBI <- function(locationinfo, data, DistinctCode=F, Grid=F, SampleDate=F, FieldReplicate=F){
  starttime <- proc.time()
  load("ibi.RData")
  ibi <- idata.frame(ibi)
  data <- IBIname_match(data=data, DistinctCode=DistinctCode)
  colnames(data)[which(colnames(data) == "FinalID")] <- "Taxa"
  colnames(data)[which(colnames(data) == "BAResult")] <- "Result"
  ###Calculate total count###
  total_count <- tapply(data$Result, data$SampleID, sum)
  ###Create sample count flag###
  sample_count_flag <- rep(NA, length(total_count))
  names(sample_count_flag) <- names(total_count)
  sample_count_flag[(which(total_count>=500))] <- "Adequate"
  sample_count_flag[(which(total_count < 500 & total_count >= 450))] <- "Within specifications"
  sample_count_flag[(which(total_count < 450))] <- "Inadequate"
  ###Subsample down to 500###
  datalength <- length(data)
  print("Starting 20 iterations of rarification")
  library(doParallel)
  registerDoParallel()
  rarificationresult <- foreach(i=1:20, .combine=cbind, .export="rarify") %dopar% {
    rarify(inbug=data, sample.ID= "SampleID", abund ="Result", 500)$Result
  }
  print("Rarification complete")
  data <- cbind(data, rarificationresult)
  colnames(data)[(datalength + 1):(datalength + 20)]<- paste("Replicate", 1:20)
  ###Number of Coleoptera taxa###
  metrics <- as.data.frame(matrix(NA, nrow = length(unique(data$SampleID)), ncol = 140))
  for(i in 1:20){
    temp <- tapply(data$SAFIT[data$distinct == "Distinct" & data[[datalength + i]]>0], data$SampleID[data$distinct == "Distinct" & data[[datalength + i]]>0],
                   function(d)sum(unique(d) %in% ibi$SAFIT[ibi$Order == "Coleoptera"]))
    metrics[[i]] <- temp[!is.na(temp)]
  }
  ###Numer of EPT taxa
  for(i in 1:20){
    temp <- tapply(data$SAFIT[data$distinct == "Distinct" & data[[datalength + i]]>0], data$SampleID[data$distinct == "Distinct" & data[[datalength + i]]>0],
                   function(d)sum(unique(d) %in% ibi$SAFIT[ibi$Order %in% c("Ephemoptera", "Plecoptera", "Trichoptera")]))
    metrics[[i+20]] <- temp[!is.na(temp)]
  }
  ###Number of Diptera taxa###
  for(i in 1:20){
    temp <- tapply(data$SAFIT[data$distinct == "Distinct" & data[[datalength + i]]>0], data$SampleID[data$distinct == "Distinct" & data[[datalength + i]]>0],
                   function(d)sum(unique(d) %in% ibi$SAFIT[ibi$Order == "Diptera"]))
    metrics[[i+40]] <- temp[!is.na(temp)]
  }
  ###Percent predators###
  for(i in 1:20){
    predatortotal <- tapply(data$SAFIT[data$distinct == "Distinct" & data[[datalength + i]]>0], data$SampleID[data$distinct == "Distinct" & data[[datalength + i]]>0], function(d)length(unique(d)))
    predatorcount <- tapply(data$SAFIT[data$distinct == "Distinct" & data[[datalength + i]]>0], data$SampleID[data$distinct == "Distinct" & data[[datalength + i]]>0],
                            function(d)sum(unique(d) %in% ibi$SAFIT[ibi$FunctionalFeedingGroup == "P"]))
    temp <- round(((predatortotal-predatorcount)/predatortotal)*100)
    metrics[[i+60]] <- temp[!is.na(temp)]
  }
  ###Percent Non-Insect taxa###
  for(i in 1:20){
    taxatotal <- tapply(data$SAFIT[data$distinct == "Distinct" & data[[datalength + i]]>0], data$SampleID[data$distinct == "Distinct" & data[[datalength + i]]>0], function(d)length(unique(d)))
    insect <-  tapply(data$SAFIT[data$distinct == "Distinct" & data[[datalength + i]]>0], data$SampleID[data$distinct == "Distinct" & data[[datalength + i]]>0],
                      function(d)sum(unique(d) %in% ibi$SAFIT[ibi$Class == "Insecta"]))
    temp <- round(((taxatotal-insect)/taxatotal)*100)
    metrics[[i+80]] <- temp[!is.na(temp)]
  }
  ###Percent Intolerant###
  data$tolerance <- ibi$MaxOfToleranceValue[match(data$Taxa, ibi$FinalID)]
  for(i in 1:20){
    toleranttotal2 <- tapply(data[data[[datalength + i]]>0 & (!is.na(data$tolerance)), datalength + i], as.character(data$SampleID[data[[datalength + i]]>0 & (!is.na(data$tolerance))]), sum)
    tol2 <- tapply(data[data[[datalength + i]]>0 , datalength + i], as.character(data$SampleID[data[[datalength + i]]>0]), 
                   function(d)sum(d[which(d <= 2)], na.rm=T))
    temp <- round(100*(tol2/toleranttotal2))
    temp <- temp[which(names(temp) != "")]
    metrics[[i+100]] <- temp
  }
  for(i in 101:120){
    metrics[which(is.na(metrics[[i]])), i] <- 0
  }
  ###Percent non-gastropod scraper###
  data$feed <- ibi$FunctionalFeedingGroup[match(data$Taxa, ibi$FinalID)]
  data$class <- ibi$Class[match(data$Taxa, ibi$FinalID)]
  for(i in 1:20){
    scrapertotal <- tapply(data[data$distinct == "Distinct" & data[[datalength + i]]>0 & (data$feed != "NONE REPORTED"), datalength + i], data$SampleID[data$distinct == "Distinct" & data[[datalength + i]]>0 & (data$feed != "NONE REPORTED")], sum)
    scraper <- tapply(data[data[[datalength + i]]>0 & data$class != "Gastropoda" & data$feed == "SC", datalength + i], data$SampleID[data[[datalength + i]]>0 & data$class != "Gastropoda" & data$feed == "SC"], 
                      function(d)sum(d, na.rm=T))
    scraper[is.na(scraper)] <- 0
    temp <- round(100*(scraper/scrapertotal))
    metrics[[i+120]] <- temp[!is.na(temp)]
  }
  ###Percent shredder###
  for(i in 1:20){
    shredtotal <- tapply(data[data[[datalength + i]]>0 & (data$feed != "NONE REPORTED"), datalength + i], data$SampleID[data[[datalength + i]]>0 & (data$feed != "NONE REPORTED")], sum)
    shred <- tapply(data[data[[datalength + i]]>0 & data$feed == "SH", datalength + i], data$SampleID[data[[datalength + i]]>0 & data$feed == "SH"], 
                    function(d)sum(d, na.rm=T))
    shred[is.na(shred)] <- 0
    temp <- round(100*(shred/shredtotal))
    metrics[[i+140]] <- temp[!is.na(temp)]
  }
  ###Identify Ecoregion###
  Ecoregion <- IBIlocation_N(locationinfo)
  data$Ecoregion <- Ecoregion[match(data$StationCode, names(Ecoregion))]
  stationlocation <- unique(data[, c("SampleID", "Ecoregion")])
  coast <- which(stationlocation[[2]]=="Coast Range")
  klamath <- which(stationlocation[[2]]=="Klamath") 
  chap <- which(stationlocation[[2]]=="Chaparral")
  missing <- which(is.na(stationlocation[[2]]))
  ###Print map###
  pdf(file = "Ecoregion_map_NorCal.pdf")
  print(IBImap(locationinfo, data, Ecoregion, type="Ecoregion", zone = "NorCal"))
  dev.off()
  ###Convert metrics to scores###
  scores <- metrics
  ###Coleoptera scores###
  for(i in 1:20){
    scores[intersect(chap, which(scores[[i]]>=8)), i] <- 10
    scores[intersect(chap, which(scores[[i]]==7)), i] <- 9
    scores[intersect(chap, which(scores[[i]]==6)), i] <- 8
    scores[intersect(chap, which(scores[[i]]==5)), i] <- 6
    scores[intersect(chap, which(scores[[i]]==4)), i] <- 5
    scores[intersect(chap, which(scores[[i]]==3)), i] <- 4
    
    scores[intersect(c(coast, klamath), which(scores[[i]]>=6)), i] <- 10
    scores[intersect(c(coast, klamath), which(scores[[i]]==5)), i] <- 9
    scores[intersect(c(coast, klamath), which(scores[[i]]==4)), i] <- 7
    scores[intersect(c(coast, klamath), which(scores[[i]]==3)), i] <- 5
    scores[intersect(c(coast, klamath), which(scores[[i]]==2)), i] <- 1
    scores[missing, i] <- NA
  }
  ###EPT scores###
  options(warn = -1)
  for(i in 1:20){
    scores[intersect(chap, which(scores[[i+20]]<=2)), i+20] <- 0
    scores[intersect(chap, which(scores[[i+20]]>=3 & scores[[i+20]]<=4)), i+20] <- 1
    scores[intersect(chap, which(scores[[i+20]]>=5 & scores[[i+20]]<=6)), i+20] <- 2
    scores[intersect(chap, which(scores[[i+20]]>=7 & scores[[i+20]]<=8)), i+20] <- 3
    scores[intersect(chap, which(scores[[i+20]]>=9 & scores[[i+20]]<=10)), i+20] <- 4
    scores[intersect(chap, which(scores[[i+20]]>=11 & scores[[i+20]]<=12)), i+20] <- 5
    scores[intersect(chap, which(scores[[i+20]]>=13 & scores[[i+20]]<=14)), i+20] <- 6
    scores[intersect(chap, which(scores[[i+20]]>=15 & scores[[i+20]]<=16)), i+20] <- 7
    scores[intersect(chap, which(scores[[i+20]]>=17 & scores[[i+20]]<=18)), i+20] <- 8
    scores[intersect(chap, which(scores[[i+20]]>=19 & scores[[i+20]]<=20)), i+20] <- 9
    scores[intersect(chap, which(scores[[i+20]]>=21)), i+20] <- 10
    scores[missing, i+20] <- NA
  }
  
  ###Diptera###
  for(i in 1:20){
    scores[scores[[i+40]]>=10, i+40] <- 10 
    scores[missing, i+40] <- NA
  }
  ###Predator###
  options(warn = -1)
  for(i in 1:20){
    scores[intersect(c(chap, coast), which(scores[[i+60]]<=1)), i+60] <- 0
    scores[intersect(c(chap, coast), which(scores[[i+60]]==2)), i+60] <- 1
    scores[intersect(c(chap, coast), which(scores[[i+60]]>=3 & scores[[i+60]]<=4)), i+60] <- 2
    scores[intersect(c(chap, coast), which(scores[[i+60]]==5)), i+60] <- 3
    scores[intersect(c(chap, coast), which(scores[[i+60]]>=6 & scores[[i+60]]<=7)), i+60] <- 4
    scores[intersect(c(chap, coast), which(scores[[i+60]]==8)), i+60] <- 5
    scores[intersect(c(chap, coast), which(scores[[i+60]]>=9 & scores[[i+60]]<=10)), i+60] <- 6
    scores[intersect(c(chap, coast), which(scores[[i+60]]==11)), i+60] <- 7
    scores[intersect(c(chap, coast), which(scores[[i+60]]>=12 & scores[[i+60]]<=13)), i+60] <- 8
    scores[intersect(c(chap, coast), which(scores[[i+60]]>=14 & scores[[i+60]]<=15)), i+60] <- 9
    scores[intersect(c(chap, coast), which(scores[[i+60]]>=16)), i+60] <- 10
    
    scores[intersect(klamath, which(scores[[i+60]]<=1)), i+60] <- 0
    scores[intersect(klamath, which(scores[[i+60]]>=2 & scores[[i+60]]<=3)), i+60] <- 1
    scores[intersect(klamath, which(scores[[i+60]]>=4 & scores[[i+60]]<=5)), i+60] <- 2
    scores[intersect(klamath, which(scores[[i+60]]>=6 & scores[[i+60]]<=7)), i+60] <- 3
    scores[intersect(klamath, which(scores[[i+60]]>=8 & scores[[i+60]]<=9)), i+60] <- 4
    scores[intersect(klamath, which(scores[[i+60]]>=10 & scores[[i+60]]<=12)), i+60] <- 5
    scores[intersect(klamath, which(scores[[i+60]]>=13 & scores[[i+60]]<=14)), i+60] <- 6
    scores[intersect(klamath, which(scores[[i+60]]==15 & scores[[i+60]]<=16)), i+60] <- 7
    scores[intersect(klamath, which(scores[[i+60]]>=17 & scores[[i+60]]<=18)), i+60] <- 8
    scores[intersect(klamath, which(scores[[i+60]]>=19 & scores[[i+60]]<=21)), i+60] <- 9
    scores[intersect(klamath, which(scores[[i+60]]>=22)), i+60] <- 10
    
    scores[missing, i+60] <- NA
  }
  ###Non-insect###
  for(i in 1:20){
    scores[scores[[i+80]]>=8 & scores[[i+80]]<=13, i+80] <- 9
    scores[scores[[i+80]]<=7, i+80] <- 10
    scores[scores[[i+80]]>= 57, i+80] <- 0
    scores[scores[[i+80]]>= 52 & scores[[i+80]]<=56, i+80] <- 1
    scores[scores[[i+80]]>= 47 & scores[[i+80]]<=51, i+80] <- 2
    scores[scores[[i+80]]>=41 & scores[[i+80]]<=46, i+80] <- 3
    scores[scores[[i+80]]>=36 & scores[[i+80]]<=40, i+80] <- 4
    scores[scores[[i+80]]>=30 & scores[[i+80]]<=35, i+80] <- 5
    scores[scores[[i+80]]>=25 & scores[[i+80]]<=29, i+80] <- 6
    scores[scores[[i+80]]>=19 & scores[[i+80]]<=24, i+80] <- 7
    scores[scores[[i+80]]>=14 & scores[[i+80]]<=18, i+80] <- 8 
    scores[missing, i+80] <- NA
  }
  ###Percent Intolerant###
  options(warn = -1)
  data$Watershed <- sapply(1:length(data$StationCode), function(i)locationinfo$WSA_lm[locationinfo$StationCode == data$StationCode[i]])
  Watershed <- tapply(data$Watershed, data$SampleID, unique)
  for(i in 1:20){
    scores[, i+100] <- sapply(1:length(Watershed), function(j)(scores[j, i+100] - (-0.089*log10(Watershed[j])+0.433)+0.296))
    scores[intersect(chap, which(scores[[i+100]]< -4)), i+100] <- 0
    scores[intersect(chap, which(scores[[i+100]]>=-4 & scores[[i+100]]<=0)), i+100] <- 1
    scores[intersect(chap, which(scores[[i+100]]>=0 & scores[[i+100]]<=2)), i+100] <- 2
    scores[intersect(chap, which(scores[[i+100]]>=3 & scores[[i+100]]<=6)), i+100] <- 3
    scores[intersect(chap, which(scores[[i+100]]>=7 & scores[[i+100]]<=9)), i+100] <- 4
    scores[intersect(chap, which(scores[[i+100]]>=10 & scores[[i+100]]<=13)), i+100] <- 5
    scores[intersect(chap, which(scores[[i+100]]>=14 & scores[[i+100]]<=16)), i+100] <- 6
    scores[intersect(chap, which(scores[[i+100]]>=17 & scores[[i+100]]<=20)), i+100] <- 7
    scores[intersect(chap, which(scores[[i+100]]>=21 & scores[[i+100]]<=23)), i+100] <- 8
    scores[intersect(chap, which(scores[[i+100]]>=24 & scores[[i+100]]<=27)), i+100] <- 9
    scores[intersect(chap, which(scores[[i+100]]>=28)), i+100] <- 10
    
    scores[intersect(c(coast, klamath), which(scores[[i+100]]< -4)), i+100] <- 0
    scores[intersect(c(coast, klamath), which(scores[[i+100]]>=-4 & scores[[i+100]]<=0)), i+100] <- 1
    scores[intersect(c(coast, klamath), which(scores[[i+100]]>=1 & scores[[i+100]]<=5)), i+100] <- 2
    scores[intersect(c(coast, klamath), which(scores[[i+100]]>=6 & scores[[i+100]]<=10)), i+100] <- 3
    scores[intersect(c(coast, klamath), which(scores[[i+100]]>=11 & scores[[i+100]]<=15)), i+100] <- 4
    scores[intersect(c(coast, klamath), which(scores[[i+100]]>=16 & scores[[i+100]]<=20)), i+100] <- 5
    scores[intersect(c(coast, klamath), which(scores[[i+100]]>=21 & scores[[i+100]]<=25)), i+100] <- 6
    scores[intersect(c(coast, klamath), which(scores[[i+100]]>=26 & scores[[i+100]]<=30)), i+100] <- 7
    scores[intersect(c(coast, klamath), which(scores[[i+100]]>=31 & scores[[i+100]]<=35)), i+100] <- 8
    scores[intersect(c(coast, klamath), which(scores[[i+100]]>=36 & scores[[i+100]]<=40)), i+100] <- 9
    scores[intersect(c(coast, klamath), which(scores[[i+100]]>=41)), i+100] <- 10
    
    scores[missing, i+100] <- NA
  }
  ###Non-gastropod scrapers###
  options(warn = -1)
  for(i in 1:20){
    scores[intersect(c(chap, klamath), which(scores[[i+120]]>=1 & scores[[i+120]]<=2)), i+120] <- 1
    scores[intersect(c(chap, klamath), which(scores[[i+120]]>=3 & scores[[i+120]]<=4)), i+120] <- 2
    scores[intersect(c(chap, klamath), which(scores[[i+120]]>=5 & scores[[i+120]]<=6)), i+120] <- 3
    scores[intersect(c(chap, klamath), which(scores[[i+120]]>=7 & scores[[i+120]]<=8)), i+120] <- 4
    scores[intersect(c(chap, klamath), which(scores[[i+120]]>=9 & scores[[i+120]]<=10)), i+120] <- 5
    scores[intersect(c(chap, klamath), which(scores[[i+120]]>=11 & scores[[i+120]]<=12)), i+120] <- 6
    scores[intersect(c(chap, klamath), which(scores[[i+120]]>=13 & scores[[i+120]]<=14)), i+120] <- 7
    scores[intersect(c(chap, klamath), which(scores[[i+120]]>=15 & scores[[i+120]]<=16)), i+120] <- 8
    scores[intersect(c(chap, klamath), which(scores[[i+120]]==17)), i+120] <- 9
    scores[intersect(c(chap, klamath), which(scores[[i+120]]>=18)), i+120] <- 10
    
    scores[intersect(coast, which(scores[[i+120]]<=4)), i+120] <- 1
    scores[intersect(coast, which(scores[[i+120]]>=5 & scores[[i+120]]<=8)), i+120] <- 1
    scores[intersect(coast, which(scores[[i+120]]>=9 & scores[[i+120]]<=12)), i+120] <- 2
    scores[intersect(coast, which(scores[[i+120]]>=13 & scores[[i+120]]<=16)), i+120] <- 3
    scores[intersect(coast, which(scores[[i+120]]>=17 & scores[[i+120]]<=20)), i+120] <- 4
    scores[intersect(coast, which(scores[[i+120]]>=21 & scores[[i+120]]<=24)), i+120] <- 5
    scores[intersect(coast, which(scores[[i+120]]>=25 & scores[[i+120]]<=28)), i+120] <- 6
    scores[intersect(coast, which(scores[[i+120]]>=29 & scores[[i+120]]<=32)), i+120] <- 7
    scores[intersect(coast, which(scores[[i+120]]>=33 & scores[[i+120]]<=36)), i+120] <- 8
    scores[intersect(coast, which(scores[[i+120]]==37 & scores[[i+120]]<=40)), i+120] <- 9
    scores[intersect(coast, which(scores[[i+120]]>=41)), i+120] <- 10
    
    scores[missing, i+120] <- NA
  }
  ###shredder###
  options(warn = -1)
  for(i in 1:20){
    scores[intersect(c(chap, klamath), which(scores[[i+140]]<= 1)), i+140] <- 0
    scores[intersect(c(chap, klamath), which(scores[[i+140]]==2)), i+140] <- 1
    scores[intersect(c(chap, klamath), which(scores[[i+140]]>=3 & scores[[i+140]]<=4)), i+140] <- 2
    scores[intersect(c(chap, klamath), which(scores[[i+140]]==5)), i+140] <- 3
    scores[intersect(c(chap, klamath), which(scores[[i+140]]>=6 & scores[[i+140]]<=7)), i+140] <- 4
    scores[intersect(c(chap, klamath), which(scores[[i+140]]==8)), i+140] <- 5
    scores[intersect(c(chap, klamath), which(scores[[i+140]]>=9 & scores[[i+140]]<=10)), i+140] <- 6
    scores[intersect(c(chap, klamath), which(scores[[i+140]]==11)), i+140] <- 7
    scores[intersect(c(chap, klamath), which(scores[[i+140]]>=12 & scores[[i+140]]<=13)), i+140] <- 8
    scores[intersect(c(chap, klamath), which(scores[[i+140]]>=14 & scores[[i+140]]<=15)), i+140] <- 9
    scores[intersect(c(chap, klamath), which(scores[[i+140]]>=16)), i+140] <- 10
    
    scores[intersect(coast, which(scores[[i+140]]<= 1)), i+140] <- 0
    scores[intersect(coast, which(scores[[i+140]]>=2 & scores[[i+140]]<=3)), i+140] <- 1
    scores[intersect(coast, which(scores[[i+140]]>=4 & scores[[i+140]]<=5)), i+140] <- 2
    scores[intersect(coast, which(scores[[i+140]]>=6 & scores[[i+140]]<=7)), i+140] <- 3
    scores[intersect(coast, which(scores[[i+140]]>=8 & scores[[i+140]]<=9)), i+140] <- 4
    scores[intersect(coast, which(scores[[i+140]]>=10 & scores[[i+140]]<=11)), i+140] <- 5
    scores[intersect(coast, which(scores[[i+140]]>=12 & scores[[i+140]]<=13)), i+140] <- 6
    scores[intersect(coast, which(scores[[i+140]]>=14 & scores[[i+140]]<=15)), i+140] <- 7
    scores[intersect(coast, which(scores[[i+140]]>=16 & scores[[i+140]]<=17)), i+140] <- 8
    scores[intersect(coast, which(scores[[i+140]]>=18 & scores[[i+140]]<=19)), i+140] <- 9
    scores[intersect(coast, which(scores[[i+140]]>=20)), i+140] <- 10
    
    scores[missing, i+140] <- NA
  }
  ###SCIBI###
  for(i in 0:19){
    scores[[i+161]] <- sapply(1:length(unique(data$SampleID)), function(j)(10/8)*(sum(c(scores[j, 1+i], scores[j, 20+i], 
                                                                                        scores[j, 40+i], scores[j, 60+i], scores[j, 80+i], scores[j, 100+i], scores[j, 120+i], scores[j, 140+i]))))
  }
  ###Calculate means for metrics and scores###
  means <- as.data.frame(matrix(NA, nrow=length(unique(data$SampleID)), ncol = 15))
  for(i in 1:8){
    means[[i]] <- apply(metrics[, (((i-1)*20)+1):(20*i)], 1, mean)
  }  
  for(i in 1:9){
    means[[i+8]] <- apply(scores[, (((i-1)*20)+1):(20*i)], 1, mean)
  }
  ###Construct output frame###
  results <- as.data.frame(matrix(NA, nrow=length(unique(data$SampleID)), ncol = 23))
  results[[1]] <- unique(data[, c("StationCode", "SampleID")])$StationCode
  results[[2]] <- unique(data$SampleID)
  results[[3]] <- unique(data[, c("SampleID", "Ecoregion")])$Ecoregion
  results[[4]] <- total_count[!is.na(total_count)]
  results[[5]] <- sample_count_flag[!is.na(sample_count_flag)]
  results[[6]] <- rep(20, times=length(unique(data$SampleID)))
  results[, 7:23] <- means
  colnames(results) <- c("StationCode", "SampleID", "Ecoregion", "Total Count", "Count Flag", "Number of Iteration", 
                         "Number of Coleoptera Taxa", "Number of EPT Taxa", "Number of Diperta Taxa", "Percent Predators", 
                         "Percent Non Insect Taxa", "Percent Intolerant", "Non gastropod scrapers", "Percent Shredder Individuals", "Coleoptera Score",
                         "EPT Score", "Diperta Score", "Predator Score", "Non Insect Score", "Intolerant Score",
                         "Non gastropod scrapers Score", "Shredder Score", "NCIBI")
  ###Representativeness flag###
  if(Grid==T){
    data$Representativeness_Flag <- rep(NA, length(data$StationCode))
    data[which(is.na(data$TotalGrids)), "TotalGrids"] <- "Missing"
    flag <- sapply(1:length(data$StationCode), function(i)
      if(data$TotalGrids[i] == "Missing"){NA}else 
        if(as.numeric(data$TotalGrids[i]) >= 3 | data$GridsAnalyzed[i] > 2 | (data$GridsVolumeAnalyzed[i]/
          as.numeric(data$TotalGrids[i])) >= .25){"Representative"}else{"Potentially nonrepresentative"})
    data$Representativeness_Flag <- flag[!is.na(flag)]
  }
  ###Optional fields###
  options(warn = -1)
  if(Grid==T){
    results$"Representativeness Flag" <- data$"Representativeness_Flag"[match(results$SampleID, data$SampleID)]
  }
  if(FieldReplicate==T){
    results$FieldReplicate <- data$FieldReplicate[match(results$SampleID, data$SampleID)]
  }
  if(SampleDate==T){
    results$SampleDate <- data$SampleDate[match(results$SampleID, data$SampleID)]
  }
  extrastuff <- sum(c(Grid, FieldReplicate, SampleDate))
  if(extrastuff>0){
    results <- results[, c(1:5, 24:(22+extrastuff), 6:24)]  
  }
  ###Print SCIBI map###
  ncibi <- results$NCIBI
  names(ncibi) <- results$StationCode
  ncibi <- ncibi[match(names(Ecoregion), names(ncibi))]
  names(ncibi) <- names(Ecoregion)
  ncibi <- unlist(ncibi)
  pdf(file="NCIBI_map_NorCal.pdf")
  print(IBImap(locationinfo, data, ncibi, type="NCIBI", zone="NorCal"))
  dev.off()
  ###Return results###
  print(proc.time() - starttime)
  return(results)
}