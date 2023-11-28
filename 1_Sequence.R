samples <- read.csv("data/phenodata.csv")
samples <- samples[,c("id", "color", "variety")]
samples$code <- paste(gsub("white", "wht", samples$color),
                      substr(
                        gsub("pinot", "pt", gsub(" ", "", samples$variety)
                             ), 1, 6),
                      "A", sep = "_")
for(i in unique(samples$code[duplicated(samples$code)])){
  idx <- which(samples$code == i)
  samples$code[idx] <- paste0(substr(samples$code[idx], 1, 11), 
                              LETTERS[1:length(idx)])
}

samples$vial <- paste0("R", rep(LETTERS[1:5], each = 8), 
                       rep(seq(8), 5))[3:(nrow(samples) + 2)]

eq <- data.frame(matrix(ncol = ncol(samples), nrow = 7))
qc <- data.frame(matrix(ncol = ncol(samples), nrow = 3))
colnames(eq) <- colnames(qc) <- colnames(samples)
eq$id <- eq$random <- eq$block <- qc$id <- qc$random <- 0
eq$color <- eq$variety <- qc$color <- qc$variety <- "x"
eq$code <- paste("X", c(rep("eqblnk", 2), rep("eqQC", 5)), "X", sep = "_")
qc$code <- paste("X", c("blank", rep("QC", 2)), "X", sep = "_")
eq$vial <- c(rep("RA1", 2), rep("RA2", 5))
qc$vial <- c("RA1", rep("RA2", 2))
eq$rep <- c(seq(2), seq(5))

set.seed(2022106)
r <- sample(seq(1000, 9999), 4*2)
n <- 0

for(p in c("POS", "NEG")){
  for(m in c("FS", "DDA", "HRp", "HRc")){
    n <- n + 1
    if(m == "FS"){
      smpl <- do.call("rbind", replicate(3, samples, simplify = FALSE))
    } else {
      smpl <- samples
    }
    set.seed(r[n])
    smpl$random <- sample(seq(nrow(smpl)), nrow(smpl), replace = F)
    smpl <- smpl[order(smpl$id, smpl$random),]
    if(m == "FS"){
      smpl$rep <- rep(seq(3), nrow(samples))
    } else {
      smpl$rep <- 1
    }
    smpl <- smpl[order(smpl$random),]
    smpl$block <- rep(seq(nrow(smpl)/10), each = 10)
    t_qc <- do.call("rbind", replicate(nrow(smpl)/10 + 1, qc, simplify = FALSE))
    t_qc <- t_qc[order(t_qc$code),]
    t_qc$rep <- c(seq(nrow(smpl)/10 + 1), seq((nrow(smpl)/10 + 1)*2))
    t_qc <- t_qc[order(as.numeric(rownames(t_qc))),]
    t_qc$block <- rep(seq(nrow(smpl)/10 + 1), each = 3)
    if(m == "FS"){
      smpl <- rbind(eq, t_qc, smpl)
    } else {
      smpl <- rbind(t_qc, smpl)
      smpl <- smpl[smpl$variety != "blank",]
    }
    smpl <- smpl[order(smpl$block, smpl$random),]
    smpl$pol <- p
    smpl$mode <- m
    smpl$order <- paste0("x", formatC(seq(nrow(smpl)), 
                                      width = 3, flag = "0"))
    
    if(exists("sq")){
      sq <- rbind(sq, smpl)
    } else {
      sq <- smpl
    }
  }
}
sq$name <- paste0(sq$order, "_", sq$code, "_", sq$rep, "_", sq$pol, "_", 
                  sq$mode)
write.csv(sq, "output/sequence_20221006.csv", row.names = FALSE)
