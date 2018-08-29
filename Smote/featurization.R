#' Illustration of Featurization
#'
#' This function takes dataset and a list of features as input and produce a features-wise dataset. The number of columns in returned dataset is equal to the number of features in featurelist.
#'
#' @param sequences provided as dataframe
#' @param string a list of aminoacids or nucleotides
#' @param seq sequence based features. by default it is true.
#' @param seqorder highest number of sequence which will be considered together
#' @param pos position specific features. by default it is true.
#' @param posorder highest number of sequence which will be considered together
#' @return a featurized dataframe
#' @export
#' @examples
#nput = list("AA", "ABCBDBEBEBBBDBDBFDFDFGGHHEEFFEECCCD")
inp= read.csv()
alphabet=c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "Y")
string = c("A","C","T","G")
featuredataAll = featurization(inp$X30mer, string, seq = TRUE,seqorder = 4, pos = TRUE, posorder = 4)

featurization <-
  function(sequences, string, seq = TRUE, seqorder = 2, pos = TRUE, posorder = 2) {
    total = 0;
    alphaMap = new.env();
    for (key in alphabet) {
      assign(key, TRUE, alphaMap);
    }
    
    nGramCount = 0;
    psfCount = 0;
    nGDipCount = 0;
    nGDipOrder=30;
    features = data.frame(1:length(sequences))
    colnames(features)[length(features)] = "Serial"
    if (seq == TRUE) {
      for (s in 1:seqorder) {
        permu = gtools::permutations(
          n = length(string), r = s, v = string, repeats.allowed = TRUE
        )
        for (i in 1:length(permu[,1])) {
          temp = countpattern(sequence = sequences, pattern = paste(permu[i,], collapse = ''))
          #cat(length(temp),permu[i,],"\n")
          if (sum(temp) >= 0) {
            features = data.frame(features, temp)
            colnames(features)[length(features)] = paste(permu[i,], collapse = '')
          }
          total = total + 1
        }
        cat(s, " order seq. features:", length(permu[,1]), ":total features = ", total, "\n")
      }
    }
    if (pos == TRUE)
    {
      minlength = min(unlist(lapply(sequences, function(s) {
        nchar(toString(s))
      })))
      for (p in 1:posorder) {
        permu = gtools::permutations(
          n = length(string), r = p, v = string, repeats.allowed = TRUE
        )
        ps = 0
        for (i in 1:(length(permu[,1]))) {
          for (j in 1:(minlength - length(permu[i,]) + 1)) {
            #cat("Checking:",paste(permu[i,],collapse = ''),"\n")
            temp = findposition(sequence = sequences, pattern = paste(permu[i,], collapse = ''), j)
            #cat(length(temp),permu[i,],"\n")
            if (sum(temp) >= 0) {
              features = data.frame(features, temp)
              colnames(features)[length(features)] = paste0(paste(permu[i,], collapse = ''), "_", j)
            }
            ps = ps + 1
            total = total + 1
          }
        }
        cat(p, "order pos. features:", ps,":total features = ", total, "\n")
      }
    }
    
    
    if (nGDipOrder > 0) {
      cat('yes')
      for (i in 1:nrow(features)) {
        strSeq = strsplit(toString(sequences[i]), "")[[1]];
        for (j in 1:length(strSeq)) {
          if (!(exists(strSeq[j], envir = alphaMap))) {
            next;
          }
          
          for (k in 1:nGDipOrder) {
            if (j+1+k > length(strSeq)) {
              break;
            }
            if (!(exists(strSeq[j+1+k], envir = alphaMap))) {
              next;
            }
            token = paste(strSeq[j], strSeq[j+1+k], sep = "");
            token = paste("G", k, token, sep = "_");
            if (!(token %in% colnames(features))) {
              # create the column on demand
              features[token] = integer(nrow(features));
              nGDipCount = nGDipCount + 1;
            }
            
            features[i,token] = features[i,token] + 1/(length(strSeq) - k -1);
          }
        }
      }
    }
    
    features$Serial = NULL
    #viennaData = viennaRNADataManipulation(sequences)
    #features = data.frame(features, viennaData)
    #cat("2 viennaRNAdata related features :total features = ", (total + 2), "\n")
    return(features)
  }


write.csv(featuredataAll,"featuredata_smote.csv")