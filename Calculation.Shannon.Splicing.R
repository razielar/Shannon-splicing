#################################### Calculation of Shannon splicing from a GTF file ####################################

### Function 1: 
fColnames.GTEx <- function(matrix.name){
  
  require(dplyr)
  
  new_names <- strsplit(colnames(matrix.name), split = ".", fixed = TRUE)
  new_names <- lapply(new_names, function(x){paste0(x,collapse = "-")})
  new_names <- new_names %>% unlist
  colnames(matrix.name) <- new_names
  
  return(matrix.name)
  
}

### Function 2: 
fConserve.more.than.one.isoform <- function(Input.matrix){
  
  require(dplyr) 
  
  # Order alphabetically the Transcript within the dataframe:
  Input.matrix <- Input.matrix[order(rownames(Input.matrix)),]
  
  # Select the isoforms by "-" and classify them:
  new_name <- strsplit(rownames(Input.matrix), split = "-", fixed = TRUE)
  new_name <- lapply(new_name, function(x){ y <-x[1:(length(x)-1 )]; paste0(y,collapse = "-")  } )
  new_name <- new_name %>% unlist()
  new_name <- data.frame(Transcript=new_name)
  new_name <- new_name %>% group_by(Transcript) %>% summarise(Frequency=n()) %>% arrange(desc(Frequency))
  new_name <- new_name[new_name$Frequency == 1,]
  
  # Change briefly the rownames of the matrix Transcript.cutoff from "( & )" to "_" in order to avoid problems using
  #grepl
  tmp.rownames <- gsub("\\(", "_", rownames(Input.matrix)) 
  tmp.rownames <- gsub("\\)", "_", tmp.rownames )
  rownames(Input.matrix) <- tmp.rownames
  
  tmp.new_names <- gsub("\\(", "_", new_name$Transcript)  
  tmp.new_names <- gsub("\\)", "_", tmp.new_names)
  new_name$Transcript <- tmp.new_names
  
  # Start the cycle: 
  for(i in 1:nrow(new_name)){
    tmp<- grepl(paste0("^",new_name$Transcript[i],"-") , rownames(Input.matrix)) %>% table()
    cat(i,"\n")
    if(tmp[2] == 1){
      Input.matrix <- Input.matrix[!grepl(paste0("^",new_name$Transcript[i],"-") , rownames(Input.matrix)),]
    }
  }
  
  # Restore the matrix Transcript.cutoff rownames from "( & )" to "_"
  Tmp.rownames <- gsub("_", "\\(", rownames(Input.matrix))
  Tmp.rownames <- gsub("_", "\\)", Tmp.rownames)
  rownames(Input.matrix) <- Tmp.rownames
  
  return(Input.matrix)
  
}

### Function 3: 
fQuality.control.isoform.ratio <- function(Input.matrix){
  
  require(dplyr)
  
  new_name <- strsplit(rownames(Input.matrix), split = "-", fixed = TRUE)
  new_name <- lapply(new_name, function(x){ y <-x[1:(length(x)-1)]; paste0(y,collapse = "-")})
  new_name <- new_name %>% unlist()
  
  q.cutoff <- data.frame(Transcript=new_name)
  q.cutoff <- q.cutoff %>% group_by(Transcript) %>% summarise(Frequency=n()) %>% arrange(desc(Frequency))
  
  tmp <- q.cutoff[q.cutoff$Frequency == 1,][,1] %>% as.data.frame()
  tmp <- as.vector(tmp[,1])
  
  if(length(tmp) >= 1){
    
    cat("There's one or more genes with only ONE ISOFORM, remove it before continue the analysis:", "\n")
    x <- q.cutoff[q.cutoff$Frequency == 1,][,1] %>% as.data.frame() 
    x <<- as.vector(x[,1])
    #cat(x,"\n")
    
  } else { cat("You can continue with the next function")  } 
  
  
}

### Function 4: 

f.Calculate_Shannon_splicing <- function(Input.matrix){
  
  new_name <- strsplit(rownames(Input.matrix), split = "-", fixed = TRUE)
  new_name <- lapply(new_name, function(x){ y <-x[1:(length(x)-1)]; paste0(y,collapse = "-")}) 
  new_name <- new_name %>% unlist()
  
  tmp.matrix <- Input.matrix %>% mutate(Transcript=new_name)
  
  Shannon_splicing <- function(x){
    
    #1) Calculate the isoform ratio:
    
    glob <- x %>% select(-Transcript) #dataframe: with the replicates
    sumvals <- colSums(glob) #numeric value: sum of the columns
    result <- sapply(1:length(sumvals), function(x){ round(glob[,x]/sumvals[x], digits = 4)}) %>% 
      do.call(rbind, .) %>% t #sapply: return a vector; dataframe
    
    #2) Calculate the Shannon's entropy:
    
    shannon_entropy <- round(-colSums(result*log(result)), digits = 4) #Shannon's formula
    shannon_entropy <- t(as.data.frame(shannon_entropy))
    shannon_entropy <- data.frame(shannon_entropy)
    
  }
  
  tmp.isoform.ratio <- tmp.matrix %>% group_by(Transcript) %>% do(Shannon_splicing(.)) 
  tmp.isoform.ratio <- tmp.isoform.ratio %>% as.data.frame()
  tmp.isoform.ratio <- fColnames.GTEx(tmp.isoform.ratio)
  
  return(tmp.isoform.ratio)
  
  
}


###################################################################################################################

