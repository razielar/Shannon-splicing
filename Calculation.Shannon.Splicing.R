#################################### Calculation of Shannon splicing from a GTF file ####################################


### Function 1: 
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

### Function 2: 
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

### Function 3: 
fCalculate.isoform.ratio.per.gene <- function(Input.Matrix){
  
  require(dplyr)
  
  new_name <- strsplit(rownames(Input.Matrix), split = "-", fixed = TRUE)
  new_name <- lapply(new_name, function(x){ y <-x[1:(length(x)-1)]; paste0(y,collapse = "-")}) 
  new_name <- new_name %>% unlist()
  
  tmp.matrix <- Input.Matrix %>% mutate(Transcript=new_name)
  
  test <- function(x){
    
    glob <- x %>% select(-Transcript)
    sumvals <- colSums(glob)
    result <- sapply(1:length(sumvals), function(x){ round(glob[,x]/sumvals[x], digits = 4)}) %>% 
      do.call(rbind, .) %>% t
    
    return(as.data.frame(result))
    
  }
  
  tmp.isoform.ratio <- tmp.matrix %>% group_by(Transcript) %>% do(test(.)) 
  
  tmp.isoform.ratio <- tmp.isoform.ratio[,-1] %>% as.data.frame() 
  rownames(tmp.isoform.ratio) <- rownames(Input.Matrix)
  
  return(tmp.isoform.ratio) 
  
}

### Function 4: 
f.Shannon.entropy <- function(Input.matrix){
  
  require(dplyr)
  
  new_name <- strsplit(rownames(Input.matrix), split = "-", fixed = TRUE)
  new_name <- lapply(new_name, function(x){ y <-x[1:(length(x)-1)]; paste0(y,collapse = "-")}) 
  new_name <- new_name %>% unlist()
  
  tmp.matrix <- Input.matrix %>% mutate(Transcript=new_name)
  
  entropy <- function(x){
    
    glob <- x%>% select(-Transcript)
    result <- round( -colSums(glob*log(glob)), digits = 4) #Shannon entropy
    result <- t(as.data.frame(result))
    result <- data.frame(result)
    #Out.put <<- rbind(result)
    
  }
  
  tmp.Shannon <- tmp.matrix %>% group_by(Transcript) %>% do(entropy(.))
  tmp.Shannon <- tmp.Shannon %>% as.data.frame()
  tmp.Shannon <- fColnames.GTEx(tmp.Shannon)
  
}


###################################################################################################################

