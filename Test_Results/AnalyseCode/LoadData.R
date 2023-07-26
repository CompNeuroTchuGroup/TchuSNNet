GetFiles <- function(test_no){
  searchPattern = paste("Test",toString(test_no),"*",sep="");
  all_files = list.files(pattern=glob2rx(searchPattern));
  #filename = list.files(path = folder_name,pattern=glob2rx(searchPattern));
  #full_filename <- paste(folder_name,"/",filename,sep="");
  #df <- data.frame("data" = full_filename);
  
  return(all_files);
}

GetFile <- function(all_files,pattern){
  full_filename = "";
  for (f in all_files) {
    if(grepl(pattern, f, fixed = TRUE))
      full_filename = f;
  }
  return(full_filename);
}

Get_Number_of_Populations <- function(all_files){
  
  full_filename = GetFile(all_files,"Data.dat");
  
  #searchPattern = paste("Test",toString(params$Test),"*_Data.dat",sep="");
  #filename = list.files(path = folder_name,pattern=glob2rx(searchPattern));
  #full_filename <- paste(folder_name,"/",filename,sep="");
  
  data  <- read.table(full_filename, header = TRUE, sep = "\t")
  data <- data[, colSums(is.na(data)) != nrow(data)]
  
  number_populations <- -3 + sqrt(8+ncol(data));
  print(number_populations);
  if(number_populations%%1!=0){
    stop("Number of columns in data file not following convention");
  }
  return(number_populations);
}

Load_data <- function(all_files,number_populations) {
  
  full_filename = GetFile(all_files,"Data.dat");
  
  #filename = list.files(path = folder_name,pattern="*_Data.dat");
  #full_filename <- paste(folder_name,"/",filename,sep="");
  
  data  <- read.table(full_filename, header = TRUE, sep = "\t")
  data <- data[, colSums(is.na(data)) != nrow(data)]
  
  var_names = c("V","r","mu_ext","mu","quenched","temporal");
  col_names = c("time");
  for (var_name in var_names){
    if(var_name == "mu"){
      for(i in 0:(number_populations-1)){
        for(j in 0:(number_populations-1)){
          s_ = paste(var_name,"_",toString(i),"_",toString(j),sep ="")
          col_names = append(col_names,s_);
        }
      }
    }
    for(i in 0:(number_populations-1)){
      s_ = paste(var_name,"_",toString(i),sep ="")
      col_names = append(col_names,s_);
    }
    
  }
  colnames(data) <- col_names;
  
  return(data);
}


PostProcess_data <- function(data,number_populations){
  
  var_names = c("V","r","mu_ext","mu","quenched","temporal");
  df = data.frame();
  
  for(var_name_ in var_names){
    
    cols = c("time");
    for(i in 0:(number_populations-1)){
      cols <- append(cols,paste(var_name_,"_",toString(i),sep=""))
    }
    
    sub_data <- data %>% select(cols) %>% gather(- time, key = "population", value = "value")
    if(ncol(df) == 0){
      df <- data.frame(sub_data);
    }else{
      df[var_name_] = sub_data["value"];
    }
    names(df)[names(df)=="value"] <- var_name_
  }

  for(i in 0:(number_populations-1)){
    df$population[df$population == paste(var_names[1],"_",toString(i),sep="")] <- i
  }
  df$population <- as.factor(df$population);
  
  return(df);
}


Load_Synapses <- function(all_files,number_populations) {
  
  full_filename = GetFile(all_files,"Synapses.dat");
  #filename = list.files(path = folder_name,pattern="*_Synapses.dat");
  #full_filename <- paste(folder_name,"/",filename,sep="");

  col_names = c("time","J","x","y","pR");
  
  if(file.exists(full_filename) && !dir.exists(full_filename)  ){
    data      <- read.table(full_filename, header = TRUE, sep = "\t", na.strings="-nan(ind)")
    data <- data[, colSums(is.na(data)) != nrow(data)]
    
    #Only get the first five columns: (1-time), (2-J),(3-x) , (4-y) , (5-pR)
    data <- data[,1:5];
    colnames(data) <- col_names;
    data <- na.omit(data);
  }else{
    data <- read.table(text="",col.names = col_names)
  }
  
  return(data);
}

Load_Potential <- function(all_files,number_populations) {
  
  full_filename = GetFile(all_files,"Potential.dat");
  
#  filename = list.files(path = folder_name,pattern="*_Potential.dat");
#  full_filename <- paste(folder_name,"/",filename,sep="");
  
  col_names = c("time");
  for(i in 0:(number_populations-1)){
    for(j in 0:(number_populations-1)){
      s_ = paste("V_",toString(i),"_",toString(j),sep ="");
      col_names = append(col_names,s_);
    }
  }
  
  if(file.exists(full_filename) && !dir.exists(full_filename)  ){
    data      <- read.table(full_filename, header = TRUE, sep = "\t")
    data <- data[, colSums(is.na(data)) != nrow(data)]
    colnames(data) <- col_names;
  }else{
    data <- read.table(text="",col.names = col_names)
  }
  
  return(data);
}


Load_Currents <- function(all_files,number_populations) {
  
  full_filename = GetFile(all_files,"Currents.dat");
  
  #filename = list.files(path = folder_name,pattern="*_Currents.dat");
  #full_filename <- paste(folder_name,"/",filename,sep="");
  
  col_names = c("time");
  for(i in 0:(number_populations-1)){
    for(j in 0:(number_populations-1)){
      s_ = paste("mu_",toString(i),"_",toString(j),sep ="");
      col_names = append(col_names,s_);
    }
  }
  
  if(file.exists(full_filename) && !dir.exists(full_filename)  ){
    data      <- read.table(full_filename, header = TRUE, sep = "\t")
    data <- data[, colSums(is.na(data)) != nrow(data)]
    colnames(data) <- col_names;
  }else{
    data <- read.table(text="",col.names = col_names)
  }
  return(data);
}

Load_Rasterplot <- function(all_files,number_populations) {
  
  full_filename = GetFile(all_files,"Rasterplot.dat");
  
  #filename = list.files(path = folder_name,pattern="*_Rasterplot.dat");
  #full_filename <- paste(folder_name,"/",filename,sep="");
  
  col_names = c("time","neuron_id","population_id");
  
  if(file.exists(full_filename) && !dir.exists(full_filename)  ){
    data_rasterplot <- read.table(full_filename, header = TRUE, sep = "\t")
    colnames(data_rasterplot) <- col_names;
  }else{
    data_rasterplot <- read.table(text="",col.names = col_names)
  }
  return(data_rasterplot);
}