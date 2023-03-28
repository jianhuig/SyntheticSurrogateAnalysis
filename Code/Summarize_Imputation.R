# read all files with .linear extension
filenames <- list.files(path = "Results/ImputationAnalysis", pattern = ".linear", full.names = T)

# get the names of the files
# only keep everything after plink2
# and before .linear

for(filename in filenames){
    # get the name of the file
    name <- gsub("(?:[^.]+\\.)([^.]+).*", "\\1", filename)
    temp <- read.table(filename, header = T)
    assign(paste0(name), temp)
}

# read SynSurr results
