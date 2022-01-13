# Print pdf versions of html slides
library(httr)
library(pagedown)

# list html files
filelist <- list.files(".", pattern = "^lect.*\\.html", 
                       full.names = TRUE, recursive = TRUE)

# Render html lectures to pdf 
if(length(filelist) > 0){
  for (f in filelist){
    chrome_print(f)
  }
}