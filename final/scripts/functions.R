library(stringr)

capitalize <- function(word) {
  word <- str_to_lower(word)
  word <- paste(str_to_upper(substring(word, 1, 1)), substring(word, 2), sep="")
  return(word)
}
 
get_taxonomies_to_keep <- function(x) {
  my_subset <- data[data$Time == x[1] & data$horse == x[2] & data$Treatment == x[3],]
  total <- sum(my_subset$Value)
  min_amount <- total * (as.numeric(args[5]) / 100)
  my_subset <- my_subset[my_subset$Value >= min_amount,]
  my_subset[,args[4]]
}