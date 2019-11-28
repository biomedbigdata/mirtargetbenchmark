#' Functionality to calculate the ensemble scores for all the tools
#' @export
#' @import rlist
#' @param tool_list A list containing the data frames containing scores for all the tools
#' @param tool_names A list containing the names of all the tools in the dataframes
#' @param performance_list A list containing performance scores of all the tools

ensemble_scores <- function(tool_list, tool_names, performance_list){

                    vector_tools <- list()
                    #convert the dataframes in the list into vectors and store in vector_tools list
                    for (x in 1:length(tool_list)) {
                      vector_tools <- list.append(vector_tools, as.vector(unlist(tool_list[[tool_names[x]]])))
                    }
                    names(vector_tools) <- tool_names

                    #scale the scores


                    #calculate the cumulative scores
                    scores = numeric(length = length(vector_tools[[tool_names[1]]]))
                    for(x in 1:length(tool_list)){
                    scores = scores + performance_list[[tool_names[x]]]*vector_tools[[tool_names[x]]]
                    }

                    return(scores)

}
