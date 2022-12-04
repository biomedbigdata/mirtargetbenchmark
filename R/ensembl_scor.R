#' @title Calculate scores for all tools
#' @description
#' Functionality to calculate the ensemble scores for all the tools
#' @import rlist
#' @param tool_list A list containing the data frames containing scores for all the tools
#' @param tool_names A list containing the names of all the tools in the dataframes
#' @param performance_list A list containing performance scores of all the tools in decreasing order
#' @export
ensemble_scores <- function(tool_list, tool_names, performance_list){

                    names(tool_list) <- tool_names

                    for (x in 1:length(tool_list)) {
                      tool_list[[tool_names[x]]] <- (tool_list[[tool_names[x]]] - min(tool_list[[tool_names[x]]]))/(max(tool_list[[tool_names[x]]]) - min(tool_list[[tool_names[x]]]))
                    }


                    #scale the scores
                    weights <- list()
                    for(x in 1:length(tool_list)){
                      weights <- list.append(weights, (length(tool_list) - x + 1)/length(tool_list))

                    }


                    names(weights) <- tool_names

                    #calculate the cumulative scores
                    scores = matrix(0,nrow = nrow(tool_list[[tool_names[1]]]),
                                    ncol = ncol(tool_list[[tool_names[1]]]))


                    for(x in 1:length(tool_list)){
                    scores = scores + weights[[tool_names[x]]]*tool_list[[tool_names[x]]]
                    }

                    scores <- (scores - min(scores))/(max(scores) - min(scores))

                    return(scores)

}
