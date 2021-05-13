#' Combine multiple lists which have the same objects
#'
#' @param lists A list of lists, where each list has the same objects
#' @param exclude.objects A character vector specifying the names of the objects within each list to
#' exclude
#' @param show.warnings Show warnings?
#' @param verbose Show progress messages?
#'
#' @return A list
#' @export
#'
combineLists <- function(lists, exclude.objects=NULL, show.warnings=T, verbose=F){
   #lists=reports

   ## Checks --------------------------------
   ## Object names
   all_object_names <- lapply(lists, function(i){
      sort(names(i))
   })

   all_object_names_order <- unique(unlist(lapply(unname(lists), names)))

   if(length(unique(all_object_names))==1){
      object_names <- all_object_names[[1]]
   } else {
      if(show.warnings){
         warning('Reports do not have the same object names. Only merging common objects.')
      }
      object_names <- Reduce(intersect, all_object_names)
   }

   ## Restore original order of objects
   object_names <- object_names[
      na.exclude(match(all_object_names_order, object_names))
   ]

   ## Exclude objects
   if(!is.null(exclude.objects)){
      object_names <- object_names[!(object_names %in% exclude.objects)]
   }

   ## Use common objects and force order
   lists <- lapply(lists,function(i){
      i[ object_names ]
   })

   ## Object types
   all_object_types <- lapply(lists, function(i){
      sapply(i,class)
   })

   object_types <- all_object_types[[1]]
   if(length(unique(all_object_types))!=1 & show.warnings){
      warning('Some objects across lists do not have the same type.\nUsing object types from 1st list')
   }

   ## Main --------------------------------
   out <- lapply(object_names, function(i){
      #i=object_names[[1]]
      #i='responses_pred'
      i_type <- object_types[[i]]
      if(any(c('data.frame','matrix') %in% i_type)){
         combine_func <- function(l){ do.call(rbind, l) }
      } else if('factor' %in% i_type){
         combine_func <- function(l){
            uniq_lvls <- unique(unlist(lapply(l, levels), use.names=F))
            v <- factor(
               unlist(lapply(l, as.character), use.names=F),
               levels=uniq_lvls
            )
            names(v) <- unlist(lapply(l, names), use.names=F)
            return(v)
         }
      } else {
         combine_func <- function(l){ do.call(c, l) }
      }

      if(verbose){ message('Merging ',i,' (class: ',paste(i_type, collapse=', '),')') }

      i_type_objects <- lapply(1:length(lists), function(j){
         #j=1
         lists[[j]][[i]]
      })

      combine_func(i_type_objects)
   })
   names(out) <- object_names
   return(out)
}

# #' Combine results from multiple prediction reports into one list object
# #'
# #' @param reports A list of reports. Each report (a list) has the one or more of the following
# #' objects with the names: probs_raw, probs_adjusted, responses_pred, responses_actual,
# #' feat_contrib, imp
# #' @param verbose Show progress messages?
# #'
# #' @return A list
# #' @export
# #'
# mergePredReports <- function(reports, verbose=T){
#
#    reports_merged <- list()
#
#    if(verbose){ 'Merging probs_raw...' }
#    reports_merged$probs_raw <- do.call(rbind, lapply(reports,`[[`,'probs_raw'))
#
#    if(verbose){ 'Merging probs_adjusted...' }
#    reports_merged$probs_adjusted <- do.call(rbind, lapply(reports,`[[`,'probs_adjusted'))
#
#    if(verbose){ 'Merging predicted responses...' }
#    reports_merged$responses_pred <- structure(
#       unlist(lapply(reports,`[[`,'responses_pred')),
#       names=unlist(lapply(reports,function(i){ names(i$responses_pred) }))
#    )
#
#    if(verbose){ 'Merging actual responses...' }
#    reports_merged$responses_actual <- structure(
#       unlist(lapply(reports,`[[`,'responses_actual')),
#       names=unlist(lapply(reports,function(i){ names(i$responses_actual) }))
#    )
#
#    if(verbose){ 'Merging feature contributions...' }
#    reports_merged$feat_contrib <- do.call(rbind, lapply(reports,`[[`,'feat_contrib'))
#    reports_merged$feat_contrib <- subset(reports_merged$feat_contrib, contrib>0)
#
#    if(verbose){ 'Merging feature importances...' }
#    reports_merged$imp <- aggregateMatrixList(lapply(reports,`[[`,'imp'), as.matrix=T)
#
#    return(reports_merged)
# }

