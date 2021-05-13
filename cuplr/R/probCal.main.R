#' Isotonic regression
#'
#' @description A wrapper for stats::isoreg()
#'
#' @param x x values
#' @param y y values
#'
#' @return A list which the objects:
#' - orig_data: Original x,y values ordered by the x values
#' - curve: The x,y values of the curve
#' @export
#'
isoReg <- function(x,y){

   coords <- data.frame(x,y)
   coords <- coords[order(coords$x),]

   ## Fit isotonic regression
   fit <- isoreg(coords$x, coords$y)
   fit$x <- sort(fit$x)
   fit$y <- sort(fit$y)

   orig_data <- data.frame(x=fit$x, y=fit$y)

   curve <- data.frame(
      x = fit$x[ fit$iKnots ],
      y = fit$yf[ fit$iKnots ]
   )

   ## Add y-values for x-values (i.e. original probabilities) of 0 and 1
   curve <- rbind(
      within(curve[1,], x <- 0),
      curve,
      within(curve[nrow(curve),], x <- 1)
   )

   ## Dedup rows with same y-values
   curve <- unique(curve)

   curve$is_uniq_first <- !duplicated(curve$y)
   curve$is_uniq_last  <- !duplicated(curve$y, fromLast=T)
   curve$is_dup_bounding <- curve$is_uniq_first | curve$is_uniq_last
   curve <- curve[curve$is_dup_bounding,]

   NULL -> curve$is_uniq_first -> curve$is_uniq_last -> curve$is_dup_bounding

   ## Output
   list(orig_data=orig_data, curve=curve)
}

####################################################################################################
#' Generate probability calibration curves
#'
#' @rdname probCal
#'
#' @param actual A vector of the actual classes
#' @param probs A matrix where rows are samples, cols are binary random forest names, and cells are
#' the prediction probabilities from each random forest
#' @param report A list with the objects with the names: prob, class_actual
#' @param output If 'coords' returns a dataframe with the x,y values of the calibration curve for
#' each `actual` class. If 'plot', returns a ggplot object
#' @param curve.coords The output from `probCalCurves(..., output='coords')`
#'
#' @description `probCalCurves()` generates calibration curves to map raw classifier
#' (pseudo-)probabilibities to real probabilities. This is done using isotonic regression, with
#' x-values being raw probs and y-values being the sample labels (1 or 0). `probCal()` is used to
#' map raw probabilities to calibrated probabilities
#'
#' @return `probCalCurves()`: A dataframe or ggplot object (see `output`). `probCal()` a matrix
#' of calibrated probabilities per class
#'
#' @export
#'
probCalCurves <- function(actual=NULL, probs=NULL, report=NULL, output=c('coords','plot')){

   # if(F){
   #    report=pred_reports$CV
   # }

   ## Init --------------------------------
   if(!is.null(report)){
      actual <- report$class_actual
      probs <- report$prob
   }

   output <- match.arg(output, c('coords','plot'))

   ## Main --------------------------------
   uniq_classes <- colnames(probs)
   isoreg_fits <- lapply(uniq_classes, function(i){
      #i='Cervix'

      ## Prep data
      df <- data.frame(
         class_target = i,
         response = as.integer(actual==i),
         class_target_prob = probs[,i],
         row.names=NULL
      )
      df <- df[order(df$class_target_prob, decreasing=F),]

      out <- isoReg(x=df$class_target_prob, y=df$response)
      lapply(out, function(j){
         j$class <- i
         return(j)
      })
   })

   orig_data <- do.call(rbind, lapply(isoreg_fits,`[[`,'orig_data'))
   curve_coords <- do.call(rbind, lapply(isoreg_fits,`[[`,'curve'))
   rm(isoreg_fits)

   ## Output --------------------------------
   if(output=='coords'){ return(curve_coords) }

   ggplot(curve_coords, aes(x=x, y=y)) +

      facet_wrap(class~.) +

      geom_abline(slope=1, intercept=0, linetype='dashed', color='grey') +

      ## Original data
      geom_point(data=orig_data, shape=21, color='red') +

      ## Curve from fit
      geom_line() +

      scale_y_continuous(name='Calibrated prob.\n(red: sample label)') +
      scale_x_continuous(name='Raw prob. (from classifier)') +
      theme_bw() +
      theme(
         panel.grid.minor=element_blank()
      )
}

#' @rdname probCal
#' @export
#'
probCal <- function(probs, curve.coords){
   # if(F){
   #    probs=pred_reports$CV$prob
   #    curve.coords=calib_curves
   # }

   ## Init --------------------------------
   if(!all(colnames(probs) %in% unique(curve.coords$class))){
      stop('`colnames(probs)` must have the same classes as in `curve.coords`')
   }

   ## Main --------------------------------
   probs_scaled <- lapply(colnames(probs), function(i){
      #i='Breast'
      lookup_table <- curve.coords[curve.coords$class==i,]
      approx(x=lookup_table$x, y=lookup_table$y, probs[,i])$y
   })

   probs_scaled <- do.call(cbind, probs_scaled)
   dimnames(probs_scaled) <- dimnames(probs)
   return(probs_scaled)
}
