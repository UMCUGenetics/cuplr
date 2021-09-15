#' Isotonic regression
#'
#' @rdname isoReg
#'
#' @description A wrapper for stats::isoreg()
#'
#' @param x Input x values
#' @param y Input y values
#'
#' @return A dataframe containing the x,y values of the curve
#' @export
#'
isoReg <- function(x,y){

   # if(F){
   #    report=pred_reports$CV
   #    sel_class <- 'HeadAndNeck_Other'
   #    x=report$prob[,sel_class]
   #    y=report$class_actual==sel_class
   # }

   coords <- data.frame(x,y)
   coords <- coords[order(coords$x),]

   ## Fit isotonic regression
   fit <- isoreg(coords$x, coords$y)
   fit$x <- sort(fit$x)
   fit$y <- sort(fit$y)

   #orig_data <- data.frame(x=fit$x, y=fit$y)

   curve <- data.frame(
      x = fit$x[ fit$iKnots ],
      y = fit$yf[ fit$iKnots ]
   )

   ## Force curve to start at (0,0) and (1,max(y))
   curve <- rbind(
      within(curve[1,], { x <- 0; y <- 0 }),
      curve,
      within(curve[nrow(curve),], x <- 1)
   )

   curve <- curve[
      !(curve$x==0 & duplicated(curve$x)) &
      !(curve$x==1 & duplicated(curve$x, fromLast=T))
   ,]

   ## Dedup rows with same y-values
   curve <- unique(curve)

   curve$is_uniq_first <- !duplicated(curve$y)
   curve$is_uniq_last  <- !duplicated(curve$y, fromLast=T)
   curve$is_dup_bounding <- curve$is_uniq_first | curve$is_uniq_last
   curve <- curve[curve$is_dup_bounding,]

   NULL -> curve$is_uniq_first -> curve$is_uniq_last -> curve$is_dup_bounding

   ## Output
   #list(orig_data=orig_data, curve=curve)
   class(curve) <- c('isoReg',class(curve))
   return(curve)
}

#' @method predict isoReg
#' @export
predict.isoReg <- function(object, newdata){
   #object=isoReg(x=df$class_target_prob, y=df$response)
   #newdata=seq(0,1,0.01)
   approx(x=object$x, y=object$y, xout=newdata)$y
}

####################################################################################################
#' Logistic regression (1 dimensional)
#'
#' @rdname logReg
#'
#' @description A wrapper for stats::glm()
#'
#' @param x Input x values
#' @param y Input y values
#'
#' @return A numeric vector with the intercept (b0) and gradient (b1)
#' @export
#'
logReg <- function(x, y){
   #x=df$class_target_prob
   #y=df$response
   coords <- data.frame(x,y)
   fit <- glm(data=coords, formula=y~x, family='binomial')
   coefs <- data.frame(
      b0=fit$coefficients[['(Intercept)']],
      b1=fit$coefficients[['x']]
   )
   class(coefs) <- c('logReg',class(coefs))
   return(coefs)
}

#' @method predict logReg
#' @export
predict.logReg <- function(object, newdata, scale01=T){
   #object=log_reg
   #newdata=seq(0,1,0.01)
   require(stats)
   intercept <- object[['b0']]
   gradient <- object[['b1']]

   main <- function(x){
      logit <- intercept+gradient*x
      1/(1+(exp(-logit)))
   }

   out <- main(newdata)

   if(scale01){
      yrange <- main(c(min=0, max=1))
      out <- (out - yrange[['min']]) / (yrange[['max']] - yrange[['min']])
   }

   return(out)
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
#' each `actual` class. If 'plot', returns a ggplot object. If 'plotdata', returns a list with the
#' raw plot data
#' @param method Type of regression to perform. Can be 'isotonic' or 'logistic'.
#' @param prob.prescale Scale probs from 0 to 1 for each sample before fitting or applying the
#' calibration curves
#' @param facet.ncol Number of facet columns when `output`=='plot'
#' @param cal.curves The output from `probCalCurves(..., output='curve')`
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
probCalCurves <- function(
   actual=NULL, probs=NULL, report=NULL,
   output=c('curve','plot','plotdata'), method=c('isotonic','logistic'), prob.prescale=T,
   facet.ncol=NULL
){
   # if(F){
   #    report=pred_reports$CV
   #    actual=report$class_actual
   #    probs=report$prob
   # }

   ## Init --------------------------------
   if(!is.null(report)){
      actual <- report$class_actual
      probs <- report$prob
   }

   output <- match.arg(output, c('curve','plot','plotdata'))
   method <- match.arg(method, c('isotonic','logistic'))

   if(prob.prescale){
      ## Adjusted probs to sum to 1
      ## In theory not necessary but leads to better calibrated prob performance
      probs <- probs / rowSums(probs)
   }

   ## Main --------------------------------
   uniq_classes <- colnames(probs)

   coords <- lapply(uniq_classes, function(i){
      #i='Skin_Other'

      ## Prep data
      df <- data.frame(
         class = i,
         x = probs[,i],
         y = as.integer(actual==i),
         row.names=NULL
      )
      df[order(df$x),]
   })
   names(coords) <- uniq_classes

   reg_func <- switch(
      method,
      isotonic=isoReg,
      logistic=logReg
   )

   fit <- lapply(uniq_classes, function(i){
      out <- reg_func(coords[[i]]$x, coords[[i]]$y)
      out$class <- i
      return(out)
   })
   names(fit) <- uniq_classes

   if(output=='curve'){
      return( do.call(rbind, unname(fit)) )
   }

   ## Plot --------------------------------
   orig_data <- do.call(rbind, unname(coords))
   if(method=='isotonic'){
      curve_coords <- do.call(rbind, fit)
   } else {
      x_values <- seq(0,1,0.01)
      curve_coords <- lapply(names(fit), function(i){
         data.frame(
            x=x_values,
            y=predict(fit[[i]], x_values),
            class=i
         )
      })
      curve_coords <- do.call(rbind, curve_coords)
   }

   if(output=='plotdata'){
      return( list(curve_coords=curve_coords, orig_data=orig_data) )
   }

   ggplot(curve_coords, aes(x=x, y=y)) +

      facet_wrap(class~., ncol=facet.ncol) +

      geom_abline(slope=1, intercept=0, linetype='dashed', color='grey') +

      ## Original data
      geom_point(data=orig_data, shape=21, color='red') +

      ## Curve from fit
      geom_line() +

      scale_y_continuous(name='Calibrated prob.\n(red dot: one sample)') +
      scale_x_continuous(name='Raw prob. (from classifier)') +
      theme_bw() +
      theme(
         panel.grid.minor=element_blank()
      )
}

#' @rdname probCal
#' @export
#'
probCal <- function(probs, cal.curves, prob.prescale=T){
   # if(F){
   #    probs=pred_reports$CV$prob
   #    cal.curves=cal_curves
   # }

   ## Init --------------------------------
   if(!all(colnames(probs) %in% unique(cal.curves$class))){
      stop('`colnames(probs)` must have the same classes as in `cal.curves`')
   }

   ## Main --------------------------------
   if(prob.prescale){
      probs <- probs / rowSums(probs)
   }

   probs_scaled <- lapply(colnames(probs), function(i){
      #i='HeadAndNeck_Other'
      predict(cal.curves[cal.curves$class==i,], probs[,i])
   })

   probs_scaled <- do.call(cbind, probs_scaled)
   dimnames(probs_scaled) <- dimnames(probs)
   return(probs_scaled)
}
