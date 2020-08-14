################################################################################
forceDfOrder <- function(df){
   as.data.frame(lapply(df, function(i){
      if(!is.numeric(i)){ i <- factor(i, unique(i)) }
      return(i)
   }))
}

################################################################################
dfToFeaturesAndResponse <- function(df, colname.response='response'){
   x <- df[,colnames(df)!=colname.response]
   x <- as.matrix(x)

   y <- as.factor(df[,colname.response])
   y_ohe <- oneHotEncode(y, sample.names=rownames(x))

   list(x=x, y=y, y_ohe=y_ohe)
}

################################################################################
univarFeatSel <- function(
   x, y,
   max.qvalue=0.01, sel.top.n.features=NULL, return.new.x=T,
   verbose=F
){
   # colname.response='response'
   # x <- df[,colnames(df)!=colname.response]
   # y <- as.factor(df[,colname.response])
   # #y <- y=='Prostate'

   if( !(is.logical(y) | is.factor(y)) ){
      stop('`y` must be a logical or factor')
   }
   #if(is.logical(y)){ y <- factor(y,c('TRUE','FALSE')) }

   main <- function(v, y.logical){
      #y.logical=y
      if(is.numeric(v)){
         #v=x$viral_ins.Hepatitis_C_virus
         v_split <- split(v, y.logical)
         wilcox.test(v_split[['TRUE']], v_split[['FALSE']])$p.value
      } else {
         #v=ifelse(x$purple.gender,'male','female')
         #v=unname(m[,'AR'])
         #y=metadata[ match(rownames(m), metadata$sample),'cancer_type' ]
         #y.logical <- y=='Prostate'

         neg_category <- levels(as.factor(v))[1]
         #v=rep('0;none',length(v))

         fisher.test(
            matrix(
               c(
                  sum(v!=neg_category & y.logical), sum(y.logical),
                  sum(v!=neg_category & !y.logical), length(y.logical)
               ),
               nrow=2
            )
         )$p.value
      }
   }

   if(is.logical(y)){
      if(verbose){ counter <- 0 }
      p_values <- apply(x, 2, function(i){
         if(verbose){
            counter <<- counter + 1
            message('[',counter,'] ', colnames(x)[[counter]] )
         }
         main(i, y)
      })
      p_values <- sort(p_values)

      q_values <- p.adjust(p_values, method='bonferroni')

      keep_features <- names(q_values)[ q_values < max.qvalue ]
      if(!is.null(sel.top.n.features)){ keep_features <- keep_features[1:sel.top.n.features] }
      keep_features <- na.exclude(keep_features)

   } else {
      y_logicals <- lapply(levels(y), function(i){ y==i })
      names(y_logicals) <- levels(y)

      if(verbose){ counter <- 0 }
      m_p_values <- do.call(cbind, lapply(y_logicals, function(y_logical){
         if(verbose){
            counter <<- counter + 1
            message('[',counter,'] ', names(y_logicals)[[counter]] )
         }
         apply(x, 2, function(feature){ main(feature, y_logical) })
      }))
      m_q_values <- apply(m_p_values, 2, p.adjust, method='bonferroni')

      keep_features <- unlist(apply(m_q_values, 2, function(i){
         names(i)[ i<max.qvalue ]
      }), use.names=F)

      keep_features <- unique(na.exclude(keep_features))
   }

   if(return.new.x){
      return(x[,keep_features])
   }
   return(keep_features)
}



################################################################################
predict.randomForestEnsemble <- function(object, newdata, type='response'){
   # object=model
   # newdata=test_data$x
   # type='response'

   m <- do.call(cbind, lapply(object, function(i){
      randomForest:::predict.randomForest(i, newdata, type='prob')[,1]
   }))

   if(type=='prob'){
      m
   } else {
      factor(
         colnames(m)[ max.col(m) ],
         levels=colnames(m)
      )
   }
}

getPredictions <- function(model, x, output.classes=F){
   ## Neural network
   if(any(grepl('keras',class(model)))){

      class_names <- sapply(0:(reticulate::py_len(model$class_names)-1),function(i){
         model$class_names[[i]]
      })

      prob <- model$predict(unname(x))
      colnames(prob) <- class_names
      rownames(prob) <- if(is.list(x) & !is.data.frame(x)){
         rownames(x[[1]])
      } else {
         rownames(x)
      }

      if(!output.classes){ return(prob) }
      #return( factor(class_names[apply(prob,1,which.max)], class_names) )
      return( factor(class_names[max.col(prob)], class_names) )
   }

   ## Random forest
   if('randomForest' %in% class(model)){
      if(output.classes){
         return(randomForest:::predict.randomForest(model, x, type='response'))
      } else {
         return(randomForest:::predict.randomForest(model, x, type='prob'))
      }
   }

   ## Random forest ensemble
   if('randomForestEnsemble' %in% class(model)){
      if(output.classes){
         return(predict.randomForestEnsemble(model, x, type='response'))
      } else {
         return(predict.randomForestEnsemble(model, x, type='prob'))
      }
   }

   stop('Model must be a randomForest or Keras model')
}

################################################################################
permutationImportance <- function(model, x, y, metric='prec', seed=NULL, verbose=T){
   #x <- train_data$x
   #y <- train_data$y

   #if(!is.matrix(x)){ stop('`x` must be a matrix') }
   if(!is.factor(y)){ stop('`y` must be a factor') }
   if(length(metric)>1){ stop('Only a single metric is allowed') }

   if(!is.null(seed)){ set.seed(seed) }


   if(verbose){ message('Predicting base classes...') }
   base_classes <- getPredictions(model, x, output.classes=T)

   if(verbose){ message('Permuting features...') }
   if(is.list(x) & !is.data.frame(x)){
      x_as_matrix <- do.call(cbind,unname(x))
   } else {
      x_as_matrix <- x
   }

   x_permuted <- lapply(1:ncol(x_as_matrix),function(i){
      #i=1
      m <- x_as_matrix
      m[,i] <- sample(m[,i])
      return(m)
   })

   if(is.list(x) & !is.data.frame(x)){
      feature_groups <- sapply(x,colnames)
      x_permuted <- lapply(x_permuted, function(i){
         lapply(feature_groups,function(j){ i[,j] })
      })
   }

   if(verbose){ message('Predicting permuted classes...') }
   perm_classes <- as.data.frame(lapply(x_permuted,function(i){
      getPredictions(model, i, output.classes=T)
   }), col.names=colnames(x_as_matrix), check.names=F)


   if(verbose){ message('Calculating decrease in performance...') }
   base_metric <- calcPerf(
      confusionMatrix(oneHotEncode(base_classes), y, simplify=T),
      metric
   )[[metric]]

   perf_decrease <- lapply(1:ncol(perm_classes),function(i){
      pert_metric <- calcPerf(
         confusionMatrix(oneHotEncode(perm_classes[,i]), y, simplify=T),
         metric
      )[[metric]]

      base_metric - pert_metric
   })
   perf_decrease <- do.call(cbind, perf_decrease)
   dimnames(perf_decrease) <- list(levels(y), colnames(x_as_matrix))

   return(perf_decrease)
}











