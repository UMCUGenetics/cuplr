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
univarFeatSel <- function(x, y, method='median', max.qvalue=0.01, sel.top.n.features=NULL, verbose=T){
   #df=training_data
   #x <- df[,colnames(df)!=colname.response]
   #y <- as.factor(df[,colname.response])

   #print(head(x))
   #print(y)

   if(method=='median'){
      whitelistByPairwiseTest <- function(m_pos, m_neg){
         med_pos <- apply(m_pos,2,median)
         med_neg <- apply(m_neg,2,median)

         med_diff <- med_pos - med_neg
         names(med_diff)[ med_diff > 0 ]
      }

   } else if(method=='wilcox'){
      whitelistByPairwiseTest <- function(m_pos, m_neg){
         #m_pos=x_split[['Prostate']]
         #m_neg=do.call(rbind, unname(x_split[names(x_split)!='Prostate']))

         p_values <- unlist(lapply(1:ncol(m_pos), function(i){
            #i=1
            wilcox.test(m_pos[,i], m_neg[,i])$p.value
         }))
         names(p_values) <- colnames(m_pos)
         p_values <- sort(p_values)
         q_values <- p.adjust(p_values,'bonferroni')

         out <- names(q_values)[ q_values < max.qvalue ]

         if(!is.null(sel.top.n.features)){
            out <- out[1:sel.top.n.features]
            out <- out[!is.na(out)]
         }

         return(out)
      }

   } else {
      stop('Invalid `method` specified')
   }

   if(is.logical(y)){
      x_split <- split(as.data.frame(x), y)

      feature_whitelist <- whitelistByPairwiseTest(x_split[['TRUE']], x_split[['FALSE']])

   } else if(is.factor(y)){
      x_split <- split(as.data.frame(x), y)
      classes <- levels(y)
      features <- colnames(x)

      l <- lapply(1:length(classes), function(i){
         #i=23
         if(verbose){ message('Processing [',i,'/',length(classes),']: ',classes[i]) }
         m_pos <- as.matrix(x_split[[i]])
         m_neg <- as.matrix(do.call(rbind, unname(x_split[-i])))

         whitelistByPairwiseTest(m_pos, m_neg)
      })

      feature_whitelist <- unique(unlist(l))

   } else {
      stop('`y` must be a logical or factor')
   }

   return( x[,colnames(x) %in% feature_whitelist] )
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











