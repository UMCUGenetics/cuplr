## Data for debugging ================================
if(F){
   base_dir <- list(
      hpc='/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/',
      mnt='/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/',
      local='/Users/lnguyen/Documents/projects/P0013_WGS_patterns_Diagn/'
   )

   for(i in base_dir){
      if(dir.exists(i)){
         base_dir <- i
         break
      }
   }

   devtools::load_all(paste0(base_dir,'/CUPs_classifier/processed/cuplr/training/'))
   library(mltoolkit)

   training_data <- read.delim(
      paste0(base_dir,'/CUPs_classifier/processed/cuplr/training/models/0.05a_addRmd/features.txt.gz'),
      check.names=F
   )

   folds <- createCvTrainTestSets(training_data)

}

####################################################################################################
trainRandomForestEnsemble <- function(
   train, test=NULL, colname.response='response',
   feat.sel.method='wilcox', feat.sel.max.qvalue=0.01, feat.sel.top.n.features=NULL,
   calc.imp=T, imp.metric='mda', ntree=200,
   seed=NULL, verbose=1
){
   # train=folds[[1]]$train
   # test=folds[[1]]$test
   # colname.response='response'
   # verbose=2

   if(!is.null(seed)){ set.seed(seed) }

   ##----------------------------------------------------------------------
   if(verbose){ message('> Preparing input features and response variable...') }
   train_data <- dfToFeaturesAndResponse(train, colname.response=colname.response)

   if(!is.null(test)){
      test_data <- dfToFeaturesAndResponse(test, colname.response=colname.response)
   }

   ##----------------------------------------------------------------------
   if(verbose){ message('> Training random forest ensemble...') }
   model <- lapply(1:ncol(train_data$y_ohe),function(i){
      #i='Colon/Rectum'

      if(verbose==2){ message( '[',i,'/',ncol(train_data$y_ohe),']: ', colnames(train_data$y_ohe)[i] ) }

      y <- unname(train_data$y_ohe[,i])

      if(!is.null(feat.sel.method)){
         if(verbose==2){ message('>> Performing feature selection...') }
         x <- univarFeatSel(
            train_data$x, y, method=feat.sel.method,
            max.qvalue=feat.sel.max.qvalue,
            sel.top.n.features=feat.sel.top.n.features
         )
      } else {
         x <- train_data$x
      }

      if(verbose==2){ message('>> Training random forest...') }
      y <- factor(y, levels=c('TRUE','FALSE'))
      randomForest::randomForest(
         x, y,
         strata=y, proximity=F, ntree=ntree,
         #importance=F,
         importance=(calc.imp & imp.metric=='mda'),
         do.trace=F
      )
   })
   names(model) <- colnames(train_data$y_ohe)
   class(model) <- c('list','randomForestEnsemble')

   out <- list()
   out$model <- model

   ##----------------------------------------------------------------------
   if(calc.imp){
      if(verbose){ message('> Calculating feature importance...') }
      exist_features <- unique(unlist(lapply(model, function(i){
         names(i$forest$ncat)
      })))
      exist_features <- colnames(test_data$x)[ colnames(test_data$x) %in% exist_features ] ## Preserve original feature order

      if(imp.metric=='mda'){
         out$imp <- do.call(rbind, lapply(model, function(i){
            #i=model[[1]]
            df <- randomForest::importance(i, type=1)
            v <- structure(df[,1],names=rownames(df))
            v <- v[exist_features]
            names(v) <- exist_features
            v[is.na(v)] <- 0
            return(v)
         }))

      } else {
         if(!is.null(test)){
            out$imp <- permutationImportance(
               model, test_data$x[,exist_features], test_data$y, metric=imp.metric,
               verbose=(verbose==2)
            )
         } else {
            stop('> Feature importance calculation requires test data')
         }
      }
   }

   ##----------------------------------------------------------------------
   if(!is.null(test)){
      if(verbose){ message('> Predicting on test set...') }
      out$test_set <- list(
         probabilities = getPredictions(model, test_data$x, output.classes=F),
         predicted = NA,
         actual = test_data$y
      )

      out$test_set$predicted <- with(out$test_set,{
         factor(
            colnames(probabilities)[max.col(probabilities)],
            levels=colnames(probabilities)
         )
      })
   }

   return(out)
}


####################################################################################################
trainRandomForest <- function(
   train, test=NULL, colname.response='response', do.feat.sel=T,
   calc.imp=T, imp.metric='f1', ntree=200,
   seed=NULL, verbose=1
){
   # train=folds[[1]]$train
   # test=folds[[1]]$test
   # colname.response='response'
   # verbose=2

   if(!is.null(seed)){ set.seed(seed) }

   ##----------------------------------------------------------------------
   if(verbose){ message('> Preparing input features and response variable...') }
   train_data <- dfToFeaturesAndResponse(train, colname.response=colname.response)

   if(!is.null(test)){
      test_data <- dfToFeaturesAndResponse(test, colname.response=colname.response)
   }

   ##----------------------------------------------------------------------
   if(do.feat.sel){
      if(verbose){ message('> Removing useless features...') }
      train_data$x <- univarFeatSel(
         train_data$x, train_data$y,
         verbose=(verbose==2)
      )

      if(!is.null(test)){
         test_data$x <- test_data$x[,colnames(test_data$x) %in% colnames(train_data$x)]
      }
   }

   ##----------------------------------------------------------------------
   if(verbose){ message('> Training...') }
   model <- randomForest::randomForest(
      train_data$x, train_data$y,
      strata=y, proximity=F, ntree=ntree,
      importance=(calc.imp & imp.metric=='mda'),
      do.trace=(verbose==2)
   )

   out <- list()
   out$model <- model

   ##----------------------------------------------------------------------
   if(calc.imp){
      if(verbose){ message('> Calculating feature importance...') }
      if(imp.metric=='mda'){

         out$imp <- suppressWarnings({
            t(randomForest::importance(model, type=1, class=levels(train_data$y)))
         })
      } else {
         if(!is.null(test)){
            out$imp <- permutationImportance(
               model, test_data$x, test_data$y, metric=imp.metric,
               verbose=(verbose==2)
            )
         } else {
            stop('> Feature importance calculation requires test data')
         }
      }
   }

   ##----------------------------------------------------------------------
   if(!is.null(test)){
      if(verbose){ message('> Predicting on test set...') }
      out$test_set <- list(
         probabilities = getPredictions(model,test_data$x, output.classes=F),
         predicted = getPredictions(model,test_data$x, output.classes=T),
         actual = test_data$y
      )
   }

   return(out)
}



####################################################################################################
trainNeuralNet <- function(
   train, test=NULL, colname.response='response',
   feature.groups=NULL, calc.imp=F, imp.metric='f1', model.out.path=NULL,
   seed=NULL, verbose=1
){
   # train=folds[[1]]$train
   # test=folds[[1]]$test
   # colname.response='response'
   # verbose=2
   # feature.groups=feature_groups
   # seed=1

   if(!is.null(seed)){ set.seed(seed) }

   ##----------------------------------------------------------------------
   if(verbose){ message('> Preparing input features and response variable...') }
   train_data <- dfToFeaturesAndResponse(train, colname.response=colname.response)

   if(!is.null(test)){
      test_data <- dfToFeaturesAndResponse(test, colname.response=colname.response)
   }

   ##----------------------------------------------------------------------
   if(verbose){ message('> Removing useless features...') }
   train_data$x <- univarFeatSel(
      train_data$x, train_data$y,
      verbose=(verbose==2)
   )

   if(!is.null(test)){
      test_data$x <- test_data$x[,colnames(test_data$x) %in% colnames(train_data$x)]
   }

   ##----------------------------------------------------------------------
   if(verbose){ message('> Splitting features for branched model...') }
   if(!is.null(feature.groups)){
      feature_names <- colnames(train_data$x)
      train_data$x <- lapply(feature.groups, function(i){
         features <- i[i %in% feature_names]
         train_data$x[,features]
      })

      if(!is.null(test)){
         test_data$x <- lapply(feature.groups, function(i){
            features <- i[i %in% feature_names]
            test_data$x[,features]
         })
      }
   }

   ##----------------------------------------------------------------------
   if(is.null(feature.groups)){
      if(verbose){ message('> Constructing sequential model...') }
      n_inputs <- ncol(train_data$x)
      n_outputs <- ncol(train_data$y_ohe)

      model <-
         keras::keras_model_sequential() %>%

         keras::layer_dense(
            input_shape=n_inputs,
            units=(n_inputs+n_outputs)/2, ## n hidden layers = mean(input, output)
            activation='relu'
         ) %>%

         keras::layer_dense(
            units=n_outputs,
            activation='softmax'
         )

   } else {
      if(verbose){ message('> Constructing multi-input model...') }
      n_inputs <- sapply(train_data$x, ncol)
      n_nodes_upper <- sapply(n_inputs, function(i){ round(i*2/3) }) ## n upper hidden layers = 2/3 * input
      n_outputs <- ncol(train_data$y_ohe)

      ## Upper layers
      layers_input <- lapply(names(n_inputs), function(i){
         keras::layer_input(n_inputs[[i]], name=paste0('input.',i))
      })
      names(layers_input) <- names(n_inputs)

      layers_hidden_upper <- lapply(names(n_inputs), function(i){
         layers_input[[i]] %>%
            keras::layer_dense(
               units=n_nodes_upper[[i]],
               activation='relu',
               name=paste0('hidden_upper.',i)
            )
      })
      names(layers_hidden_upper) <- names(n_inputs)

      model <- layer_concatenate(unname(layers_hidden_upper))

      ## Lower layers
      model <- model %>%

         keras::layer_dense(
            input_shape=n_inputs,
            units=(sum(n_nodes_upper)+n_outputs)/2, ## n lower hidden layers = mean(total hidden layers, outputs)
            activation='relu'
         ) %>%

         keras::layer_dense(
            units=n_outputs,
            activation='softmax'
         )

      model <- keras::keras_model(inputs=layers_input, outputs=model)
   }
   #deepviz::plot_model(model)
   #keras::k_clear_session()

   model <- model %>%
      ## Compile
      keras::compile(
         loss='categorical_crossentropy',
         optimizer='adam',
         metrics='categorical_accuracy'
      )

   model$class_names <- colnames(train_data$y_ohe)

   ##----------------------------------------------------------------------
   if(verbose){ message('> Training...') }
   history <- keras::fit(
      model,
      unname(train_data$x), train_data$y_ohe,
      epoch=100, batch_size=32, validation_split=0.2,
      callbacks=callback_early_stopping(min_delta=0.001, patience=5),
      verbose=(verbose==2)
   )

   ## Initialize outputs
   if(!is.null(model.out.path)){
      keras::save_model_hdf5(model, model.out.path)
   }

   out <- list()
   out$history <- history

   ##----------------------------------------------------------------------
   if(calc.imp){
      if(!is.null(test)){
         if(verbose){ message('> Calculating feature importance...') }
         out$imp <- permutationImportance(
            model, test_data$x, test_data$y, metric=imp.metric,
            verbose=(verbose==2)
         )
      } else {
         stop('> Feature importance calculation requires test data')
      }
   }

   ##----------------------------------------------------------------------
   if(!is.null(test)){
      if(verbose){ message('> Predicting on test set...') }
      out$test_set <- list(
         probabilities = getPredictions(model,test_data$x, output.classes=F),
         predicted = getPredictions(model,test_data$x, output.classes=T),
         actual = test_data$y
      )
   }

   # tab <- table(actual=out$test_set$actual, predicted=out$test_set$predicted)
   # tab <- t(apply(tab,1,function(i){ i/sum(i) }))
   # tab <- round(tab,2)
   # plotMatrixHeatmap(tab, show.labels=T, y.lab='Actual',x.lab='Predicted')

   #keras::k_clear_session()
   return(out)
}

####################################################################################################
crossValidate <- function(
   df, train.func, train.func.args=list(), colname.response='response',
   k=5, seed=1, multi.core=F, verbose=1
){
   #df=training_data
   #colname.response='response'
   #k=5
   #verbose=2
   #train.func=trainNeuralNet
   #train.func.args=list(feature.groups=feature_groups)

   if(!is.list(train.func.args)){ stop('`train.func.args` must be a list') }

   folds <- createCvTrainTestSets(df, stratify.by.col=colname.response, k=k)

   if(verbose){ message('Performing CV...') }
   main <- function(i){
      #i=1
      if(verbose){ message('\n## Fold: ',i) }
      fold <- folds[[i]]

      fold_seed <- as.integer(paste0(seed,i))
      set.seed(fold_seed)

      out <- do.call(
         train.func,
         c(
            list(train=fold$train, test=fold$test, verbose=verbose),
            train.func.args
         )
      )

      out$seed <- fold_seed
      return(out)
   }

   #out <- main(1)
   if(!multi.core){
      lapply(1:length(folds), main)
   } else {
      require(foreach)
      require(parallel)
      require(doParallel)

      #cores <- Sys.getenv("SLURM_NTASKS_PER_NODE")
      cores <- parallel::detectCores()
      if(verbose){ message('Multicore: ',cores,' cores') }

      registerDoParallel(cores=cores)
      foreach(i=1:length(folds)) %dopar% { main(i) }
      #stopCluster(cl)
   }
}

# cv_out <- crossValidate(
#    training_data, train.func=trainNeuralNet,
#    train.func.args=list(feature.groups=feature_groups, calc.imp=T),
#    k=3, verbose=2
# )
#
# cv_out <- crossValidate(
#    training_data, train.func=trainRandomForest,
#    train.func.args=list(calc.imp=T),
#    k=3, verbose=2
# )
# cv_out[[1]]$imp
