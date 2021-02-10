#Prometheus functions
{
  #funzioni di calcolo
  select_regressors <- function(y,
                                x,
                                maxPowerSingleTerms,
                                maxPowerCrossTerms = maxPowerSingleTerms,
                                maxPowerCrossTermsTerms = maxPowerCrossTerms - 1,
                                maxNumberTerms = NULL,
                                criterion = "AIC",
                                formula = NULL,
                                distance = NULL,
                                kf = 2,
                                testSet = "CV",
                                rowsToTest = NULL,
                                rowsNumber = NULL,
                                includeBackward = FALSE,
                                decomposition = "QR",
                                complexity = 0,
                                verbose = TRUE,
                                glmFamily = gaussian(link = "identity"),
                                weights = rep(1,length(y)),
                                initialModel = NULL){
    start.time <- Sys.time()
    x <- as.matrix(x)
    folds <- NULL
    # creating the folds if cross validation is chosen
    if (criterion == "PREVISION" & testSet == "CV"){
      n <- length(y)
      #calculating the module of the division number of rows/number of groups
      mod <- n %% kf
      #calculating the basic dimension of each groups
      foldSize <- rep(trunc(n / kf,0),kf)
      #the module is given randomly to the groups
      index <- sample((1:kf),mod,replace = FALSE)
      foldSize[index] <-  foldSize[index] + 1
      permutation <- sample(1:n)
      st <- 1
      #folds will contain the groups
      folds <- c()
      for (i in 1:kf){
        folds[i] <- list(permutation[st:(st+(foldSize[i]-1))])
        st <- st + foldSize[i]
      }
    }
    # building up zero model
    {
      if(criterion == "PREVISION" & testSet == "DYNAMIC"){
        n <- length(y)
        rowsToTest <- sample(1:n, rowsNumber, replace = FALSE)
      }
      if (criterion == "PREVISION" & testSet == "CV"){
        tool <- function(parameter){
          return(compute_criterion(criterion = criterion,
                                   y = y,
                                   distance = distance,
                                   rowsToTest = parameter,
                                   decomposition = decomposition,
                                   glmFamily = glmFamily,
                                   weights = weights))
        }
        interceptInfo <- mean(sapply(folds, tool))
      } else {
        interceptInfo <- compute_criterion(criterion = criterion,
                                           y = y,
                                           distance = distance,
                                           rowsToTest = rowsToTest,
                                           formula = formula,
                                           decomposition = decomposition,
                                           glmFamily = glmFamily,
                                           weights = weights)
      }
      #information <- c(interceptInfo + abs(interceptInfo), interceptInfo)
      information <- c(interceptInfo + 10, interceptInfo)
      candidate <- diag(c(rep(1, dim(x)[2])))
      colnames(candidate) <- colnames(x)
      chosenMatrix <- c()
      if(is.null(initialModel)){
        chosenDesignMatrix <- matrix(1, nrow = dim(x)[1], ncol = 1)
      } else {
        for(i in 1:nrow(initialModel)){
          chosenMatrix <- rbind(chosenMatrix, initialModel[i, ])
          candidate <- rbind(initialModel[i, ], candidate)
          candidate <- unique(candidate)
          candidate <- update_candidate(candidate,1,chosenMatrix,maxPowerSingleTerms,maxPowerCrossTerms,maxPowerCrossTermsTerms)
        }
        chosenDesignMatrix <- get_design_matrix(chosenMatrix, x)
      }
      step <- 0
      stop <- 0
      backwardIndicator <- c()
      back_activation <- 0
      numFit <- 0
      numFitB <- 0
      if (interceptInfo > 0){
        information <- c(interceptInfo*(1+complexity) + 10, interceptInfo)
        threshold <- information[step + 1]
      } else {
        information <- c(interceptInfo*(1-complexity) + 10, interceptInfo)
        threshold <- information[step + 1] 
      }
    }
    if (verbose) {
      info = format(round(interceptInfo,2),nsmall=2)
      message = paste0("Initial Model evaluated. ", toupper(criterion), ": ",info,"  Starting selection algorithm.")
      cat(paste0(message,"\n"))
    }
    # adaptive forward selection algorithm
    while (information[step + 2] < threshold  &
           (is.null(maxNumberTerms) ||
            step + 1 - back_activation <= maxNumberTerms) & stop != 1){
      
      if (information[step + 1] > 0){
        threshold = (1 + complexity) * information[step + 1]
      } else {
        threshold = (1 - complexity) * information[step + 1] 
      }
      
      # Forward step
      forwardStep <- estimate_forward(candidate = candidate,
                                      x = x,
                                      y = y,
                                      chosenDesignMatrix = chosenDesignMatrix,
                                      information = information,
                                      criterion = criterion,
                                      numFit = numFit,
                                      distance = distance,
                                      formula = formula,
                                      folds = folds,
                                      testSet = testSet,
                                      rowsNumber = rowsNumber,
                                      rowsToTest = rowsToTest,
                                      decomposition = decomposition,
                                      glmFamily = glmFamily,
                                      weights = weights)
      #update step and information
      step <- step + 1
      information <- forwardStep$information
      numFit <- forwardStep$numFit
      # update the chosen matrix if the conditions of the while cycle are still valid
      if (information[step + 2] < threshold) {
        # Adds the candidate chosen in the current iteration to the matrix of previous chosen predictors
        chosenDesignMatrix <- cbind(chosenDesignMatrix, extract_term(candidate[forwardStep$chosen,], x))
        chosenMatrix <- rbind(chosenMatrix, candidate[forwardStep$chosen,])
        
        # update_candidate gives the candidates for the next iteration
        candidate <- update_candidate(
          candidate,
          forwardStep$chosen,
          chosenMatrix,
          maxPowerSingleTerms,
          maxPowerCrossTerms,
          maxPowerCrossTermsTerms
        )
        # Backward Step
        if (includeBackward) {
          removable <- select_removable(chosenMatrix)
          backwardStep <- estimate_backward(removable = removable,
                                            y = y,
                                            chosenDesignMatrix = chosenDesignMatrix,
                                            information = information,
                                            criterion = criterion,
                                            numFitB = numFitB,
                                            distance = distance,
                                            formula = formula,
                                            folds = folds,
                                            testSet = testSet,
                                            rowsNumber = rowsNumber,
                                            rowsToTest = rowsToTest,
                                            decomposition = decomposition,
                                            glmFamily = glmFamily,
                                            weights = weights)
          backwardIndicator[step] <- !is.null(backwardStep$chosen)
          numFitB <- backwardStep$numFitB
          if (!is.null(backwardStep$chosen)) {
            back_activation <- sum(backwardIndicator)
            information <- backwardStep$information
            candidate <- deupdate_candidate(candidate, chosenMatrix[backwardStep$chosen,])
            chosenMatrix <- chosenMatrix[-backwardStep$chosen,]
            chosenDesignMatrix <- chosenDesignMatrix[,-(backwardStep$chosen + 1)]
          }
        }
        if (verbose) {
          selected_term = paste(paste(
            colnames(chosenMatrix)[which(chosenMatrix[step-back_activation,] != 0)],
            "^", chosenMatrix[step-back_activation, which(chosenMatrix[step-back_activation,] != 0)],sep=""
          ), collapse = " ")
          info = format(round(information[step + 2],2),nsmall=2)
          if (criterion=="RSQUARED"|criterion=="ADJRSQUARED") info = format(round(-information[step + 2],5),nsmall=2)
          if (includeBackward) {
            if (!is.null(backwardStep$chosen)) {
              message <- paste0("Iteration: ", sprintf("%03d",step), "  ", toupper(criterion), ": ",info,"  Selected Term: ", selected_term, "   Backward step performed")
            }  else {
              message <- paste0("Iteration: ", sprintf("%03d",step), "  ", toupper(criterion), ": ",info,"  Selected Term: ", selected_term)
            }
          } else {
            message <- paste0("Iteration: ", sprintf("%03d",step), "  ", toupper(criterion), ": ",info,"  Selected Term: ", selected_term)
          }
          cat(paste0(message,"\n"))
        }
        # Prematurely exits the cycle if there are no more candidates
        if (length(candidate) == 0) {
          stop <- 1
        }
      }
    }
    #exit reason
    exitReason <- ""
    if (information[step + 2] >= threshold){
      exitReason <- "Criterion satisfied"
    }
    if (step + 1 - back_activation > maxNumberTerms){
      exitReason <- "Max number of terms reached"
    } 
    if (stop == 1) {
      exitReason <- "No candidate left"
    } 
    #information managing
    initialInfo <- information[2]
    if (information[length(information)] > information[length(information) - 1]) {
      information <- information[-c(1, 2, length(information))]
    } else {
      information <- information[-c(1, 2)]
    }
    if (criterion =="RSQUARED"|criterion=="ADJRSQUARED"){
      information <- -information
    } 
    #time measuring
    end.time <- Sys.time()
    time.taken <- round(end.time - start.time,2)
    time.taken <- round(end.time - start.time,2)
    if (units(time.taken)=="secs") time <- paste(time.taken,"seconds")
    if (units(time.taken)=="mins") time <- paste(time.taken,"minutes")
    if (units(time.taken)=="hours") time <- paste(time.taken,"hours")
    if (units(time.taken)=="days") time <- paste(time.taken,"days")
    if (units(time.taken)=="weeks") time <- paste(time.taken,"weeks")
    #verbose
    if (verbose) {
      if (includeBackward) {
        message = paste0("Selection procedure completed at iteration ", length(information), ". Backward step activated ", sum(backwardIndicator), " times."," Execution time: ",time,".",sep="")
      } else {
        message = paste0("Selection procedure completed at iteration ", length(information), ".", " Execution time: ", time, ".", sep="")
      }
      cat(paste0(message,"\n"))
    }
    #output
    result <- list(
      "initialInfo" = initialInfo,
      "information" = information,
      "chosenMatrix" = chosenMatrix,
      "backwardIndicator" = backwardIndicator,
      "step" = step,
      "time" = time,
      "numForwardFit" = numFit[-1],
      "numBackwardFit" = numFitB[-1],
      "exitReason" = exitReason
    )
    return(result)
  }
  estimate_forward <- function(candidate,
                               x,
                               y,
                               chosenDesignMatrix,
                               information,
                               criterion,
                               numFit,
                               distance = NULL,
                               formula = NULL,
                               folds = NULL,
                               testSet = "CV",
                               rowsToTest = NULL,
                               rowsNumber = NULL,
                               decomposition = "QR",
                               glmFamily = gaussian(link = "identity"),
                               weights = rep(1,length(y))){
    newInfo <- c()
    # computes all candidate models
    for (l in 1:dim(candidate)[1]) {
      # extracts the potential new term in matricial view
      selectedCandidate <- candidate[l, ,drop=FALSE]
      # extracts the value of the potential new term in each scenario
      newTerm <- extract_term(selectedCandidate, x)
      # adds the vector above as a new column in the design matrix
      designMatrix <- cbind(chosenDesignMatrix, newTerm)
      # computes the statistics of this model
      if(criterion == "PREVISION" & testSet == "DYNAMIC"){
        n <- length(y)
        rowsToTest = sample(1:n,rowsNumber,replace = FALSE)
      }
      if (criterion == "PREVISION" & testSet == "CV"){
        tool <- function(parameter){
          return(compute_criterion(criterion = criterion,
                                   y = y,
                                   x = designMatrix,
                                   distance = distance,
                                   rowsToTest = parameter,
                                   decomposition = decomposition,
                                   glmFamily = glmFamily,
                                   weights = weights))
        }
        newInfo[l] <- mean(sapply(folds, tool))
      } else {
        newInfo[l] <- compute_criterion(criterion = criterion,
                                        y = y,
                                        x = designMatrix,
                                        distance = distance,
                                        rowsToTest = rowsToTest,
                                        formula = formula,
                                        decomposition = decomposition,
                                        glmFamily = glmFamily,
                                        weights = weights)
      }
    }
    chosen <- which.min(newInfo)
    information <- c(information, newInfo[chosen])
    numFit <- c(numFit,dim(candidate)[1])
    results <- list(
      "chosen" = chosen,
      "information" = information,
      "numFit" = numFit
    )
    return(results)
  }
  
  estimate_backward <- function(removable,
                                y,
                                chosenDesignMatrix,
                                information,
                                criterion,
                                numFitB,
                                distance = NULL,
                                formula = NULL,
                                folds = NULL,
                                testSet = "CV",
                                rowsToTest = NULL,
                                rowsNumber = NULL,
                                decomposition = "QR",
                                glmFamily = gaussian(link = "identity"),
                                weights = rep(1,length(y))){
    newInfo <- c()
    chosen <- NULL
    nActivation <- 0
    for (l in removable){
      designMatrix <- as.matrix(chosenDesignMatrix[, - (l + 1)])
      if(criterion == "PREVISION" & testSet == "DYNAMIC"){
        n <- length(y)
        rowsToTest = sample(1:n, rowsNumber, replace = FALSE)
      }
      if (criterion == "PREVISION" & testSet == "CV"){
        tool <- function(parameter){
          return(compute_criterion(criterion = criterion,
                                   y = y,
                                   x = designMatrix,
                                   distance = distance,
                                   rowsToTest = parameter,
                                   decomposition = decomposition,
                                   glmFamily = glmFamily,
                                   weights = weights))
        }
        newInfo[l] <- mean(sapply(folds, tool))
      }else{
        newInfo[l] <- compute_criterion(criterion = criterion,
                                        y = y,
                                        x = designMatrix,
                                        distance = distance,
                                        rowsToTest = rowsToTest,
                                        formula = formula, 
                                        decomposition = decomposition,
                                        glmFamily = glmFamily,
                                        weights = weights)
      }
    }
    numFitB <- c(numFitB, length(removable))
    if (min(newInfo[!is.na(newInfo)]) < information[length(information)]){
      chosen <- which.min(newInfo)
      information[length(information)] <- newInfo[chosen]
    }
    backward <- list(
      "chosen" = chosen,
      "information" = information,
      "numFitB" = numFitB
    )
    return(backward)
  }
  
  compute_criterion <- function(criterion,
                                y,
                                x = NULL,
                                distance = NULL,
                                formula = NULL,
                                rowsToTest = NULL,
                                decomposition = "QR",
                                glmFamily = gaussian(link = "identity"),
                                weights = rep(1,length(y))){
    require(RcppEigen)
    require(fastglm)
    
    rowsToTest <- unlist(rowsToTest)
    criterion <- toupper(criterion)
    decomposition <- toupper(decomposition)
    #checks for errors
    
    # manage fastLmVariable
    if (decomposition == "QR") {
      decompositionType <- 1
    } else if (decomposition == "CHOLESKY") {
      decompositionType <- 2
    } else if (decomposition == "SVD") {
      decompositionType <- 4
    }
    
    if (glmFamily$family == "gaussian" && glmFamily$link == "identity"){
      #if to manage the intercept case
      if (is.null(x)) {
        k <- 1
        res <- y - mean(y)
      } else{
        k <- ncol(x)
        fit <- RcppEigen::fastLmPure(x, y, method = decompositionType)
        res <- fit$residuals
      }
      #some variables used by all criterions
      n <- length(y)
      ssr <- as.numeric(t(res) %*% res)
      #output according to the chosen criterion
      if (criterion == "AIC") {
        computedCriterion <- 2 * (k + 1) + n * (log(ssr / n) + log(2 * pi) + 1)
      } else if (criterion == "BIC") {
        computedCriterion <- log(n) * (k + 1) + n * (log(ssr / n) + log(2 * pi) + 1)
      } else{
        s <- y - mean(y)
        sst <- as.numeric(t(s) %*% s)
        rsquared <- 1 - ssr / sst
        if (criterion == "ADJRSQUARED") {
          computedCriterion <- -(rsquared - (1 - rsquared) * (k - 1) / (n - k))
        } else if (criterion == "RSQUARED") {
          computedCriterion <- -rsquared
        } else if (criterion == "CUSTOM") {
          #here we can add all the variables that we want the user to use
          LogLH <- n * log(1 / n * sst)
          computedCriterion = eval(parse(text = formula))
        } else if (criterion == "PREVISION") {
          rows <- c(1:n)
          #extract the other rows to test the model
          rowsToTrain <- rows[c(-rowsToTest)]
          yTrain <- as.vector(y[rowsToTrain])
          yTest <- as.vector(y[rowsToTest])
          if (is.null(x) == FALSE) {
            #create the matrix to fit the model and the matrix to test the model
            xTrain <- (x[rowsToTrain, , drop = FALSE ])
            xTest <- (x[rowsToTest, , drop = FALSE ])
            #fit the model and compute the coefficients
            fitRMSEP <- RcppEigen::fastLmPure(xTrain, yTrain, method = decompositionType)
            beta <- fitRMSEP$coefficients
            #compute the predicted valus with the xTest matrix
            yPred <- beta %*% t(xTest)
          } else{
            yPred <- rep(mean(yTrain), n - length(rowsToTrain))
          }
          #compute the rmsep
          diff <- yPred - yTest
          err <- as.numeric(diff)
          if (distance == "RMSEP") {
            computedCriterion <- as.numeric(sqrt((t(err) %*% err) / length(yTest)))
          } else if (distance == "MAE") {
            computedCriterion <- mean(abs(err))
          } else if (distance == "WRMAE") {
            computedCriterion <- weighted.mean(abs(err), yTest / sum(yTest))
          } else if (distance == "MINIMAX") {
            computedCriterion <- max(abs(err))
          }
        }
      }
    } else {
      if(is.null(x)){
        x <- as.matrix(rep(1,length(y)))
      }
      fit <- fastglm::fastglmPure(x,y, family = glmFamily, method = decompositionType, weights = weights)
      k <- ncol(x)
      n <- length(y)
      #number of estimated scale parameters (used in aic and bic calculation) depending on the family
      if (glmFamily$family == "Gamma" | glmFamily$family == "inverse.gaussian"){
        s <- 1
      } else {
        s <- 0
      }
      if (criterion == "AIC") {
        computedCriterion <- fit$family$aic(y = y, 
                                            n = length(y), 
                                            mu=  fit$family$linkinv(fit$linear.predictors),
                                            wt = weights,
                                            dev = fit$deviance)+2*length(fit$coefficients)
      } else if (criterion == "BIC") {
        AICcriterion <- fit$family$aic(y = y, 
                                       n = length(y), 
                                       mu=  fit$fitted.values,
                                       wt = weights,
                                       dev = fit$deviance)+2*length(fit$coefficients)
        LogLH <- ((AICcriterion - 2*(k+s))/2)*(-1)
        computedCriterion <- log(n)*(k+s) -2*LogLH
      }else if (criterion == "DEVIANCE"){
        computedCriterion <- sum(fit$family$dev.resids(y = y, mu = fit$fitted.values, wt = weights))
      } else if (criterion == "CUSTOM") {
        #here we can add all the variables that we want the user to use
        AICcriterion <- fit$family$aic(y = y, 
                                       n = length(y), 
                                       mu=  fit$family$linkinv(fit$linear.predictors),
                                       wt = weights,
                                       dev = fit$deviance)+2*length(fit$coefficients)
        LogLH <- ((AICcriterion - 2*(k+s))/2)*(-1)
        computedCriterion = eval(parse(text = formula))
      } else if (criterion == "PREVISION") {
        rows <- c(1:n)
        #extract the other rows to test the model
        rowsToTrain <- rows[c(-rowsToTest)]
        yTrain <- as.vector(y[rowsToTrain])
        yTest <- as.vector(y[rowsToTest])
        wTrain <- as.vector(wt[rowsToTrain])
        wTest <- as.vector(wt[rowsToTest])
        #create the matrix to fit the model and the matrix to test the model
        xTrain <- (x[rowsToTrain, , drop = FALSE ])
        xTest <- (x[rowsToTest, , drop = FALSE ])
        #fit the model and compute the coefficients
        fitRMSEP <- fastglm::fastglmPure(xTrain, yTrain, family = glmFamily, method = decompositionType, weights = wTrain)
        #compute the predicted valus with the xTest matrix
        yPred <- fitRMSEP$fitted.values
        #compute the rmsep
        diff <- yPred - yTest
        err <- as.numeric(diff)
        if (distance == "RMSEP") {
          computedCriterion <- as.numeric(sqrt((t(err) %*% err) / length(yTest)))
        } else if (distance == "MAE") {
          computedCriterion <- mean(abs(err))
        } else if (distance == "WRMAE") {
          computedCriterion <- weighted.mean(abs(err), yTest / sum(yTest))
        } else if (distance == "MINIMAX") {
          computedCriterion <- max(abs(err))
        } else if (distance == "DEVIANCE"){
          computedCriterion <- sum(fitRMSEP$family$dev.resids(y = yTest, mu = yPred, wt = wTest))
        }
      }
    }
    return(computedCriterion)
  }
  check_marginality <- function(v, b, side){
    side <- toupper(side)
    if (side == "LEFT")
      f <- function(a, b)
        isTRUE(all(a <= b))
    if (side == "RIGHT")
      f <- function(a, b)
        isTRUE(all(a >= b))
    result <- apply(
      X = v,
      FUN = f,
      MARGIN = 1,
      b = b
    )
    return(result)
  }
  select_removable <- function(mat){
    rownames(mat) <- 1:dim(mat)[1]
    stepFinal <- c()
    while ((dim(mat)[1] > 1) && (max(rowSums(mat)) != 1) && (length(unique(as.numeric(rowSums(mat)))) > 1)) {
      remove <- c()
      row <- as.numeric(rowSums(mat))
      M <- mat[which(row == max(row)), , drop = FALSE]
      N <- mat[which(row != max(row)), , drop = FALSE]
      v <- matrix(ncol = dim(N)[1], nrow = dim(M)[1])
      for (i in 1:dim(M)[1]) {
        v[i, ] <- check_marginality(N, M[i, ], side = "LEFT")
      }
      for (i in 1:dim(v)[2]) {
        remove[i] <- any(v[, i])
      }
      remove1 <- rownames(N[which(remove == TRUE), , drop = FALSE])
      remove2 <- rownames(M)
      mat <- mat[!rownames(mat) %in% c(remove1, remove2), , drop = FALSE]
      stepFinal <- rbind(stepFinal, M)
    }
    removables <- rbind(mat, stepFinal)
    removableRows <- as.numeric(rownames(removables))
    return(sort(removableRows))
  }
  extract_term <- function(selectedTerm, x){
    require(Rfast)
    terms <- c()
    for (i in 1:max(selectedTerm)) {
      factors <- which(selectedTerm == i)
      terms <- cbind(terms, x[, factors] ^ i)
    }
    computedTerm <- Rfast::rowprods(terms)
    return(computedTerm)
  }
  update_candidate <- function(candidate,chosen,chosenMatrix,maxPowerSingleTerms,maxPowerCrossTerms,maxPowerCrossTermsTerms){
    # Extracts the chosen vector from the matrix
    chosenVector <- candidate[chosen,]
    # For each new possible candidate, check if all conditions are met
    # If yes it is added to the matrix containing the old candidates
    for (i in 1:nrow(chosenMatrix)) {
      potential <- chosenVector + chosenMatrix[i,]
      zero <- matrix(0, nrow = 1, ncol = length(potential))
      # Checks that the term hasn't already been chosen
      if (!(any(apply(chosenMatrix, 1, function(x, want)
        isTRUE(all.equal(x, want)), potential)))) {
        # Checks that the term is admissible according to the chosen parameters
        if (max(potential) <= maxPowerSingleTerms) {
          # count values different from 0 in the potential candidate to see if it is a single term or a cross term
          nonZeros <- sum(potential != 0)
          # if non_zeros <= 1 it is a single term therefore we do not check cross term conditions
          if (nonZeros <= 1 | (
            nonZeros > 1 &
            max(potential) <= maxPowerCrossTermsTerms &
            sum(potential) <= maxPowerCrossTerms
          ))
          {
            # Builds a matrix containing all terms that have to be already chosen for the new term to be added (if the term is x^3, it would be x^1 and x^2)
            potentialMatrix <- matrix(potential,
                                      nrow = 1,
                                      ncol = length(potential))
            val <- list()
            # Creates all possible powers for each risk factor there is one vector for each risk factor
            for (j in 1:length(potential)) {
              val[[j]] <- seq(from = 0, to = potential[j])
            }
            # Builds a matrix with zero in the first row, the potential new candidate in the second row, followed by all possible combinations of powers
            potentialMatrix <- rbind(zero,
                                     potentialMatrix, as.matrix(expand.grid(val)))
            # deletes duplicates and the first two rows therefore ensuring that the potential_matrix does neither include the intercept nor the potential new candidate
            potentialMatrix <- unique(potentialMatrix)
            potentialMatrix <- potentialMatrix[-c(1, 2), , drop = FALSE]
            # Checks if the terms considered have already been chosen aux is true if they have, false if they have not
            aux <- c()
            for (r in 1:nrow(potentialMatrix)) {
              if (any(apply(chosenMatrix, 1, function(x, want)
                isTRUE(all.equal(x, want, check.attributes = FALSE)),
                potentialMatrix[r,]))) {
                aux[r] <- TRUE
              } else{
                aux[r] <- FALSE
              }
            }
            # If all conditions are met the term is added to candidate
            if (all(aux) & length(aux) > 0) {
              candidate <- rbind(candidate, potential)
            }
          }
        }
      }
    }
    # Removes the chosen term from the candidates
    candidate <- candidate[-chosen,,drop=FALSE]
    # Returns the upgraded candidate matrix
    return(candidate)
  }
  deupdate_candidate <- function(candidate, selected){
    marginalityIndicator <- check_marginality(candidate, selected, side = "RIGHT")
    marginals <- which(marginalityIndicator == TRUE)
    if (length(marginals) != 0) {
      candidate <- candidate[-marginals, , drop = FALSE]
    }
    return(candidate)
  }
  
  #funzioni interne
  get_design_matrix <- function(chosenMatrix, x){
    designMatrix <- matrix(ncol = (nrow(chosenMatrix)), nrow = nrow(x))
    designMatrix <-
      apply(
        X = chosenMatrix,
        MARGIN = 1,
        FUN = extract_term,
        x
      )
    return(designMatrix)
  }
  get_selected_terms <- function(chosenMatrix){
    terms <- c()
    for (i in 1:nrow(chosenMatrix)) {
      terms <- rbind(terms, paste(paste(
        colnames(chosenMatrix)[which(chosenMatrix[i,] != 0)],
        "^", chosenMatrix[i, which(chosenMatrix[i,] != 0)],sep=""
      ), collapse = " "))
    }
    return(terms)
  }
  get_regression_model <- function(y, x, chosenMatrix){
    selectedTerms <- get_selected_terms(chosenMatrix)
    designMatrix <- get_design_matrix(chosenMatrix, x)
    colnames(designMatrix) <- selectedTerms
    fullMatrix <- as.data.frame(cbind(y, designMatrix))
    lm <- lm(y ~ ., data = fullMatrix, model = TRUE)
    return(lm)
  }
  
  check_data <- function(y,
                         x,
                         maxPowerSingleTerms,
                         maxPowerCrossTerms = maxPowerSingleTerms,
                         maxPowerCrossTermsTerms = maxPowerCrossTerms - 1,
                         maxNumberTerms = NULL,
                         criterion = "AIC",
                         formula = NULL,
                         distance = NULL,
                         kf = 2,
                         testSet = "CV",
                         rowsToTest = NULL,
                         rowsNumber = NULL,
                         includeBackward = FALSE,
                         decomposition = "QR",
                         complexity = 0,
                         verbose = TRUE,
                         glmFamily = gaussian(link = "identity"),
                         weights = rep(1,length(y))){
    #checks on Y
    if (!(is.numeric(unlist(y))))
      stop("y must be numeric") 
    
    #checks on X
    if (!(is.matrix(x)))
      stop("x must be a matrix") 
    if (!(is.numeric(x)))
      stop("x must be a numeric matrix")
    if (!(nrow(x)>ncol(x)))
      stop("x must have a number of rows greater than the number of columns")
    if (!(nrow(x)==length(y)))
      stop("x and y must have the same number of observations")
    
    #checks on maxPowerSingleTerms
    if(!(is.numeric(maxPowerSingleTerms) & maxPowerSingleTerms >= 1))
      stop("maxPowerSingleTerms must be a strictly positive number")
    
    #checks on maxPowerCrossTerms
    if(!(is.numeric(maxPowerCrossTerms) & maxPowerCrossTerms >= 0)) 
      stop("maxPowerCrossTerms must be a positive number")
    
    #checks on maxPowerCrossTermsTerms
    if(!(is.numeric(maxPowerCrossTermsTerms) & maxPowerCrossTermsTerms >= 0))
      stop("maxPowerCrossTermsTerms must be a positive number")
    
    #checks on maxNumberTerms
    if(!is.null(maxNumberTerms))
      if(!(is.numeric(maxNumberTerms) & maxNumberTerms >= 1 ))
        stop("maxNumberTerms must be a strictly positive number, or NULL")
    
    #checks on criterion
    criterion <- toupper(criterion)
    if(!criterion %in% c("AIC", "BIC", "ADJRSQUARED", "RSQUARED", "CUSTOM", "PREVISION"))
      stop("criterion must be one of the following: AIC, BIC, ADJRSQUARED, RSQUARED, CUSTOM, PREVISION")
    
    #checks on CUSTOM da triggerare se criterion=CUSTOM? tbd
    
    #checks on complexity
    if(!(is.numeric(complexity) & complexity >= 0))
      stop("complexity must be a positive number")
    
    #checks on includeBackward
    if (!is.logical(includeBackward))
      stop("includeBackward must be a logical value")
    
    #checks on decomposition
    decomposition<-toupper(decomposition)
    if(!decomposition %in% c("QR", "CHOLESKY", "SVD"))
      stop("decomposition must be one of the following: QR, CHOLESKY, SVD ")
    
    #checks on verbose
    if (!is.logical(verbose))
      stop("verbose must be a logical value")
  }
  
  #prometheus
  pr <- function(y,
                 x,
                 maxPowerSingleTerms,
                 maxPowerCrossTerms = maxPowerSingleTerms,
                 maxPowerCrossTermsTerms = maxPowerCrossTerms - 1,
                 maxNumberTerms = NULL,
                 criterion = "AIC",
                 formula = NULL,
                 distance = NULL,
                 kf = 2,
                 testSet = "CV",
                 rowsToTest = NULL,
                 rowsNumber = NULL,
                 includeBackward = FALSE,
                 decomposition = "QR",
                 complexity = 0,
                 verbose = TRUE,
                 glmFamily = gaussian(link = "identity"),
                 weights = rep(1,length(y)),
                 initialModel = NULL){
    #add check function here
    check_data(y,
               x,
               maxPowerSingleTerms,
               maxPowerCrossTerms,
               maxPowerCrossTermsTerms,
               maxNumberTerms,
               criterion,
               formula,
               distance,
               kf,
               testSet,
               rowsToTest,
               rowsNumber,
               includeBackward = includeBackward,
               decomposition,
               complexity,
               verbose = verbose,
               glmFamily,
               weights)
    
    myRegressors <- select_regressors(y,
                                      x,
                                      maxPowerSingleTerms,
                                      maxPowerCrossTerms,
                                      maxPowerCrossTermsTerms,
                                      maxNumberTerms,
                                      criterion,
                                      formula,
                                      distance,
                                      kf,
                                      testSet,
                                      rowsToTest,
                                      rowsNumber,
                                      includeBackward = includeBackward,
                                      decomposition,
                                      complexity,
                                      verbose = verbose,
                                      glmFamily = glmFamily,
                                      weights = weights,
                                      initialModel)
    myMatrix <- myRegressors$chosenMatrix
    
    if(!is.null(myMatrix)){
      
      selectedTerms <- get_selected_terms(myMatrix)
      
      designMatrix <- get_design_matrix(myMatrix, x)
      
      myModel <- get_regression_model(y, x, myMatrix)
      
    } else {
      
      selectedTerms <- NULL
      
      designMatrix <- NULL
      
      myModel <- NULL
    }

    
    myInitialSettings <- list(y=y,
                              x=x,
                              maxPowerSingleTerms= maxPowerSingleTerms,
                              maxPowerCrossTerms = maxPowerCrossTerms,
                              maxPowerCrossTermsTerms = maxPowerCrossTermsTerms,
                              maxNumberTerms = maxNumberTerms,
                              criterion = criterion,
                              formula = formula,
                              distance = distance,
                              kf = kf,
                              testSet = testSet,
                              rowsToTest = rowsToTest,
                              rowsNumber = rowsNumber,
                              includeBackward = includeBackward,
                              decomposition = decomposition,
                              complexity = complexity,
                              verbose =  verbose,
                              glmFamily = glmFamily,
                              weights = weights)
    
    elements <- list("selected_terms" = selectedTerms,
                     "designMatrix" = designMatrix,
                     "initialInfo" = myRegressors$initialInfo,
                     "information" = myRegressors$information,
                     "chosenMatrix" = myRegressors$chosenMatrix,
                     "backwardIndicator" = myRegressors$backwardIndicator,
                     "step" = myRegressors$step,
                     "time" = myRegressors$time,
                     "numForwardFit" = myRegressors$numForwardFit,
                     "numBackwardFit" = myRegressors$numBackwardFit,
                     "initialSettings" = myInitialSettings,
                     "exitReason" = myRegressors$exitReason
    )
    prometheus_obj <- append(x = myModel, values = elements, after = 0)
    class(prometheus_obj) <- "lm"
    return(prometheus_obj)
  }
}

######################################
########## DRAFT FUNCTIONS ##########
######################################


#summary_function (this is still a draft)
summary.pr <- function(prom, input = TRUE){
  
  #check x must be "prometheus"
  
  #title
  title <- "Polynomial Regression automatic procedure"
  
  #input details
  inp1 <- "input details:"
  inp2 <- paste0("Data matrix: ",length(prom$initialSettings$x[,1])," rows, ", length(prom$initialSettings$x[1,])," columns")
  inp3 <- paste0("Criterion: ",prom$initialSettings$criterion)
  #questo non so se funziona
  if (prom$initialSettings$criterion=="custom") inp3b <- paste0(prom$initialSettings$formula)
  #da capire anche questo
  if (prom$initialSettings$criterion=="PREVISION"){
    inp3b <- paste0("Prevision distance: ", prom$initialSettings$distance)
    inp3c <- paste0("Test set method: ", prom$initialSettings$testSet)
    #only for CV
    if (is.null(prom$initialSettings$rowsToTest)!=FALSE) inp3x <- paste0("Rows to test: Manual")
    else if (prom$initialSettings$testSet == "CV") inp3ca <- paste0("K-fold: ", prom$initialSettings$kf)
    else if (prom$initialSettings$testSet== "DYNAMIC") inp3cb <- paste0("Rows number: ", prom$initialSettings$rowsNumber)
  }
  inp3z <- paste0("Decomposition method: ", prom$initialSettings$decomposition) 
  if (prom$initialSettings$includeBackward==TRUE) inp4 <- "Procedure: stepwise" else inp4 <- "Procedure: forward"
  inp5 <- paste0("Max number of terms: ",prom$initialSettings$maxNumberTerms)
  inp6 <- paste0("Max single terms power: ",prom$initialSettings$maxPowerSingleTerms)
  inp7 <- paste0("Max cross terms power: ",prom$initialSettings$maxPowerCrossTerms)
  inp8 <- paste0("Max cross terms terms power: ",prom$initialSettings$maxPowerCrossTermsTerms)
  
  #output details
  out1 <- "output details:"
  out2 <- paste0("The procedure ended after ", prom$step, " iterations. Reason: ", tolower(prom$exitReason), sep="")
  out3 <- paste0("Number of terms selected: ", length(prom$selected_terms[,1]))
  out4 <- paste0("Total number of forward regressions: ", sum(prom$numForwardFit))
  out5 <- paste0("Computational time: ", prom$time)
  if (prom$initialSettings$includeBackward == TRUE) out6 <- paste0("Number of backward step performed: ", sum(prom$backwardIndicator))
  if (prom$initialSettings$includeBackward == TRUE) out7 <- paste0("Total number of backward regressions: ", sum(prom$numBackwardFit))
  
  #out
  titleText <- cat("\n",title, sep="\n")
  if (prom$initialSettings$includeBackward == FALSE) 
    outputText <- cat("\n", out1, out2, out3, out4, out5, sep = "\n")
  if (prom$initialSettings$includeBackward == TRUE) 
    outputText <- cat("\n", out1, out2, out3, out4, out5, out6, out7, sep = "\n")
  if (input==TRUE){
    if (prom$initialSettings$criterion == "custom") 
      inputText <- cat("\n",inp1, inp2, inp3, inp3b, inp3z, inp4, inp5, inp6, inp7, inp8, sep = "\n")
    if (prom$initialSettings$criterion == "PREVISION"&&is.null(prom$initialSettings$rowsToTest) != FALSE) 
      inputText <- cat("\n",inp1, inp2, inp3, inp3x, inp3z, inp4, inp5, inp6, inp7, inp8, sep = "\n")
    else if (prom$initialSettings$criterion == "PREVISION"&&prom$initialSettings$testSet == "CV") 
      inputText <- cat("\n",inp1, inp2, inp3, inp3b, inp3c, inp3ca, inp3z, inp4, inp5, inp6, inp7, inp8, sep = "\n")
    else if (prom$initialSettings$criterion == "PREVISION"&&prom$initialSettings$testSet == "DYNAMIC") 
      inputText <- cat("\n",inp1, inp2, inp3, inp3b, inp3c, inp3cb, inp3z, inp4, inp5, inp6, inp7, inp8, sep = "\n")
    else inputText <- cat("\n",inp1, inp2, inp3, inp3z, inp4, inp5, inp6, inp7, inp8, sep = "\n")
  }
}

#explore function (to compute other criterions on the pr object, it is still a draft)
explore.pr <- function(prom,
                       criterion,
                       distance,
                       formula,
                       folds,
                       testSet = "CV",
                       rowsToTest = NULL,
                       rowsNumber = NULL,
                       decomposition = "QR"){
  x <- prom$initialSettings$x
  y <- prom$initialSettings$y
  chosenMatrix <- prom$chosenMatrix
  designMatrix <- matrix(1, nrow = dim(x)[1], ncol = 1)
  Info <- matrix(nrow=length(chosenMatrix[,1]),ncol=length(criterion))
  # computes criteion on all selected terms
  for (l in 1:length(chosenMatrix[,1])) {
    Step <- chosenMatrix[l, ,drop=FALSE]
    newTerm <- extract_term(Step, x)
    designMatrix <- cbind(designMatrix, newTerm)
    for (k in 1:length(criterion)){
      if(criterion[k] == "PREVISION" & testSet == "DYNAMIC"){
        n <- length(y)
        rowsToTest = sample(1:n,rowsNumber,replace = FALSE)
      }
      if (criterion[k] == "PREVISION" & testSet == "CV"){
        tool <- function(parameter){
          return(compute_criterion(criterion = criterion[k],
                                   y = y,
                                   x = designMatrix,
                                   distance = distance,
                                   rowsToTest = parameter,
                                   decomposition = decomposition))
        }
        Info[l,k] <- mean(sapply(folds, tool))
      } else {
        Info[l, k] <- compute_criterion(criterion = criterion[k],
                                        y = y,
                                        x = designMatrix,
                                        distance = distance,
                                        rowsToTest = rowsToTest,
                                        formula = formula,
                                        decomposition = decomposition)
      }
      if (criterion[k] == "RSQUARED"|| criterion[k] =="ADJRSQUARED") Info[l,k] <- -Info[l,k]
    }
    
  }
  colnames(Info) <- criterion
  return(Info)
}
