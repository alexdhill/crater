


binary_classifier <- function(y, X, nrepeats = 10, nfolds = 5, final_fold = NA, seed = 1337) {
    `%do%` <- foreach::`%do%`
    set.seed(seed)

    libraries = list("dplyr", "foreach", "caret", "glmnet", "pROC")
    has_libraries <- check_packages(libraries)
    if (!all(has_libraries)) {
        stop(paste("Missing required packages:\n", paste(libraries[!has_libraries], collapse = "\n")))
    }

    if (!(class(y) %in% c('factor', 'numeric', 'integer'))) {
        stop("y must be a factor, numeric, or integer vector")
    }
    if (!is.data.frame(X) & !is.matrix(X)) {
        stop("X must be a data frame or matrix")
    }
    if (length(y) != nrow(X)) {
        stop("Length of y must match number of rows in X")
    }
    if (!(is.numeric(nrepeats) & (nrepeats > -1) & (nrepeats %% 1 == 0))) {
        stop("nrepeats must be a positive integer")
    }
    if (!(is.numeric(nfolds) & (nfolds > 1) & (nfolds %% 1 == 0))) {
        stop("nfolds must be an integer greater than 1")
    }
    if (!(is.na(final_fold) | class(final_fold) %in% c("logical", "numeric", "integer"))) {
        stop("final_fold must be a logical vector")
    }

    if (is.na(final_fold)) {
        final_fold <- caret::createFolds(y, k = nfolds, list = TRUE, returnTrain = TRUE)
    } else if (class(final_fold) %in% c("logical", "numeric", "integer")) {
        if (length(final_fold) != length(y)) {
            stop("final_fold must be the same length as y")
        }
    }

    feature_reduction <- foreach::foreach(1:nrepeats, .combine = 'c') %do% {
        lasso_models <- glmnet::cv.glmnet(
            y = y, x = as.matrix(X), alpha = 1,
            family = "binomial", type.measure = "mse", nfolds = nfolds
        )

        return(list(list(
            mse = lasso_models$cvm[lasso_models$index[1]],
            models = lasso_models
        )))
    }

    best_feature_set <- feature_reduction[[
        which.min(sapply(feature_reduction, function(x) x$mse))
    ]]$models
    best_features <- (best_feature_set$model$beta[, best_feature_set$model$index[1]]) != 0
    if (sum(best_features) < 2) {
        warning("Less than 2 features selected by LASSO; ignoring this step.")
        best_features <- rep(TRUE, ncol(X))
    }

    parameter_tuning <- foreach::foreach(alpha = seq(0, 1, by = 0.01), .combine = 'c') %do% {
        X_best <- X[, best_features]

        tuning_models <- glmnet::cv.glmnet(
            y = y, x = as.matrix(X_best), alpha = alpha, lambda = seq(0.01, 1, by = 0.01),
            family = "binomial", type.measure = "mse", nfolds = nfolds
        )

        return(list(list(
            mse = tuning_models$cvm[tuning_models$index[1]],
            alpha = alpha, lambda = tuning_models$lambda.min,
            models = tuning_models
        )))
    }

    best_parameters <- parameter_tuning[[
        which.min(sapply(parameter_tuning, function(x) x$mse))
    ]]

    X_final <- X[final_fold, best_features]
    final_model <- glmnet::glmnet(
        y = y[final_fold], x = as.matrix(X_final),
        alpha = best_parameters$alpha, lambda = best_parameters$lambda,
        family = "binomial"
    )

    return(list(
        model = final_model,
        features = colnames(X)[best_features],
        alpha = best_parameters$alpha,
        lambda = best_parameters$lambda,
        training_set = final_fold
    ))
}
