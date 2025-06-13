#' @export
fit_with_tuning_real <- function(x, y, alpha = 1.0, nlambda = 50,
                                 lambda.min.ratio = ifelse(nobs <= nvars, 0.01, 1e-04),
                                 lambda = NULL, standardize = TRUE, intercept = TRUE, h=0.25, tau=0.5,
                                 val_ratio=0.2, penalty.factor = rep(1,nvars)){
  np <- dim(x)
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])
  train <- sample(1:nobs, (1-val_ratio)*nobs)
  val <- setdiff(1:nobs,train)
  y <- drop(y)
  lambda_batch_index <- 1
  beta0 <- double(np[2])
  val_errors = rep(Inf, nlambda)
  if(is.null(lambda)){
    lambda_max = 0.0
    ju <- rep(1L, nvars)
    x_list = list("Dense", np[1], np[2], x)
    vp <- as.double(penalty.factor)
    .Call("max_lambda",alpha, x_list, y, as.double(rep(1/nobs, nobs)),
          as.integer(intercept), as.integer(standardize), ju, vp,
          "quantile", h, tau, NULL, F, lambda_max)

    ratio <- lambda.min.ratio^(1/(nlambda-1))
    lambdas = double(nlambda)
    lambdas[1] = lambda_max
    for(i in 2:nlambda){
      lambdas[i] = lambdas[i-1]*ratio
    }
  }else{
    lambdas = lambda
    nlambda = length(lambda)
  }
  fit <- myglmnet(x[train,], y[train], family="quantile", alpha=alpha, nlambda=nlambda,
                lambda=lambdas, penalty.factor = penalty.factor,
                standardize=standardize, intercept=intercept,
                beta0=beta0, h=h, tau=tau)
  for(j in 1:nlambda){
    pred <- x[val, ]%*%fit$beta[, j] + fit$a0[j]
    val_errors[j] = mean(0.5*abs(y[val]-pred)+(tau-0.5)*(y[val]-pred))
    print(sprintf("Validation error at lambda = %f is %f", fit$lambda[j], val_errors[j]))
  }
  best_lambda_index <- which.min(val_errors)
  fit <- myglmnet(x, y, family="quantile", alpha=alpha, nlambda=1,
               lambda=lambdas[best_lambda_index], penalty.factor = penalty.factor,
               standardize=standardize, intercept=intercept,
               beta0=fit$beta[,best_lambda_index], h=h, tau=tau)
  return(fit)
}

fit_with_tuning_snp <- function(genotype.pfile, phenotype_data, covariates, phenotype,
                                alpha = 1.0, nlambda = 50, lambda.min.ratio = ifelse(nobs < nvars, 0.01, 1e-03),
                                lambda = NULL, standardize = TRUE, intercept = TRUE, h=0.25, tau=0.5,
                                val_ratio=0.2, lambda_batch_size=5, penalty.factor = rep(1,nvars)){

  vars <- dplyr::mutate(dplyr::rename(data.table::fread(cmd=paste0('zstdcat', ' ',
                        paste0(genotype.pfile, '.pvar.zst'))), 'CHROM'='#CHROM'),
                        VAR_ID=paste(ID, ALT, sep='_'))$VAR_ID
  psam.ids <- pgenlibr::readIDsFromPsam(paste0(genotype.pfile, '.psam'))
  nobs = length(psam.ids)
  nvars = length(vars) #+ length(covariates)
  # make sure the phe.master has the same individual ordering as in the genotype data
  # so that we don't have error when opening pgen file with sample subset option.
  phenotype_data.sorted <- phenotype_data %>%
    dplyr::left_join(data.frame(ID = psam.ids, stringsAsFactors=F) %>%
      dplyr::mutate(sort_order = 1:n()),
      by='ID'
    ) %>%
    dplyr::arrange(sort_order) %>% dplyr::select(-sort_order) %>%
    data.table::as.data.table()
  train_ids <- sample(psam.ids, (1-val_ratio)*nobs)
  val_ids <- setdiff(psam.ids,train_ids)
  subset_train <- sort(match(train_ids,psam.ids))
  subset_val <- sort(match(val_ids,psam.ids))
  y <- as.vector(phenotype_data.sorted[, ..phenotype][[1]])
  pgen_train <- pgenlibr::NewPgen(paste0(genotype.pfile, '.pgen'),
                                  sample_subset=subset_train)
  m2_train <- ReadList(pgen_train, 1:8000, meanimpute = T)
  snp_matrix_train <- pgenlibr::prepareFeatures(pgen_train,vars,vars)
  snp_matrix_train <- pgenlibr::setcovs(snp_matrix_train,
                    as.matrix(phenotype_data.sorted[subset_train,..covariates]))

  pgen_val <- pgenlibr::NewPgen(paste0(genotype.pfile, '.pgen'),
                                sample_subset=subset_val)
  snp_matrix_val <- pgenlibr::prepareFeatures(pgen_val,vars,vars)
  snp_matrix_val <- pgenlibr::setcovs(snp_matrix_val,
                    as.matrix(phenotype_data.sorted[subset_val,..covariates]))
  snp_matrix_val@xim <- snp_matrix_train@xim

  pgen_all <- pgenlibr::NewPgen(paste0(genotype.pfile, '.pgen'),
                                sample_subset=1:nobs)
  x <- pgenlibr::prepareFeatures(pgen_all,vars,vars)
  m2 <- ReadList(pgen_all, 1:8000, meanimpute = T)
  x <- pgenlibr::setcovs(x,as.matrix(phenotype_data.sorted[1:nobs,..covariates]))

  penalty.factor <- rep(1.0, nvars)
  #penalty.factor[1:length(covariates)] = 0.0

  lambda_batch_index <- 1
  beta0 <- double(nvars)
  val_errors = rep(Inf, nlambda)
  lambda_max = 0.0
  ju <- rep(1L, nvars)
  vp <- as.double(penalty.factor)
  if(x@ncov > 0){
    covmat <- x@covs
    if(intercept) {
      # Maybe the mean and sd should be computed using the
      # observation weights?
      covmean <- apply(covmat, 2, mean)
      covmat <- sweep(covmat, 2, covmean)
    }
    if(standardize) {
      covsd <- apply(covmat, 2, sd)
      covmat <- sweep(covmat, 2, covsd, "/")
      if(!is.null(beta0)){
        beta0[1:x@ncov] <- beta0[1:x@ncov] * covsd
      }
    }
    x_list = list("Plink", dim(x)[1], dim(x)[2], x@ptr, x@xim, x@ncov, covmat)
  } else {
    x_list = list("Plink", dim(x)[1], dim(x)[2], x@ptr, x@xim, 0L)
  }
  .Call("max_lambda",alpha, x_list, y, as.double(rep(1/dim(x)[1], dim(x)[1])),
        as.integer(intercept), as.integer(standardize), ju, vp,
        "quantile", h, tau, NULL, F, lambda_max)

  ratio <- lambda.min.ratio^(1/(nlambda-1))
  lambdas = double(nlambda)
  lambdas[1] = lambda_max
  for(i in 2:nlambda){
    lambdas[i] = lambdas[i-1]*ratio
  }
  lambda_batch_index = 1
  while(lambda_batch_index <= nlambda/lambda_batch_size){
    start_index <- lambda_batch_size*(lambda_batch_index-1)+1
    end_index <- lambda_batch_size*lambda_batch_index
    fit <- myglmnet(snp_matrix_train, y[subset_train], family="quantile", alpha=alpha,
                    nlambda=lambda_batch_size,lambda=lambdas[start_index:end_index],
                    penalty.factor = penalty.factor,
                    standardize=F, intercept=intercept, beta0=beta0, h=h, tau=tau)
    beta0 <- fit$beta[,lambda_batch_size]
    lambda_batch_index <- lambda_batch_index + 1
    pred_result <- PlinkPredict(fit,snp_matrix_val)
    for(j in 1:lambda_batch_size){
      pred <- pred_result[,j]
      val_errors[start_index+j-1] = mean(0.5*abs(y[subset_val]-pred)+(tau-0.5)*(y[subset_val]-pred))
      print(sprintf("Validation error at lambda = %f is %f", fit$lambda[j], val_errors[start_index+j-1]))
    }
    #if(val_errors[end_index]>min(val_errors)){
      #break
    #}
  }
  best_lambda_index <- which.min(val_errors)
  fit <- myglmnet(x, y, family="quantile", alpha=alpha, nlambda=1,
                  lambda=lambdas[best_lambda_index],
                  penalty.factor = penalty.factor,
                  standardize=standardize, intercept=intercept,
                  beta0=beta0, h=h, tau=tau)
  return(fit)
}

