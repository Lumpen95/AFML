ewma <- function(x, alpha = NULL, span = NULL) {
  
  if (!is.null(alpha) && !is.null(span)) {
    stop("alpha 와 span 중 하나만 입력해야 합니다.")
  }
  
  if (is.null(alpha) && is.null(span)) {
    stop("alpha 또는 span 중 하나를 반드시 입력해야 합니다.")
  }
  
  if (!is.null(span)) {
    alpha <- 2 / (span + 1)
  }
  
  n <- length(x)
  res <- numeric(n)
  
  for (t in 1:n) {
    w <- (1 - alpha)^(t:1 - 1) # 가중치 벡터
    res[t] <- sum(w * x[1:t]) / sum(w)
  }
  return(res)
}

ewmsd <- function(x, alpha = NULL, span = NULL, bias = FALSE) {
  
  # --- 입력 검증 ---
  if (!is.null(alpha) && !is.null(span)) {
    stop("alpha 와 span 중 하나만 입력해야 합니다.")
  }
  if (is.null(alpha) && is.null(span)) {
    stop("alpha 또는 span 중 하나를 반드시 입력해야 합니다.")
  }
  if (!is.null(span)) {
    alpha <- 2 / (span + 1)
  }
  
  n <- length(x)
  ewma_mean <- ewma(x, alpha = alpha)
  res <- numeric(n)
  
  for (t in 1:n) {
    w <- (1 - alpha)^(t:1 - 1) # 가중치 벡터
    centered_sq <- (x[1:t] - ewma_mean[t])^2
    
    w_sum  <- sum(w)
    w_sum2 <- sum(w^2)
    
    if (bias) {
      denom <- w_sum
    } else {
      denom <- w_sum - w_sum2 / w_sum
    }
    res[t] <- sqrt(sum(w * centered_sq) / denom)
  }
  return(res)
}

### 대칭 CUSUM 필터
# 누적 편차가 임계값을 넘어선 경우 이벤트(구조적 변화) 발생
# x(Numeric vector): 국지적 정상성을 가진 IID 관측값(로그 수익률)
# h(Numeric vector): 구조적 변화를 판단할 수 있는 동적 임계값(지수가중 이동표준편차)
# return(Numeric vector): 이벤트 발생 Index
getTEvents <- function(x, h) {
  
  if (!is.numeric(x) || !is.vector(x)) stop("x는 정수형 벡터여야 합니다.")
  if (!is.numeric(h) || !is.vector(h)) stop("h는 정수형 벡터여야 합니다.")
  if (length(x)!=length(h)) stop("x와 h는 길이가 같아야 합니다.")
  
  tEvents <- c()
  sPos <- 0;sNeg <- 0
  diff <- diff(x)
  
  for (i in 2:length(diff)) {
    if (is.na(h[i])) next
    
    sPos <- max(0, sPos + diff[i])
    sNeg <- min(0, sNeg + diff[i])
    if (sNeg < -h[i]) {
      sNeg <- 0
      tEvents <- c(tEvents, i)
    } else if (sPos > h[i]) {
      sPos <- 0
      tEvents <- c(tEvents, i)
    }
  }
  return(tEvents)
}

### 수직 배리어(t1) 함수
# 이벤트 발생(t0) 이후 특정일 반환
# t0(Date scalar): 이벤트 발생일
# index(Date vector): 전체 이벤트 발생일
# num_days(Integer scalar): 최대 보유 기간
add_vertical_barrier <- function(t0, index, num_days = 5) {
  t1_temp <- t0 + days(num_days)
  t1 <- index[index >= t1_temp][1]
  if (is.na(t1)) t1 <- max(index)
  return(t1)
}

### 메타 레이블링 함수
# ptsl(Numeric vector, length=2): 상단/하단 배리어 너비
apply_pt_sl_on_t1 <- function(ret, trgt, side, ptsl) {
  
  if (side == 0) return(NA_integer_)
  
  pt <- ptsl[1] * trgt
  sl <- -ptsl[2] * trgt
  
  cum_ret <- cumsum(ret * side)
  
  pt_hit <- which(cum_ret >= pt)
  sl_hit <- which(cum_ret <= sl)
  
  pt_idx <- if (length(pt_hit) > 0) pt_hit[1] else Inf
  sl_idx <- if (length(sl_hit) > 0) sl_hit[1] else Inf
  
  if (pt_idx < sl_idx) {
    return(1) 
  } else {
    return(0)
  }
  
  final_ret <- tail(cum_ret,1)
  if (final_ret > 0) return(1) else return(0)
}

getConcurrentBar <- function(closeidx, events) {
  # closeidx: 전체 시계열 index vector
  # events: t0, t1 칼럼으로 구성된 data frame
  t_start <- min(events$t0)
  t_end <- max(events$t1)
  
  idx <- closeidx[closeidx >= t_start & closeidx <= t_end] # 첫 이벤트와 마지막 이벤트 기간
  count <- rep(0L, length(idx))
  names(count) <- idx
  
  for (i in seq_len(nrow(events))) {
    overlap <- (idx >= events$t0[i]) & (idx <= events$t1[i])
    count[overlap] <- count[overlap] + 1L
  }
  return(count)
}

getAvgUniqueness <- function(events, numCoEvents) {
  wght <- numeric(nrow(events))
  names(wght) <- events$t0
  
  for (i in seq_len(nrow(events))) {
    t_in <- events$t0[i]
    t_out <- events$t1[i]
    
    ct_slice <- numCoEvents[
      names(numCoEvents) >= t_in & names(numCoEvents) <= t_out
    ]
    
    ct_slice <- ct_slice[ct_slice > 0]
    wght[i] <- mean(1 / ct_slice)
  }
  
  return(wght)
}

getIndMatrix <- function(index, events) {
  # index: 전체 시계열 index vector
  # events: t0, t1 칼럼으로 구성된 data frame
  I <- nrow(events)
  indM <- matrix(
    0,
    nrow = length(index),
    ncol = I,
    dimnames = list(as.character(index), seq_len(I))
  )
  
  for (i in seq_len(I)) {
    t0 <- events$t0[i]
    t1_i <- events$t1[i]
    
    indM[index >= t0 & index <= t1_i, i] <- 1
  }
  return(indM)
}

getAvgUniqueness_indM <- function(indM) {
  c <- rowSums(indM) # 시점별 공존 레이블 개수
  u <- indM / c # 시점별 고유도 행렬
  u[indM == 0] <- NA
  avgU <- colMeans(u, na.rm = T)
  return(avgU)
}

### 수익률 기여도 함수
getReturnAttrWeight <- function(index, events, log_ret, indM) {
  
  # 1. 동시 이벤트 개수
  c_t <- rowSums(indM)
  # 2. 이벤트 개수
  I <- nrow(events)
  # 3. 이벤트별 w 계산
  w_raw <- numeric(I)
  
  for (i in seq_len(I)) {
    t0 <- events$t0[i]
    t1 <- events$t1[i]
    
    idx <- which(index>=t0 & index<=t1)
    w_raw[i] <- sum(log_ret[idx]/c_t[idx],na.rm=TRUE)
  }
  # 4. 절대값
  w_raw <- abs(w_raw)
  # 5. 정규화
  w <- w_raw * I / sum(w_raw,na.rm=TRUE)
  
  return(w)
}

### 시간 감쇠 함수
## 사용자 정의 파라미터 c
## 가장 오래된 이벤트가 최신 이벤트 대비 어느 정도 중요도를 갖는지 의미
## 가장 최근 이벤트 가중값은 1이고, 가장 오래된 이벤트의 가중값이 c
# c = 1: 시간 감쇠 없음
# 0 < c < 1: 시간에 따라 가중값 선형 감쇠, 모든 가중값이 양수
# c = 0: 가중값이 선형으로 0에 수렴(가장 오래된 이벤트 가중치 0)
# c < 0: 오래된 이벤트 일부 완전 삭제
getTimeDecay <- function(avgU, c = 1) {
  x <- cumsum(avgU)
  T <- max(x)
  
  if (c >= 0) {
    slope <- (1 - c) / T 
  } else {
    slope <- 1 / ((c + 1) * T)
  }
  
  intercept <- 1 - slope * T
  
  d <- intercept + slope * x
  d[d < 0] <- 0
  return(d)
}

# fractionally diff
getWeights_ffd <- function(d, thres = 1e-5) {
  # d: 분수 미분 계수
  # thres: 임계값
  w <- c(1) # 초기 가중값은 무조건 1이다.
  k <- 1 # 차수
  
  repeat {
    w_k <- -w[k] * (d - k + 1) / k
    if (abs(w_k) < thres) break
    w <- c(w, w_k)
    k <- k + 1
  }
  
  return(w)
}

fracdiff_ffd <- function(series, d, thres = 1e-5) {
  w <- getWeights_ffd(d, thres)
  width <- length(w) - 1
  
  out <- rep(NA, length(series))
  
  for (i in (width + 1):length(series)) {
    window <- series[(i - width):i]
    out[i] <- sum(rev(w) * window)
  }
  return(out)
}

pcaWeights <- function(cov,riskDist=NULL,riskTarget=1) {
  eig <- eigen(cov)
  idx <- order(eig$values,decreasing = T)
  eVal <- eig$values[idx]
  eVec <- eig$vectors[,idx]
  
  if (is.null(riskDist)) {
    riskDist <- rep(0,dim(cov)[1])
    riskDist[length(riskDist)] <- 1
  }
  
  loads <- riskTarget*sqrt(riskDist/eVal)
  wghts <- eVec%*%loads
  return(wghts)
}

PurgedKFold <- function(events, n_splits=3, pctEmbargo=0) {
  stopifnot(all(c("t0","t1") %in% colnames(events)))
  
  N <- nrow(events)
  indices <- seq_len(N)
  
  mbrg <- floor(N*pctEmbargo)
  
  folds <- split(
    indices,
    cut(indices, breaks = n_splits, labels = FALSE)
  )
  
  cv <- vector("list", length(folds))
  
  for (k in seq_along(folds)) {
    test_idx <- folds[[k]]
    test_t0 <- events$t0[min(test_idx)]
    test_t1_max <- max(events$t1[test_idx])
    
    train_idx_left <- which(events$t1 < test_t0)
    
    embargo_start <- which(events$t0 > test_t1_max)
    if (length(embargo_start) > 0) {
      embargo_start <- embargo_start[1] + mbrg
      train_idx_right <- if (embargo_start <= N) {
        seq(embargo_start, N)
      } else {
        integer(0)
      }
    } else {
      train_idx_right <- integer(0)
    }
    
    train_idx <- c(train_idx_left, train_idx_right)
    
    cv[[k]] <- list(
      train = train_idx,
      test = test_idx
    )

  }
  names(cv) <- paste0("fold",seq_along(cv))
  return(cv)
}

sequentialBootstrap <- function(indM, sample_size=NULL) {
  
  T <- nrow(indM)
  N <- ncol(indM)
  
  phi <- integer(0)
  c_t <- numeric(T)
  
  selected <- integer(0)
  all_cols <- seq_len(ncol(indM))
  
  for (k in seq_len(sample_size)) {
    
    incrU <- numeric(N)
    
    for (j in seq_len(N)) {
      overlap <- indM[,j] == 1
      if (!any(overlap)) {
        incrU[j] <- 0
      } else {
        incrU[j] <- mean(1 / (c_t[overlap] + 1))
      }
    }
    
    prob <- incrU / sum(incrU)
    new_j <- sample(seq_len(N), size=1, prob=prob)
    
    phi <- c(phi, new_j)
    c_t <- c_t + indM[,new_j]
  }
  return(phi)
}

getBootstrappedTrainIdx <- function(cv, indM, sample_size=NULL, num.threads) {
  
  plan(multisession, workers = num.threads)
  
  boot_indices <- future_map(
    seq_along(cv),
    function(k) {
      train_idx <- cv[[k]]$train
      
      indM_tr <- indM[,train_idx,drop=FALSE]
      
      if (is.null(sample_size)) {
        sample_size <- length(train_idx)
      }
      
      boot_cols <- sequentialBootstrap(
        indM = indM_tr,
        sample_size = sample_size
      )
      
      train_idx[boot_cols]
    },
    .options = furrr_options(seed = TRUE)
  )
  plan(sequential)
  
  names(boot_indices) <- paste0("fold", seq_along(cv))
  return(boot_indices)
}

### 특징 중요도 분석
## 성능지표로 mn_log_loss 사용
# mn_log_loss는 accuracy는 특정 임계값 이상이면 모두 동일하게 보는 것과 다르게 확신의 정도를 고려
# 0에 가까울수록 완벽한 모델(작을수록 좋음)
featImportance <- function(data, 
                           cv,
                           cv_replace = TRUE,
                           num_trees = 1000,
                           sample_fraction = 1,
                           case_weights,
                           method = c("MDI", "MDA", "SFI")
) {
  
  imp_fit <- ranger::ranger(bin~., 
                            data = data, 
                            num.trees = num_trees, 
                            mtry = 1, # 선택 변수 1개로 제한
                            sample.fraction = 1, 
                            replace = TRUE,
                            case.weights = case_weights, 
                            importance = "none", 
                            num.threads = 1, 
                            oob.error = TRUE, 
                            probability = TRUE, 
                            node.stats = TRUE)
  
  oob <- imp_fit$prediction.error
  
  ## 특성 중요도 계산
  if (method == "MDI") {
    imp <- featImpMDI(imp_fit)
    oos <- cvScore(data = data, cv = cv, cv_replace = cv_replace, num_trees = num_trees, case_weights = case_weights)
  } else if (method == "MDA") {
    res <- featImpMDA(data = data, cv = cv, cv_replace = cv_replace, num_trees = num_trees, case_weights = case_weights)
    imp <- res$imp
    oos <- res$oos
  } else if (method == "SFI") {
    imp <- featImpSFI(data, cv, case_weights)
    oos <- cvScore(data = data, cv = cv, cv_replace = cv_replace, num_trees = num_trees, case_weights = case_weights)
  }
  
  return(list(
    Importance = imp, 
    OutofBag = oob, 
    OutofSample = oos))
}

### OutofSample 성능 지표 함수
### PurgedKFold 스케줄 제공 필요
cvScore <- function(data, cv, cv_replace, num_trees, sample_fraction = 1, case_weights) {
  
  score <- numeric(length(cv))
  
  for (i in seq_along(cv)) {
    
    train_idx <- cv[[i]]$train_idx
    test_idx <- cv[[i]]$test_idx
    
    train <- data[train_idx,]
    test <- data[test_idx,]
    
    train_w <- case_weights[train_idx]
    test_w <- case_weights[test_idx]
    
    fit <- ranger::ranger(bin~., 
                          data = train, 
                          num.trees = num_trees, 
                          mtry = 1, # 선택 변수 1개로 제한
                          sample.fraction = sample_fraction, 
                          replace = cv_replace,
                          case.weights = train_w,
                          num.threads = 1,
                          probability = TRUE)
    
    pred <- predict(fit, test)$predictions[,"1"]
    eps <- 1e-15
    pred <- pmin(pmax(pred, eps),1-eps)
    
    eval_df <- tibble(
      actual = test$bin,
      pred = pred,
      w = test_w
    )
    score[i] <- -yardstick::mn_log_loss(eval_df, truth = actual, pred, case_weights = w)$.estimate
  }
  return(mean(score))
}

### MDI(Mean Decrease Impurity)
## 각 특성이 얼마나 불순도를 감소시키는지 측정
## 값이 클수록 중요한 특성으로 간주
featImpMDI <- function(rf_fit) {
  
  featNames <- rf_fit$forest$independent.variable.names
  numTrees <- rf_fit$num.trees
  # 특성 중요도 행렬 생성 
  imp_mat <- matrix(0,
                    nrow = numTrees,
                    ncol = length(featNames),
                    dimnames = list(
                      paste0("tree_", seq_len(numTrees)),
                      featNames
                    ))
  # tree별 특징 불순도 감소량 합 계산
  for (tid in seq_len(numTrees)) {
    ti <- ranger::treeInfo(rf_fit,tid)
    ti <- ti[!ti$terminal,]
    
    if (nrow(ti) == 0) next
    
    agg <- aggregate(
      splitStat ~ splitvarName,
      data = ti,
      sum
    )
    imp_mat[tid, agg$splitvarName] <- agg$splitStat
  }
  # 불순도 감소량이 0 인 경우는 특징으로 선택도지 않은 경우이기에 계산에서 제외
  imp_mat[imp_mat == 0] <- NA
  
  imp <- data.frame(
    mean = colMeans(imp_mat, na.rm = TRUE),
    std = apply(imp_mat, 2, sd, na.rm = TRUE) / sqrt(numTrees)
  )
  # 정규화
  imp <- imp / sum(imp$mean, na.rm = TRUE)
  imp <- rownames_to_column(imp, var = "feature")
  
  return(imp)
}

### MDA(Mean Decrease Accuracy, Permutation Importance)
## 각 특성별로 순열(무작위 셔플) 후 계산된 성능과 기존 성능을 비교하여 중요도 결정
# 단일 특성을 permuation 했을 때, 기존 성능이 얼마나 나빠지는가
# 값이 클수록 정보기여도가 큰 특성
# 음수이면 모델을 해치는 특성
featImpMDA <- function(data, cv, case_weights, num_trees = 1000, sample_fraction = 1, cv_replace = TRUE) {
  
  features <- setdiff(colnames(data), "bin")
  
  base_score <- numeric(length(cv))
  perm_score <- matrix(NA, nrow = length(cv), ncol = length(features))
  colnames(perm_score) <- features
  
  for (i in seq_along(cv)) {
    
    train_idx <- cv[[i]]$train_idx
    test_idx <- cv[[i]]$test_idx
    
    train <- data[train_idx,]
    test <- data[test_idx,]
    
    train_w <- case_weights[train_idx]
    test_w <- case_weights[test_idx]
    
    fit <- ranger::ranger(bin~., 
                          data = train, 
                          num.trees = num_trees, 
                          mtry = 1, # 선택 변수 1개로 제한
                          sample.fraction = sample_fraction, 
                          replace = cv_replace,
                          case.weights = train_w,
                          num.threads = 1,
                          probability = TRUE)
    
    base_pred <- predict(fit, test)$predictions[,"1"]
    
    # log loss 무한대 방지 장치(max: 1-eps, min: eps)
    eps <- 1e-15
    base_pred <- pmin(pmax(base_pred, eps),1-eps)
    
    base_df <- tibble(
      actual = test$bin,
      pred = base_pred,
      w = test_w
    )
    
    base_score[i] <- -yardstick::mn_log_loss(base_df, truth = actual, pred, case_weights = w)$.estimate
    
    ## permutation
    for (j in seq_along(features)) {
      
      test_perm <- test
      test_perm[[features[j]]] <- sample(test_perm[[features[j]]]) # 선택된 특성 j 데이터 무작위 샘플링
      
      perm_pred <- predict(fit, test_perm)$predictions[,"1"]
      perm_pred <- pmin(pmax(perm_pred, eps),1-eps)
      
      perm_df <- tibble(
        actual = test_perm$bin,
        pred = perm_pred,
        w = test_w
      )
      
      perm_score[i, j] <- -yardstick::mn_log_loss(perm_df, truth = actual, pred, case_weights = w)$.estimate
    }
  }
  
  ## aggregation
  # base_score: 순열 전 oos score
  # perm_score: 특정 칼럼 순열 후 oos score
  # imp = (perm_score - base_score) / abs(base_score)
  imp_norm <- sweep(-perm_score, 1, base_score, FUN = "+") / -perm_score
  
  imp <- tibble(
    feature = colnames(imp_norm),
    mean = colMeans(imp_norm, na.rm = TRUE),
    std = apply(imp_norm, 2, sd, na.rm = TRUE) / sqrt(nrow(imp_norm))
  )
  return(list(imp = imp, oos = mean(base_score)))
}

featImpSFI <- function(data, cv, case_weights, num_trees = 1000, sample_fraction = 1, cv_replace = TRUE) {
  
  features <- setdiff(colnames(data), "bin")
  
  sfi_mat <- matrix(NA, nrow = length(cv), ncol = length(features))
  colnames(sfi_mat) <- features
  
  for (j in seq_along(features)) {
    
    X_single <- data[, c("bin", features[j])]
    
    scores <- cvScore(
      data = X_single,
      cv = cv,
      cv_replace = cv_replace,
      num_trees = num_trees,
      sample_fraction = sample_fraction,
      case_weights = case_weights
    )
    sfi_mat[, j] <- scores
  }
  
  imp <- tibble(
    feature = colnames(sfi_mat),
    mean = colMeans(sfi_mat, na.rm = TRUE),
    std  = apply(sfi_mat, 2, sd, na.rm = TRUE) / sqrt(nrow(sfi_mat))
  )
  return(imp)
}

apply_etf_trick <- function(
    prices,
    rebal_freq=c("months","quarters","years"),
    riskDist=NULL,
    lookback=252
) {
  
  assets <- colnames(prices)
  
  # 1. 로그수익률
  ret <- diff(log(prices)) |> na.omit()
  
  ## ETF NAV xts
  K <- xts(rep(NA,nrow(ret)),order.by = index(ret))
  colnames(K) <- "NAV"
  K[1] <- 1 # 최초AUM=1
  
  h <- matrix(0,nrow=nrow(ret),ncol=length(assets))
  colnames(h) <- assets
  
  ## 리밸런싱 날짜 index
  rebal_dates <- index(ret)[endpoints(ret,on = rebal_freq)]
  rebal_dates <- rebal_dates[rebal_dates>index(ret)[lookback]]
  
  rebal_w <- lapply(rebal_dates, function(d) {
    window <- ret[paste0("/",d)] |> tail(lookback)
    cov_mat <- cov(window)
    w <- c(pcaWeights(cov_mat,riskDist = riskDist))
  })
  
  rebal_xts <- xts(do.call(rbind, rebal_w),order.by = rebal_dates)
  
  for (t in 2:nrow(ret)) {
    date <- index(ret)[t]
    date1 <- index(ret)[t-1]
    
    if (date %in% rebal_dates) {
      w <- as.numeric(rebal_xts[date,])
      w <- w/sum(abs(w))
      h[t,] <- (K[t-1]%*%w)/as.numeric(prices[date,assets])
    } else {
      h[t,] <- h[t-1,]
    }
    
    p_t0 <- coredata(prices[date,assets])
    p_t1 <- coredata(prices[date1,assets])
    
    pnl <- sum(
      h[t-1,]*(p_t0-p_t1)
    )
    
    K[t] <- K[t-1]+pnl
    
  }
  return(K)
}