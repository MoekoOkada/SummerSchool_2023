computePval <- function(i, j, k,
                        ratioList, ratioSD, ctrlExpMean, objExpMean,
                        ctrlSD, objSD, sampleList, verbose = FALSE) {
  sampleLen <- length(sampleList)
  stdRPval <- dnorm(i, mean = ratioList, sd = ratioSD)
  if (verbose) {
    cat("i:", i, "\t")
    cat("mean:", ratioList, "\t")
    cat("sd:", ratioSD, "\t")
    cat("RPval", stdRPval, "\n")
  }

  # ctrlのread数のp-value
  ctrlPval <- dnorm(j, mean = ctrlExpMean * i, sd = ctrlSD)
  if (verbose) {
    cat("j:", j, "\t")
    cat("mean:", ctrlExpMean * i, "\t")
    cat("sd:", ctrlSD, "\t")
    cat("CtrlPval:", ctrlPval, "\n")
  }

  # objのread数のp-value
  objPval <- dnorm(k, mean = objExpMean * i, sd = objSD)
  if (verbose) {
    cat("k:", k, "\t")
    cat("mean:", objExpMean * i, "\t")
    cat("sd:", objSD, "\t")
    cat("objPval:", objPval, "\n")
  }

  return(stdRPval * ctrlPval * objPval)
}

estimateVar <- function(fit, logMean, i) {
  corMean <- logMean
  corMean[logMean < 0] <- 0
  return(predict(fit, lp(corMean, i)))
}

calcRevisedVar <- function(fit, logMean, expMean, predictVar, i, residue, sampleLen, verbose = FALSE) {
  predictVar <- estimateVar(fit, logMean, i)
  if (verbose) {
    cat("Predict:", logMean, " &", i, " ->", predictVar, "\n")
    cat("ExpPredictVar:", exp(predictVar), "\n")
  }
  predictVar <- exp(predictVar) + residue
  # revisedVar <- sapply(predictVar, function(x) max(x, 1.0))
  # revisedVar <- sapply(1:sampleLen, function(x) max(predictVar[x], expMean[x]*i[x]*0.1,1))
  revisedVar <- sapply(1:sampleLen, function(x) max(predictVar[x], MINVALVAR))
  if (verbose) {
    cat("LogMean:", logMean, "\n")
    cat("i:", i, "\n")
    cat("predictVar", predictVar, "\n")
    cat("revisedVar", revisedVar, "\n")
  }
  return(revisedVar)
}

sigChangeMH <- function(sampleList, numSamples = 1000, verbose = FALSE) {
  validIdx <- is.finite(OCRatioList[sampleList]) &
    is.finite(OCRatioVarList[sampleList]) &
    is.finite(ctrlLoggeomeans[sampleList]) &
    is.finite(ODRatioList[sampleList])
  curSampleList <- sampleList[validIdx]
  sampleLen <- length(curSampleList)
  # 基準の確率を求める．
  i <- ORGRatioList[curSampleList] # OCRatioList[curSampleList] # (0:100)*0.01
  # ctrl のread数
  j <- exp(ctrlLoggeomeans[curSampleList]) * OCRatioList[curSampleList] # 0:exp(ctrlLoggeomeans)
  # objのread数
  k <- exp(objLoggeomeans[curSampleList]) * ODRatioList[curSampleList] # 0:exp(objLoggeomeans)
  # if(is.nan(i) | is.nan(j) | is.nan(k)) return(NaN)
  ratioList <- OCRatioList[curSampleList]
  ratioVar <- sapply(OCRatioVarList[curSampleList], function(x) max(x, MINPVAR))
  ratioSD <- sqrt(ratioVar)
  ctrlLogMean <- ctrlLoggeomeans[curSampleList]
  ctrlExpMean <- exp(ctrlLogMean)
  objLogMean <- objLoggeomeans[curSampleList]
  objExpMean <- exp(objLogMean)
  ctrlVarList <- OCVarList[curSampleList]
  predVar <- estimateVar(fit, ctrlLogMean, ratioList)
  predVar <- exp(predVar)
  for (x in 1:sampleLen) {
    ctrlVarList[x] <- ifelse(ctrlVarList[x] > predVar[x] * 2, predVar[x] * 2,
      ifelse(ctrlVarList[x] < predVar[x] / 2, predVar[x] / 2, ctrlVarList[x])
    )
  }
  # ctrlVarList <- sapply(1:sampleLen, function(x) max(ctrlVarList[x], ctrlExpMean[x]*ratioList[x]*0.1,1))
  objRatioList <- ODRatioList[curSampleList]
  objVarList <- ODVarList[curSampleList]
  predVar <- estimateVar(fit, objLogMean, objRatioList)
  predVar <- exp(predVar)
  for (x in 1:sampleLen) {
    objVarList[x] <- ifelse(objVarList[x] > predVar[x] * 2, predVar[x] * 2,
      ifelse(objVarList[x] < predVar[x] / 2, predVar[x] / 2, objVarList[x])
    )
  }
  # objVarList <- sapply(1:sampleLen, function(x) max(objVarList[x], objExpMean[x]*objRatioList[x]*0.1,1))
  overallRatioList <- ORGRatioList[curSampleList]
  overallRatioSDList <- ORGRatioSDList[curSampleList]

  if (verbose) {
    for (x in curSampleList) {
      cat("Original Value:", x, ":", OCscaledCounts[x, ], ", Var:", OCVarList[x], "\n")
      cat("Object Value:", x, ":", ODscaledCounts[x, ], ", Var:", ODVarList[x], "\n")
      cat("Orig Ratio:", x, ":", counts[x, OCcolumns] / ctrlCounts[x, ], ", Var:", OCRatioVarList[x], "\n")
      cat("Obj  Ratio:", x, ":", counts[x, ODcolumns] / objCounts[x, ], ", Var:", ODRatioVarList[x], "\n")
    }
    cat("CtrlRatio:", ratioList, "[", ctrlLogMean, "]", "\n")
    cat("ObjRatio: ", objRatioList, "[", objLogMean, "]", "\n")
  }

  # fittedVar <- sapply(1:sampleLen, function(x) predict(fit, lp(ctrlLogMean[x], ratioList[x])))
  ctrlFittedVar <- estimateVar(fit, ctrlLogMean, ratioList)
  if (verbose) {
    cat("Predict:", ctrlLogMean, " &", ratioList, " ->", ctrlFittedVar, "\n")
  }

  if (verbose) {
    cat("===\n")
    cat("ctrlLogMean", ctrlLogMean, "\n")
    cat("ctrlExpMean", ctrlExpMean, "\n")
    cat("ctrlRatioList", ratioList, "\n")
    cat("ctrlVarList", ctrlVarList, "\n")
    cat("ctrlFittedVar", exp(ctrlFittedVar), "\n")
  }
  ctrlResidue <- ctrlVarList - exp(ctrlFittedVar)
  if (verbose) {
    cat("ctrlResidue", ctrlResidue, "\n")
  }
  # cat("ctrlExpMean", ctrlExpMean, "\n")
  # ctrlResidue <- sapply(1:sampleLen, function(x) sign(ctrlResidue[x])*min(abs(ctrlResidue[x]), ctrlExpMean[x]))
  # ctrlResidue <- sapply(ctrlResidue, function(x) max(x, 1.0) )
  # cat("ctrlResidue", ctrlResidue, "\n")
  # cat("Residue:", residue, "\n")
  objFittedVar <- estimateVar(fit, objLogMean, objRatioList)
  if (verbose) {
    cat("Predict:", objLogMean, " &", objRatioList, " ->", objFittedVar, "\n")
  }
  if (verbose) {
    cat("objLogMean", objLogMean, "\n")
    cat("objExpMean", objExpMean, "\n")
    cat("objRatioList", objRatioList, "\n")
    cat("objVarList", objVarList, "\n")
    cat("objFittedVar", exp(objFittedVar), "\n")
  }
  objResidue <- objVarList - exp(objFittedVar)
  if (verbose) {
    cat("objResidue", objResidue, "\n")
  }
  # objResidue <- sapply(1:sampleLen, function(x) sign(objResidue[x])*min(abs(objResidue[x]), objExpMean[x]))
  # objResidue <- sapply(objResidue, function(x) max(x, 1.0) )
  if (verbose) {
    # cat("objResidue", objResidue, "\n")
    cat("===\n")
  }

  #
  thresholdPval <- computePval(i, j, k,
    overallRatioList, overallRatioSDList, ctrlExpMean, objExpMean,
    sqrt(ctrlVarList), sqrt(objVarList), curSampleList,
    verbose = verbose
  )
  if (verbose) {
    cat("Threshold:", thresholdPval, "\n", sep = ",")
  }

  # cat("Iteration\n")

  curNum <- 0
  pValList <- c()
  countSignificant <- rep(0, sampleLen)

  while (curNum < numSamples) {
    curNum <- curNum + 1

    i <- rep(NaN, sampleLen)
    invalidIdx <- rep(TRUE, sampleLen)
    # cat("IDX:", invalidIdx, "\n")
    while (sum(invalidIdx) > 0) {
      i[invalidIdx] <- rnorm(sum(invalidIdx),
        mean = overallRatioList[invalidIdx],
        sd = overallRatioSDList[invalidIdx]
      )
      invalidIdx[invalidIdx] <- i[invalidIdx] < 0 | i[invalidIdx] > 1
    }

    # iの値に応じたアップデート (ctrl)
    ctrlPredictRevisedVar <- calcRevisedVar(fit, ctrlLogMean, ctrlExpMean,
      ctrlPredictVar, i, ctrlResidue, sampleLen,
      verbose = verbose
    )
    ctrlPredictRevisedSD <- sqrt(ctrlPredictRevisedVar)

    # cat("--Compute j\n")
    j <- rep(NaN, sampleLen)
    invalidIdx <- rep(TRUE, sampleLen)
    while (sum(invalidIdx) > 0) {
      j[invalidIdx] <- rnorm(sum(invalidIdx),
        mean = ctrlExpMean[invalidIdx] * i[invalidIdx],
        sd = ctrlPredictRevisedSD[invalidIdx]
      )
      invalidIdx[invalidIdx] <- j[invalidIdx] < 0 | j[invalidIdx] > ctrlExpMean[invalidIdx]
    }

    # iの値に応じたアップデート (obj)
    # cat("Compute obj predict var\n")
    # objPredictVar <- predict(fit, lp(objLogMean + log(i), i))
    # objPredictVar <- exp(objPredictVar) + objResidue
    # objPredictRevisedVar <- sapply(objPredictVar, function(x) max(x, 1.0))
    # objPredictRevisedVar <- sapply(1:sampleLen, function(x) max(objPredictVar[x], objExpMean[x]*i[x]*0.1,1))
    objPredictRevisedVar <- calcRevisedVar(fit, objLogMean, objExpMean,
      objPredictVar, i, objResidue, sampleLen,
      verbose = verbose
    )
    objPredictRevisedSD <- sqrt(objPredictRevisedVar)

    k <- rep(NaN, sampleLen)
    # leftidx <- curSampleList[is.finite(objLogMean) & is.finite(i) & is.finite(j)]
    invalidIdx <- rep(TRUE, sampleLen)
    while (sum(invalidIdx) > 0) {
      k[invalidIdx] <- rnorm(sum(invalidIdx),
        mean = objExpMean[invalidIdx] * i[invalidIdx],
        sd = objPredictRevisedSD[invalidIdx]
      )
      # cat("k: last:", sum(invalidIdx), "\n")
      invalidIdx[invalidIdx] <- k[invalidIdx] < 0 | k[invalidIdx] > objExpMean[invalidIdx]
    }

    newPval <- computePval(i, j, k,
      overallRatioList, overallRatioSDList, ctrlExpMean, objExpMean,
      ctrlPredictRevisedSD, objPredictRevisedSD, curSampleList,
      verbose = verbose
    )
    if (verbose) {
      cat("Pval:", newPval, "\n")
      cat("--\n")
    }

    countSignificant <- countSignificant + ifelse(newPval <= thresholdPval, 1, 0)
  }
  finalCountSig <- rep(1.0, length(sampleList))
  finalCountSig[validIdx] <- countSignificant / curNum
  # cat("CurNum: ", curNum, "\n")

  return(finalCountSig)
}
