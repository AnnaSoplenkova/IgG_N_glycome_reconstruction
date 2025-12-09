#оптимизирует порог, по которому следует отбирать пики в дживас 

#Использовали оптимизацию "brent" для поиска порога, который минимизирует дисперсию 
#остатков в линейной регрессии (регрессия вида R2 пика ~ принадлежность пика 
#к "хорошей"/"плохой" группе). Чем лучше порог разделяет значения, тем меньше дисперсия остатков.


R_squares = c(0.55881503, 0.68609639, 0.86630444, 0.82258543, 0.72081415, 0.58726064, 0.72690404, 0.72189759, 0.77729200, 0.81012458, 0.79019724, 0.76815433, 0.61119996, 0.57670609, 0.79019838, 0.32869866, 0.06897064, 0.15330284, 0.31556924, 0.70230255, 0.25854845)

thresholder <- function(threshold, vector) {
  var(residuals(lm(vector ~ ifelse(vector > threshold, "A", "B"))))
}

find_threshold <- function(vector, original_idx = NULL) {
  if (var(vector) != 0 && !is.na(var(vector))) {
    thres <- optim(
      fn = function(th) thresholder(th, vector),
      p = 0,
      method = "Brent",
      lower = min(vector),
      upper = max(vector)
    )
    idx <- which(vector > thres$par)
    return(idx)
  } else {
    # Return original indices if no expansion is possible
    return(original_idx)
  }
}

find_threshold(R_squares)
