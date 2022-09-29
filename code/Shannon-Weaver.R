
Shannon.Weaver.QT <- function(dados, Min = NULL, Max = NULL, Mean = NULL, Sd = NULL, formula = 1){
if(formula == 1){
dados <- as.numeric(dados)

if(is.null(Mean) && is.null(Sd)) {
  Mean <- mean(dados, na.rm = T)
  Sd <- sd(dados, na.rm = T)
}

if(is.numeric(Mean) && is.numeric(Sd)){
  if(is.numeric(dados)){
    P1 <- (sum(dados <= (Mean - Sd), na.rm = T))/(length(dados) - sum(is.na(dados)))

    P2 <- (sum(dados <= (Mean - Sd), na.rm = T) - sum(dados < (Mean - Sd), na.rm = T))/(
      length(dados) - sum(is.na(dados)))

    P3 <- (sum(dados <= (Mean), na.rm = T) - sum(dados < (Mean - Sd/2), na.rm = T))/(
      length(dados) - sum(is.na(dados)))

    P4 <- (sum(dados <= (Mean + Sd/2), na.rm = T) - sum(dados < (Mean), na.rm =T))/(
      length(dados) - sum(is.na(dados)))

    P5 <- (sum(dados <= (Mean + Sd), na.rm = T) - sum(dados < (Mean + Sd/2), na.rm = T))/(
      length(dados) - sum(is.na(dados)))

    P6 <- (length(dados) - sum(dados < (Mean + Sd), na.rm = T) - sum(is.na(dados)))/(
      length(dados) - sum(is.na(dados)))

    Ps <- rbind(P1, P2, P3, P4, P5, P6)
    Ps[Ps == 0] <- (10^-50)
    Agalinha <- -sum(Ps*log(Ps))/log(6)
    return(Agalinha)
  } else {
    print("Please make sure your data has numeric class")
    }
  } else {
  print("Please make sure Mean and Sd arguments are numeric")
  }
}
  else if(formula == 2){
  dados <- as.numeric(dados)

  if(is.null(Min) && is.null(Max)) {
    Min <- min(dados, na.rm = T)
    Max <- max(dados, na.rm = T)
  }

  if(is.numeric(Min) && is.numeric(Max)){
    if(is.numeric(dados)){

  ClasApt <- (Max - Min)/6
  P1 <- (sum(dados <= (Min + 1*ClasApt), na.rm = T))/(length(dados) - sum(is.na(dados)))

  P2 <- (sum(dados <= (Min + 2*ClasApt), na.rm = T) - sum(dados < (Min + 1*ClasApt), na.rm = T))/(
           length(dados) - sum(is.na(dados)))

  P3 <- (sum(dados <= (Min + 3*ClasApt), na.rm = T) - sum(dados < (Min + 2*ClasApt), na.rm = T))/(
           length(dados) - sum(is.na(dados)))

  P4 <- (sum(dados <= (Min + 4*ClasApt), na.rm = T) - sum(dados < (Min + 3*ClasApt), na.rm = T))/(
           length(dados) - sum(is.na(dados)))

  P5 <- (sum(dados <= (Min + 5*ClasApt), na.rm = T) - sum(dados < (Min + 4*ClasApt), na.rm = T))/(
           length(dados) - sum(is.na(dados)))

  P6 <- (sum(dados <= (Max), na.rm = T) - sum(dados < (Min + 5*ClasApt), na.rm = T))/(
           length(dados) - sum(is.na(dados)))

  Ps <- rbind(P1, P2, P3, P4, P5, P6)
  Ps[Ps == 0] <- (10^-50)
  Agalinha <- -sum(Ps*log(Ps))/log(6)
  return(Agalinha)
    } else {
      print("Please make sure your data has numeric class")
    }
  } else {
    print("Please make sure Min and Max arguments are numeric")
  }

} else {print("This fuction offers two options of class delimitaion
limits to estimate Shannon-Weaver index.
Please choose '1' to use mean and sd to delimited Classes or
choose '2' to use the amplitude to delimited the classes")
}
}

Shannon.Weaver.QL <- function(dados, nclass){
  pr <- table(as.character(dados))
  pr <- pr[!is.na(pr)]
  pr2 <- pr/sum(pr)
  Agalinha <- -sum(pr2*log(pr2))/log(nclass)
return(Agalinha)}
