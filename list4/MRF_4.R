#-----------------*************** Lista 4 MRF ***************-------------------

#_______________________________________________________________________________
#
#------------------------------Zadanie 4.1--------------------------------------
#
#_______________________________________________________________________________
# Załóżmy, że funkcja stopy terminowej miała postać 
#
#         r(t) = 0,01 + 0,005t,
#
# gdy dziesięcioletnia obligacja zerokuponowa została zakupiona. 
# Pół roku później, gdy sprzedano obligację funkcja stopy terminowej 
# miała postać 
#
#         r(t) = 0,02 + 0,001t.
#
# Jaka była stopa zwrotu z tej inwestycji?
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

t=10

r0 <- function(t){
  0.01 + 0.005 * t
}

r1 <- function(t){
  0.02 + 0.001 * t
}
calka_r0 <- integrate(r0,0,10)$value

calka_r1 <- integrate(r1,0,9.5)$value
  
B0 <- exp(-1*calka_r0) #cena w momencie 0
B1 <- exp(-1*calka_r1)

stopa <- (B1-B0)/ B0
stopa




#_______________________________________________________________________________
#
#------------------------------Zadanie 4.2--------------------------------------
#
#_______________________________________________________________________________
# Załączony plik TermStructureZero.csv zawiera symulowane dane 
# dotyczące 20 obligacji zerokuponowych. 
# W kolejnych polach znajdują się:
#
# Maturity  = liczba lat do wykupu,
# Quote     = kurs obligacji.   
#
# Napisz kod w R, który uzupełni plik o pola:
#
# SpotRate    = stopa spot dla n lat (zerokuponowa)
# ForwardRate = stopa terminowa dla n-tego roku.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

dane <- read.csv("C:\\Users\\olgaf\\OneDrive\\Pulpit\\MRF_2025\\lista 4\\TermStructureZero.csv", header = TRUE, sep = ",")

dane$SpotRate <- (100/dane$Quote)^(1/dane$Maturity)-1

dane$ForwardRate[1] <- dane$SpotRate[1]
for (i in 2:nrow(dane)) {
  dane$ForwardRate[i]<-
    (((1+dane$SpotRate[i])^dane$Maturity[i])/((1+dane$SpotRate[i-1])^dane$Maturity[i-1]) -1)
}

write.csv(dane, "C:\\Users\\olgaf\\OneDrive\\Pulpit\\MRF_2025\\lista 4\\TermStructureZero1.csv", row.names = FALSE)


#_______________________________________________________________________________
#
#------------------------------Zadanie 4.3--------------------------------------
#
#_______________________________________________________________________________
# Załączony plik TermStructure.csv zawiera symulowane dane 
# dotyczące 20 obligacji o kuponach wypłacanych raz w roku 
# (pierwsza wypłata za rok). W kolejnych polach znajdują się:
#
# Maturity  = liczba lat do wykupu,
# FV        = wartość nominalna,   
# Coupon    = wartość kuponu,
# Price     = cena obligacji,
#
# Napisz kod w R, który uzupełni plik o pola:
#
# SpotRate      = stopa spot dla n lat (zerokuponowa)
# ForwardRate   = stopa terminowa dla n-tego roku.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
dane_1 <- read.csv("C:\\Users\\olgaf\\OneDrive\\Pulpit\\MRF_2025\\lista 4\\TermStructure.csv", header = TRUE, sep = ",")

C <- dane_1$Coupon[1]
F <- dane_1$FV[1]
P <- dane_1$Price[1]
n <- dane_1$Maturity[1]
dane_1$SpotRate[1] <- (C + F)/P - 1

for (i in 2:nrow(dane_1)) {
  C <- dane_1$Coupon[i]
  F <- dane_1$FV[i]
  P <- dane_1$Price[i]
  n <- dane_1$Maturity[i]
  PV <- 0
  for (j in 1:(n-1)) {
    PV <- PV + C / (1 + dane_1$SpotRate[j])^j
  }
  dane_1$SpotRate[i] <- ((C + F)/(P - PV))^(1/n) - 1
}

dane_1$ForwardRate[1] <- dane_1$SpotRate[1]
##print(dane_1$SpotRate[10])

for (i in 2:nrow(dane)) {
  dane_1$ForwardRate[i]<-
    (((1+dane_1$SpotRate[i])^dane_1$Maturity[i])/((1+dane_1$SpotRate[i-1])^dane_1$Maturity[i-1]) -1)
}
##print(dane_1$ForwardRate[10])
write.csv(dane_1, "C:\\Users\\olgaf\\OneDrive\\Pulpit\\MRF_2025\\lista 4\\TermStructure1.csv", row.names = FALSE)

#_______________________________________________________________________________
#
#------------------------------Zadanie 4.4--------------------------------------
#
#_______________________________________________________________________________
# Dla notowań spółek giełdowych A,B,C z pliku NotowaniaABC.csv wyznacz
#
# a) oczekiwane zwroty i wariancje tych zwrotów dla każdej ze spółek A,B,C;
# b) oczekiwany zwrot i wariancję zwrotów dla portfela z wagami 
#
#       w_A=0.2,  w_B=0.3,  w_C=0.5.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# podpunkt 4.4a

dane_noto <- read.csv("C:\\Users\\olgaf\\OneDrive\\Pulpit\\MRF_2025\\lista 4\\NotowaniaABC.csv", header = TRUE, sep = ",")
dane_noto$zwrotA <- NA
dane_noto$zwrotB <- NA
dane_noto$zwrotC <- NA

for(i in 2:nrow(dane_noto)){
  dane_noto$zwrotA[i] <- (dane_noto$A[i]-dane_noto$A[i-1])/dane_noto$A[i-1]
  print(dane_noto$zwrotA[i])
  dane_noto$zwrotB[i] <- (dane_noto$B[i]-dane_noto$B[i-1])/dane_noto$B[i-1]
  dane_noto$zwrotC[i] <- (dane_noto$C[i]-dane_noto$C[i-1])/dane_noto$C[i-1]
}

mean(dane_noto$zwrotA, na.rm = TRUE)
mean(dane_noto$zwrotB, na.rm = TRUE)
mean(dane_noto$zwrotC, na.rm = TRUE)

sqrt(var(dane_noto$zwrotA, na.rm = TRUE))
sqrt(var(dane_noto$zwrotB, na.rm = TRUE))
sqrt(var(dane_noto$zwrotC, na.rm = TRUE))

# podpunkt 4.4b

w_A=0.2
w_B=0.3
w_C=0.5

ocz_zwrot <- mean(dane_noto$zwrotA, na.rm = TRUE) * w_A + 
  mean(dane_noto$zwrotB, na.rm = TRUE) * w_B +
  mean(dane_noto$zwrotC, na.rm = TRUE) * w_C

wariancja <- sqrt(var(dane_noto$zwrotA, na.rm = TRUE) * w_A^2 + 
  var(dane_noto$zwrotB, na.rm = TRUE) * w_B^2 +
  var(dane_noto$zwrotC, na.rm = TRUE) * w_C^2 + 
  2*w_A*w_B*cov(dane_noto$zwrotA,dane_noto$zwrotB,use="complete.obs")+
  2*w_A*w_C*cov(dane_noto$zwrotA,dane_noto$zwrotC,use="complete.obs")+
  2*w_C*w_B*cov(dane_noto$zwrotB,dane_noto$zwrotC,use="complete.obs"))

#_______________________________________________________________________________
#
#------------------------------Zadanie 4.5--------------------------------------
#
#_______________________________________________________________________________
# Dla notowań spółek giełdowych A,B,C z pliku NotowaniaABC.csv wyznacz,
# korzystając z odpowiednich twierdzeń z wykładu 4:
#
# a) wagi portfela o minimalnej wariancji;
# b) wagi portfela o minimalnej wariancji i oczekiwanym zwrocie 15%.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dane_ABC <- as.matrix(read.csv("list4/NotowaniaABC.csv"))
#head(dane_ABC)
#str(dane_ABC)
# podpunkt 4.5a
r <- diff(dane_ABC)/dane_ABC[-nrow(dane_ABC), ] #proste stopy zwrotu
C <- cov(r) 
u <- t(rep(1, 3))
C_inv <- solve(C)
w <- (u %*% C_inv) / as.numeric(u %*% C_inv %*% t(u))
w

# podpunkt 4.5b
mi <-  t(colMeans(r))       
miv <- 0.15
miu <- mi %*% C_inv %*% t(u)
mimi <- mi %*% C_inv %*% t(mi) 
umi <- u %*% C_inv %*% t(mi) 
uu <- u %*% C_inv %*% t(u) 

licznik1 <- det(matrix(c(1,miv,miu,mimi), nrow = 2,ncol = 2))
licznik2 <- det(matrix(c(uu,umi,1,miv), nrow = 2,ncol = 2))
licznik <- licznik1* u %*% C_inv + licznik2*mi %*% C_inv
mianownik <- det(matrix(c(uu,umi,miu,mimi), nrow = 2,ncol = 2))
wyn <- licznik/mianownik
wyn



#_______________________________________________________________________________
#
#------------------------------Zadanie 4.6--------------------------------------
#
#_______________________________________________________________________________
# Korzystając z funkcji solve.QP z pakietu quadprog wykonaj zadanie 4.5.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
library(quadprog)
# podpunkt 4.5a za pomocą solve.QP z pakietu quadprog

Dmat <- 2*C #wynika ze wzoru, chcemy minimalizowac w.T * C * w, a w quadprog musimy miec 1/2 przed
dvec <- rep(0, 3) #0 bo nie mamay u nas formy liniowej a jedynie to co wyzej


u <- rep(1, 3)
Amat <- cbind(u)
bvec <- c(1)
meq <- 1

sol <- solve.QP(Dmat, dvec, Amat, bvec, meq)
w <- sol$solution
w



# podpunkt 4.5b za pomocą solve.QP z pakietu quadprog

library(quadprog)

Dmat <- 2*C #wynika ze wzoru, chcemy minimalizowac w.T * C * w, a w quadprog musimy miec 1/2 przed
dvec <- rep(0, 3) #0 bo nie mamay u nas formy liniowej a jedynie to co wyzej

u <- rep(1, 3)
mi <-  colMeans(r)    
Amat <- cbind(u, mi)
bvec <- c(1, miv)
meq <- 2  


sol <- solve.QP(Dmat, dvec, Amat, bvec, meq)
w <- sol$solution
w


#_______________________________________________________________________________
#
#------------------------------Zadanie 4.7--------------------------------------
#
#_______________________________________________________________________________
# Czy (i ewentualnie jak?) zmieni się odpowiedź do zadania 4.6, przy założeniu,
# że krótka sprzedaż nie jest dostępna? Co w przypadku, gdy oczekiwany zwrot 
# będzie równy 10% a nie 15%? 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# podpunkt 4.5a za pomocą solve.QP z pakietu quadprog i bez krótkiej sprzedaży







# podpunkt 4.5b za pomocą solve.QP z pakietu quadprog i bez krótkiej sprzedaży







# podpunkt 4.5b za pomocą solve.QP z pakietu quadprog i bez krótkiej sprzedaży
# i przy oczekiwanym zwrocie równym 10% zamiast 15%






