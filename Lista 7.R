#-----------------*************** Lista 7 MRF ***************-------------------


#_______________________________________________________________________________
#
#------------------------------Zadanie 7.1--------------------------------------
#
#_______________________________________________________________________________
# Cena kontraktu futures na jedną akcję nie wypłacającą dywidend z terminem
# rozliczenia za 3 miesiące wynosi 100. Cena europejskiej opcji kupna 
# z terminem wygaśnięcia za 3 miesiące i ceną wykonania 100 wynosi 5.
# Cena europejskiej opcji sprzedaży z terminem wygaśnięcia za 3 miesiące 
# i ceną wykonania 100 wynosi 6. Wskaż możliwy arbitraż. 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# arbitraż jest gdy Vx_0 =/= Vx_t
X=100 # cena wykonania
Ft=100 # cena kontraktu
t=3
Ce=5
Pe=6
r=1

Vx_0=Ce-Pe
Vx_t=(Ft - X)*exp(-r*t/12)
arbitraz = Vx_t - Vx_0
arbitraz
if (Vx_0 != Vx_t) {
  cat("Możliwy arbitraż: ",arbitraz)
} else {
  print("Nie ma arbitrażu")
}

# Należy kupić opcję kupna oraz sprzedać opcję sprzedaży oraz zająć długą 
# pozycję na kontrakcie futures 

#_______________________________________________________________________________
#
#------------------------------Zadanie 7.2--------------------------------------
#
#_______________________________________________________________________________
# Cena europejskiej opcji kupna o cenie wykonania 60 wygasającej 
# za 6 miesięcy wynosi 4. Cena akcji, na którą została wystawiona opcja, to 58.
# Na akcję wypłacana będzie dywidenda w wysokości 1 za 2 miesiące i w tej samej
# wysokości za 5 miesięcy. Stopa procentowa dla wszystkich okresów to 10%
# (kapitalizacja ciągła). Ile wynosi cena europejskiej opcji sprzedaży 
# o cenie wykonania 60 wygasającej za 6 miesięcy?
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# wzór W8 twierdzenie 2, 3 kropa
Ce = 4
X = 60
t1=6
cena_akcji = 58 # S(0)
r=0.1
#dane do dywidend
dywidenda = 1
t2=2
t3=5

#Szukana Pe z ceną wykonania 60
# Ce - Pe = S(0) - D0 - exp(-rT)*X

# wzór na D) dla kapitalizacji ciąglej
# D0 = D1*exp(-rt1) + D2*exp(-rt2)

D0 = dywidenda*exp(-r*t2/12) + dywidenda*exp(-r*t3/12)
Pe = Ce - cena_akcji + D0 + exp(-r*6/12)*X
Pe

#_______________________________________________________________________________
#
#------------------------------Zadanie 7.3--------------------------------------
#
#_______________________________________________________________________________
# Napisz funkcję (oznaczenia jak na wykładzie)
#
# BlackScholesGreeks(S, X, time, r, sigma),
#
# która zwraca wektor 
#
# c(CE, PE, DeltaCall, DeltaPut, Gamma, Vega), 
#
# gdzie
#
# CE        = cena europejskiej opcji kupna wg wzoru Blacka-Scholesa,
# PE        = cena europejskiej opcji sprzedaży wg wzoru Blacka-Scholesa,
# DeltaCall = współczynnik grecki delta dla europejskiej opcji kupna,
# DeltaPut  = współczynnik grecki delta dla europejskiej opcji sprzedaży,
# Gamma     = współczynnik grecki delta dla opcji europejskiej,
# Vega      = współczynnik grecki vega dla opcji europejskiej.
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# wzory wykład 10

BlackScholesGreeks=function(S, X, time, r, sigma){
  
  d1 <- (log(S/X) + time*(r + 0.5*sigma^2))/(sigma*sqrt(time))
  d2 <- d1 - sigma*sqrt(time)
  CE <- S*pnorm(d1) - X*exp(-r*time)*pnorm(d2)
  PE <- X*exp(-r*time)*(1-pnorm(d2)) - S*pnorm(-d1)
  DeltaCall <- pnorm(d1)
  DeltaPut <- 1-pnorm(d1) 
  Gamma <- dnorm(d1)/(S*sigma*sqrt(time))
  Vega <- S* dnorm(d1)*sqrt(time)
  
  return(c(CE,PE,DeltaCall,DeltaPut, Gamma, Vega))
}

BlackScholesGreeks(S=80, X=100, time=2, r=0.1, sigma=0.05)

#_______________________________________________________________________________
#
#------------------------------Zadanie 7.4--------------------------------------
#
#_______________________________________________________________________________
# Napisz kod funkcji (oznaczenia jak na wykładzie)
#
# CRRPrice(TypeFlag, S, X, u, d, r, N),
#
# obliczającej ceny opcji europejskich kupna (TypeFlag = "ce") 
# i sprzedaży (TypeFlag = "pe") wprost ze wzorów Coxa-Rossa-Rubinsteina. 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# model dwumianowy

CRRPrice<-function(TypeFlag, S, X, u, d, r, n){
  
  p_gwiazdka <- (r-d)/(u-d)
  wyplata=0
  cena<-0
  
  for(k in range(0:n)){
    
    S_n <- S * (1 + u)^k * (1 + d)^(n - k)
    if(TypeFlag=="ce"){wyplata<-max(0,S_n-X)}
    if(TypeFlag=="pe"){wyplata<-max(0,X-S_n)}
    prawdopodobienstwo = choose(n,k)*p_gwiazdka^k*(1-p_gwiazdka)^(n-k)
    cena=cena+prawdopodobienstwo*wyplata
  }
  return(cena/(1 + r)^n) #dyskontujemy bo chcemy cenę dziś

}

CRRPrice(TypeFlag = "pe", S=80, X=80, u=0.1, d=-0.05, r=0.05, n=2) 

#_______________________________________________________________________________
#
#------------------------------Zadanie 7.5--------------------------------------
#
#_______________________________________________________________________________
# Napisz funkcję (oznaczenia jak na wykładzie)
#
# BinomialTreeOption(TypeFlag, S, X, d, r, u, N)
#
# implementując algorytmy wyceny 
#
# - europejskich opcji kupna (TypeFlag = "ce") i sprzedaży (TypeFlag = "pe")oraz
# - amerykańskich opcji kupna (TypeFlag = "ca") i sprzedaży (TypeFlag = "pa")
#
# na drzewie dwumianowym omawiane na wykładzie. 
# W wyniku realizacji kodu powinniśmy otrzymać macierz trójkątną górną wymiaru
# (N+1) x (N+1) oznaczoną Tree, w której kolumnach znajdą się wartości opcji
# w kolejnych chwilach.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


BinomialTreeOption<-function(TypeFlag, S, X, d, r, u, n){
  
  
  
  
  
  
  
  
  
  
  
}



#_______________________________________________________________________________
#
#------------------------------Zadanie 7.6--------------------------------------
#
#_______________________________________________________________________________
# Napisz funkcję
#
# BinomialTreePlot(Tree)
#
# do graficznego przedstawienia drzewa Tree cen opcji otrzymanego 
# za pomocą funkcji z zadania 7.5.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

BinomialTreePlot <- function(Tree){
  
  
  
  
  
  
  
  
  
  
}







#_______________________________________________________________________________
#
#------------------------------Zadanie 7.7--------------------------------------
#
#_______________________________________________________________________________
# Broker wystawił opcje kupna na 100 akcji nie wypłacających dywidendy.
# Parametry opcji: czas do wygaśnięcia = 91 dni, cena wykonania X=42, zmienność
# implikowana akcji sigma=0.5 w skali roku, stopa wolna od ryzyka r=8% 
# w skali roku (kapitalizacja ciągła).
#
# (a) Podaj skład portfela zabezpieczającego (delta hedging) 
#     czyli portfela o wartości równej zeru, składającego się z 
#     w/w opcji (-100), akcji i gotówki (lub długu) na rachunku pieniężnym,
#     jeśli jedna akcja w momencie wystawienia opcji kosztowała 40.
#
# Uwaga: uzupełniając portfel o gotówkę na oprocentowanym rachunku pieniężnym 
#        (w taki sposób, aby wartość początkowa portfela była równa 0) 
#        uwzględniamy koszty finansowania portfela.
#
# (b) Jaka będzie wartość tego portfela po dwóch dniach jeśli cena akcji 
#     wzrosła do 42, a pozostałe parametry (poza czasem do wykupu)
#     się nie zmieniły?
#
# (c) Jaka byłaby wartość portfela niezabezpieczonego po tych dwóch dniach?
#
#     Uwaga: portfel niezabezpieczony składa się z krótkiej pozycji 
#     na opcjach kupna na 100 akcji i gotówki otrzymanej ze sprzedaży tych opcji
#     na rachunku pieniężnym (oprocentowanym wg stopy wolnej od ryzyka).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# podpunkt a









# podpunkt b










# podpunkt c









#_______________________________________________________________________________
#
#------------------------------Zadanie 7.8--------------------------------------
#
#_______________________________________________________________________________
# Napisz funkcję
#
# BrownianMotion(n=1, N, t, mu, sigma, S),
#
# zwracającą macierz U rozmiaru N x n, której wierszami są niezależne trajektorie
# ruchu Browna (na przedziale [0,t] i startujacego z punktu S)
# z parametrami mu oraz sigma.
#
# W szczególności, k-ta wartość dla j-tej trajektorii (tj. S_j(kt/n))
# jest generowana zgodnie z następującym schematem Eulera
#
# S_j(0)=S, S_j(k*t/n) = S_j((k-1)*t/n) + mu * t/n + sigma * \sqrt{t/n} * w_j^(k)
#
# gdzie
#
# k=1,...,n oraz j=1,...,N, natomiast [w_j^{(1)},...,w_j^{(n)}] 
# jest (dla ustalonego j) próbą rozmiaru n ze standardowego rozkładu normalnego.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
BrownianMotion=function(N, n, t, mu, sigma, S){
  macierz <- matrix(0, nrow=N, ncol=n)
  macierz[,1]<-S
  ##S_j(k*t/n) = S_j((k-1)*t/n) + mu * t/n + sigma * \sqrt{t/n} * w_j^(kolumna)
  ## pozostala czesc macierzy
  for (wiersz in 1:N) {
     S_wczesniejsza <-S
    for (kolumna in 2:n) {
      w <- rnorm(1)   # jest (dla ustalonego j) próbą rozmiaru n ze standardowego rozkładu normalnego.
      macierz[wiersz, kolumna] <- S_wczesniejsza + 
                                  mu * t / n + sigma * sqrt(t / n) * w
      S_wczesniejsza<- macierz[wiersz, kolumna]}}
  
  return(macierz)
}
BrownianMotion(N=50, n=7, t=10, mu=0.5, sigma=1, S=50)
#_______________________________________________________________________________
#
#------------------------------Zadanie 7.9--------------------------------------
#
#_______________________________________________________________________________
# Napisz funkcję
#
# MC.euro(TypeFlag=c("call","put"),t,r,sigma,S,X,N=500,n=100)
#
# wyznaczającą przybliżoną cenę europejskiej opcji kupna ("call") 
# lub sprzedaży("put") w chwili 0 za pomocą metody Monte Carlo.
# Odpowiednie symulacje Monte Carlo wykonywane są na podstawie N 
# (domyślnie np. N=500) trajektorii geometrycznego ruchu Browna, 
# startującego z punktu S i określonego na przedziale [0,time].
#
# Każda z N trajektorii jest wygenerowana na podstawie n (domyślnie np. n=100)
# wartości wspomnianego geometrycznego ruchu Browna.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


MC.euro=function(TypeFlag=c("call","put"),t,r,sigma,S,X,N=500,n=100){
  
  
  
  
  
  
  
  
  
}





