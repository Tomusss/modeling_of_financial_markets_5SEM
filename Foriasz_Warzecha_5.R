
#-----------------*************** Lista 5 MRF ***************-------------------


#_______________________________________________________________________________
#
#------------------------------Zadanie 5.1--------------------------------------
#
#_______________________________________________________________________________
# Na rynku występują jedynie akcje dwóch spółek Alfa i Omega.
# Cena jednej akcji spółki Alfa wynosi 30 PLN. 
# Liczba akcji spółki Alfa równa jest 10 000. 
# Współczynnik beta akcji spółki Alfa wynosi 1,1. 
# Cena jednej akcji spółki Omega wynosi 10 PLN. 
# Liczba akcji spółki Omega jest równa 20 000. 
# Wyznacz wartość współczynnika beta dla akcji spółki Omega.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

cena_akcji_a = 30
liczba_a = 10000
beta_a = 1.1
cena_akcji_o=10
liczba_o=20000

Va <- cena_akcji_a*liczba_a
Vo <- liczba_o*cena_akcji_o
V <- Va+Vo

wa <- Va/V
wo<- Vo/V

#beta = beta_a*wa + beta_o*wo, beta jeden bo wzor Cov(K,K) = Var(K)/Var(K)=1

beta_o <- (1-(beta_a*wa))/wo
beta_o

#_______________________________________________________________________________
#
#------------------------------Zadanie 5.2--------------------------------------
#
#_______________________________________________________________________________
# Załóżmy, że współczynnik korelacji pomiędzy stopą zwrotu z akcji A,
# a stopą zwrotu z portfela rynkowego wynosi 0.9.
# Jaką część całkowitego ryzyka stopy zwrotu z akcji A
# stanowi tutaj ryzyko niesystematyczne?
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# ryzyko systematyczne - odnosi się do całego rynku
# ryzyko niesystematyczne - odnosi się do akcji A

# ryzyko = ryzyko systematyczne + ryzyko niesystematyczne
# sigma_v^2 = beta_v^2*sigma_m^2 + Var(E_v)

# calkowite ryzyko oznacza 100% wariancji

# czesc systematycznego = beta_v^2*sigma_m^2/sigma_v^2 = wspolczynnik^2

wspolczynnik = 0.9
wspolczynnik2 = 0.9^2

# cala czesc - 1

czesc_niesystematyczna <- 1 - wspolczynnik2
czesc_niesystematyczna 


#_______________________________________________________________________________
#
#------------------------------Zadanie 5.3--------------------------------------
#
#_______________________________________________________________________________
# Oczekiwana stopa zwrotu z portfela rynkowego wynosi 18% a wariancja 
# tej stopy zwrotu jest równa 7%. Oczekiwana stopa zwrotu z akcji 
# spółki A wynosi 24%. Stopa zwrotu z aktywów wolnych od ryzyka to 4%. 
# Wyznacz kowariancję między stopą zwrotu z portfela rynkowego i stopa zwrotu 
# z portfela akcji spółki A.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

mi_a = 0.18
wariancja_mi_a = 0.07
mi_v = 0.24
r_f = 0.04

beta = (mi_v - r_f)/(mi_a - r_f)
kowariancja = beta* wariancja_mi_v
kowariancja

#_______________________________________________________________________________
#
#------------------------------Zadanie 5.4--------------------------------------
#
#_______________________________________________________________________________
# Portfel inwestycyjny składa się z 5 akcji (od A1 do A5). 
# Współczynnik beta dla portfela wynosi 1,7. 
# Każda z akcji ma równy wartościowy udział w portfelu.
# Współczynnik beta dla akcji A1 wynosi 1,8. 
# Akcja A1 zostałą sprzedana i w jej miejsce zakupiono akcję X. 
# Ile powinien wynieść współczynnik beta dla nowo zakupionej akcji X, 
# aby współczynnik beta dla portfela osiągnął poziom 1,65?
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 








#_______________________________________________________________________________
#
#------------------------------Zadanie 5.5--------------------------------------
#
#_______________________________________________________________________________
# Współczynnik beta akcji spółki ETA wynosi 1,2, a oczekiwana stopa zwrotu 
# z akcji tej spółki to 8%. Oczekiwana stopa zwrotu z portfela rynkowego 
# wynosi 7%. Korzystając z CAPM wyznacz wartość stopy wolnej od ryzyka. 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

beta = 1.2
mi_v = 0.08
mi_m = 0.07
r_f = (mi_v-beta*mi_m)/(1-beta)
r_f


#_______________________________________________________________________________
#
#------------------------------Zadanie 5.6--------------------------------------
#
#_______________________________________________________________________________
# Kurs akcji A w dniu 1 stycznia 2023 to 120. Posiadacz akcji A 
# otrzymuje dywidendę wypłacana 1 lipca 2023 w wysokości 1 (na każdą akcję) 
# i ponownie w dniu 1 października 2023 w wysokości 2. 
# Stopa procentowa w roku 2023 wynosi 12% w skali roku dla każdego okresu
# (kapitalizacja ciągła). Cena akcji A w kontrakcie terminowym “forward” 
# z dostawą w dniu 1 listopada 2023 to 131. W obliczeniach załóż, że każdy 
# miesiąc to dokładnie 1/12 roku.
#
# (a) Uzasadnij, że w dniu 1 stycznia 2023 istniała możliwość arbitrażu.
#
# (b) Oblicz zysk z tego arbitrażu na dzień 1 listopada 2023. 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# podpunkt 5.5(a)





# podpunkt 5.5(b)





#_______________________________________________________________________________
#
#------------------------------Zadanie 5.7--------------------------------------
#
#_______________________________________________________________________________
# Jest 01.12.2024. Kurs akcji S wynosi 400. Wiadomo, że na każdą akcję S w
# dniu 1.01.2025 będzie wypłacana dywidenda w wysokości 15. Stopa wolna od
# ryzyka (kapitalizacja ciągła) wynosi 6% dla wszystkich okresów.
#
# (a) Jaka jest teoretyczna cena jednorocznego kontraktu forward na tę akcję?
# 
# (b) Pośrednik kwotuje kontakty na forward na tę akcję po 420
#     (prowizja wynosi 1 za kontrakt, płatne z góry). Jesteś arbitrażystą. 
#     Możesz się finansować po stopie 6,25% w skali roku i lokować wolne środki
#     finansowe wg stopy 5,75% w skali roku. Jak wygląda Twoja strategia?
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# podpunkt 5.6(a)





# podpunkt 5.6(b)




#_______________________________________________________________________________
#
#------------------------------Zadanie 5.8--------------------------------------
#
#_______________________________________________________________________________
# Cena spot srebra wynosi 48 USD za uncję. 
# Koszty przechowywania wynoszą 0.64 USD za uncję na rok 
# płatne kwartalnie z góry. Zakładając, że stopa procentowa w USD wynosi 
# 4.5% na rok (kapitalizacja ciągła) dla wszystkich okresów oblicz cenę srebra 
# w kontraktach futures z dostawą za 9 miesięcy.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 







#_______________________________________________________________________________
#
#------------------------------Zadanie 5.9--------------------------------------
#
#_______________________________________________________________________________
# Czy dla podanych kursów EUR/PLN możliwy jest arbitraż? 
# Jeśli tak, to podaj zysk z niego przy inwestycji 1 000 000 PLN.
#
#             Kupno (Bid)   Sprzedaż (Ask) 
#
#   Rynek 1     4.3157        4.3227 
#
#   Rynek 2     4.3257        4.3360 
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# sprzedaż - kupno > 0
# my kupujemy gdy oni sprzedają i my sprzedajemy gdy oni kupują

zysk =  max(4.3157,4.3257) - min(4.3227, 4.3360)
zysk

# zatem możliwy jest arbitraż

inwestycja = 1000000
zysk_i = (inwestycja/min(4.3227, 4.3360) )*max(4.3157,4.3257)
zysk = zysk_i - inwestycja
zysk


#_______________________________________________________________________________
#
#------------------------------Zadanie 5.10--------------------------------------
#
#_______________________________________________________________________________
# Wiadomo, że kursy walut wynoszą
#
#   USD/EUR = 0.9405,
#   USD/JPY = 157.771,
#   EUR/JPY = 168.708.
#
# Czy jest tutaj możliwy arbitraż dla osoby posiadającej dolary? 
# Jeśli tak, to jaki byłby z niego zysk przy inwestycji 100 000 USD?
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# 1 opcja: kupujemy EUR, potem JPY
# 2 opcja: kupujemy JPY



