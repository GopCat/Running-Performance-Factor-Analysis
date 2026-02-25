# Running Performance Factor Analysis
# Author: Aleksandra Ivanova
# Statistical dimensionality reduction project

install.packages("readxl")
install.packages("MASS")
install.packages("Hmisc")
install.packages("ggfortify")
install.packages("corrplot")
install.packages("REdaS")
install.packages("NbClust")
install.packages("fpc")
install.packages("factoextra")
install.packages("clValid")
install.packages("openxlsx")
install.packages("psych")
install.packages("rattle")
install.packages("GPArotation")
install.packages("polycor")

library(readxl)
library(DataExplorer)
library(MASS)
library(Hmisc)
library(ggfortify)
library(corrplot)
library(REdaS)
library(NbClust)
library(fpc)
library("FactoMineR")
library(clValid)
library(openxlsx) 
library(GPArotation)
library(polycor)


library(psych)
library(rattle)
library(factoextra)

# Import runners.csv
wtr <- runners
head(wtr)

introduce(wtr)
sapply(wtr,class) # typy promennych
sum(duplicated(wtr)) 
colSums(is.na(wtr))

# pouze kvantitativni data
wtr_x<-wtr[,2:8]
head(wtr_x)
m_tr=ncol(wtr_x) # pocet sloupcu
m_tr 

# Mahalabisova vzdalenost
means <- sapply(wtr_x, mean, na.rm = TRUE)
mahal_wtr <- mahalanobis(wtr_x, means, cov(wtr_x), inverted = FALSE)      
summary(mahal_wtr)
plot(mahal_wtr)
plot(density(mahal_wtr, bw = 0.3), main="Squared Mahalanobis distances") ; rug(mahal_wtr)
alfa = 0.001
out_wtr <- which(mahal_wtr > qchisq(1-alfa,df=ncol(wtr_x)))
out_wtr

# vylouceni odlehlych pozorovani 			
wtr2 <- wtr[-c(which(mahal_wtr > qchisq(1-alfa,df=ncol(wtr_x)))),]
wtr_x2 <- wtr[-c(which(mahal_wtr > qchisq(1-alfa,df=ncol(wtr_x)))),2:8]

##### overeni faktorovatelnosti #####

#a) korelacni matice
corrplot(cor(wtr_x2), "number")
det(cor(wtr_x2)) 
cor_wtr_inv <- solve(cor(wtr_x2)) 
VIF_wtr <- diag(cor_wtr_inv) 
VIF_wtr > 5	
#skoro vsechny korelacni koeficienty jsou vyssi nez 0,5 => lze provest FA

#b) Anti-image correlation matrix diagonals
anti_image_wtr <- diag(1, nrow = ncol(wtr_x2)) 
- diag(1 / diag(cor_wtr_inv)) %*% cor_wtr_inv %*% diag(1 / diag(cor_wtr_inv))
round(anti_image_wtr, 3)
print(round(diag(anti_image_wtr),3), row.names = FALSE)
#diagonalni prvky anti-image korelacni matice by mely byt vetsi nez 0,5 -> tyto diagonalni prvky jsou rovny 1 => lze provest FA

#c) Kaiser Meyer Olkin test
KMOS(wtr_x2, use = "complete.obs")
#KMO test je vetsi nez 0,5, vychazi cca 0,84 => lze provest FA

#d) Bartletuv test sfericity
n=nrow(wtr_x2)
cortest.bartlett(cor(wtr_x2),n)
#p-value je vyznamne nizsi nez 0,05 => korelacni matice neni jednotkova matice, existuje vztah mezi promennymi.

#e) rozsah datoveho souboru - melo by byt minimalne 5x vice objektu nez znaku
k=ncol(wtr_x2)
Pomer<- n/k
Pomer # ~7,3x; dostatecne velky rozsah datoveho souboru pro ucely FA

##### optimalni pocet faktoru #####

# kovariancní a korelacni matice
cov_wtr<-cov(wtr_x2)
cor_wtr<-cor(wtr_x2)

# vlastni cisla a vlastni vektory
vlc_cor_wtr <- eigen(cor_wtr)
vlc_cor_wtr$values
vlc_cor_wtr$vectors
# Kaiserovo kriterium: zachovat factory s vl.cisly >1 (v nasem pripade bychom meli pouze 1 factor... pouzijeme 2 nebo 3)

# analyza hlavnich komponent
# optimalni pocet hlavnich komponent
sum(vlc_cov_wtr$values)/m_tr	#minimalni "akceptovatelna" hodnota vlastniho cisla pro hlavni komponenty, v pripade analyzy zalozene na kovariancni matici
D_e_tr = 1-(det(cor(wtr_x2))^(1/m_tr))  # efektivni rozptyl
k_tr = m_tr*(0.8-0.5*D_e_tr) 
k_tr # pocet vyznamnych hlavnich komponent

pca_wtr <- prcomp(wtr_x2, center = TRUE,scale = TRUE) #scale=TRUE -> nejsou stejne merne jednotky, nutno znormovat data
summary(pca_wtr)
pca_wtr$rotation #komponentni zateze/korelace (loadings)
pca_wtr$x #komponentni skore

# grafy
plot(pca_wtr,type="lines")
screeplot(pca_wtr) 
fa.parallel(wtr_x2, fa = "fa", n.iter = 100, show.legend = TRUE) # porovnava skutecna vlastni cisla s temi z nahodnych dat. 
# pocet faktoru nebo komponent je urcen jako pocet vlastnich cisel skutecných dat, ktera jsou vetsi nez odpovidajici vlastni cisla nahodnych dat.
autoplot(pca_wtr)
biplot(pca_wtr) # spodni osa(komponentni skore PC1), leva osa (komponentni skore PC2),
# horni osa (komponentní zatez PC1), prava osa (komponentni zatez PC2).


##### faktorova analyza 2 faktory #####
fanal0_2 <- factanal(wtr_x2, 2, rotation="none",scores = "Bartlett")
fanal_2 <- factanal(wtr_x2, 2, rotation="varimax",scores = "Bartlett")

fanal_2$loadings

round(cbind(fanal0_2$loadings[,1:2],fanal_2$loadings[,1:2]),2)
round(cbind(fanal0_2$loadings[,1:2]-fanal_2$loadings[,1:2]),2)

fanal_2$uniquenesses #zbytek rozptylu promenne, ktera neni vysvetlena vybranymi faktory
100*(1-fanal_2$uniquenesses) #komunality (vysvetlena cast rozptylu vybranymi faktory)
#komunalita i-te promenne vyjadruje miru promenlivosti a je vahou, s jakou jednotlive spolecne faktory
#prispivaji do rozptylu dane promenne.

cor(fanal_2$scores) 
fanal_2$rotmat #transformacni matice
fanal_2$PVAL # H0 zvoleny pocet faktoru je dostatecny pro zachyceni dimenze ulohy (H0 zamitame pro PVAL<0,05)
# PVAl=0,048 => malo faktoru, mame zvetsit pocet
print(fanal_2, digits=3, cutoff=.3, sort=TRUE)

##### faktorova analyza 3 faktory #####
fanal0 <- factanal(wtr_x2, 3, rotation="none",scores = "Bartlett")
fanalv <- factanal(wtr_x2, 3, rotation="varimax",scores = "Bartlett")
fanalo <- factanal(wtr_x2, 3, rotation="oblimin",scores = "Bartlett")

# Varimax rotace
fanalv$loadings

round(cbind(fanal0$loadings[,1:3],fanalv$loadings[,1:3]),3)
round(cbind(fanal0$loadings[,1:3]-fanalv$loadings[,1:3]),3)

fanalv$uniquenesses #zbytek rozptylu promenne, ktera neni vysvetlena vybranymi faktory
100*(1-fanalv$uniquenesses)

cor(fanalv$scores) 
fanalv$rotmat #transformacni matice
fanalv$PVAL # H0 zvoleny pocet faktoru je dostatecny pro zachyceni dimenze ulohy (H0 zamitame pro PVAL<0,05)
# PVAL = 0,447 => pocet faktoru OK
print(fanalv, digits=3, cutoff=.3, sort=TRUE)

# plot factor 1 by factor 2
loadv <- fanalv$loadings[,1:2]
plot(loadv,type="n", main = "Varimax rotace") # set up plot
text(loadv,labels=names(wtr_x2),cex=.7) # add variable names


##### faktorova analyza 3 faktory (2. postup) #####
# bez rotace
fa_0 <- fa(wtr_x2, nfactors = 3, rotate = "none", fm = "ml")
print(fa_0)

# varimax rotace
fa_var <- fa(wtr_x2, nfactors = 3, rotate = "varimax", fm = "ml")
print(fa_var)

# oblimin rotace
fa_obl <- fa(wtr_x2, nfactors = 3, rotate = "oblimin", fm = "ml")
print(fa_obl)

# Vizualizace faktoru
fa.diagram(fa_var, main = "Varimax Rotace")
fa.diagram(fa_obl, main = "Oblimin Rotace")

# Interpretace faktoru: Faktory
print(fa_var$loadings)
print(fa_obl$loadings)
