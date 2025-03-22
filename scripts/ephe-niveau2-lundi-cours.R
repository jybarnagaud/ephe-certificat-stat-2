setwd("D:/AAAAAAAAAAAAA/certificat/2023/materiel_cours/niveau2/supports/donnees")
library(ade4)
library(vegan)

#--------------------------------#
#### Indices de dissimilarité ####
#--------------------------------#

## Calculer l'indice de Jaccard sous R
#-------------------------------------

# données
libL=matrix(c(0,1,1,0,1,0,1,1,0,0,0,1),byrow=F,ncol=3)

colnames(libL)=c("Libellula depressa","Calopteryx splendens","Calopteryx virgo")
rownames(libL)=c("Site 1","Site 2","Site 3","Site 4")

# distance de Jaccard avec ade4
library(vegan)
jaccard.sim=vegdist(libL,method="jaccard",binary=T)
jaccard.sim

# dendrogramme sur distance de Jaccard
arbre=hclust(jaccard.sim)
plot(arbre,main="Ressemblance entre communautés de libellules",xlab="",ylab="hauteur")

#-----------------------------------------------#
#### Analyse factorielle des correspondances ####
#-----------------------------------------------#

# analyse des correspondances sur une matrice sites x espéces (reptiles)
msp=read.table("LR_reptiles_sites_especes.txt",header=T,sep="\t")
msp[msp>0]=1

head(msp)

# matrice de dissimilarité du Chi²
rept.diss = vegdist(msp,method="chisq")
heatmap(as.matrix(rept.diss))

# calcul des vecteurs et valeurs propres
co.rept=dudi.coa(msp,scannf=F,nf=2)

# valeurs propres = "éboulis des variances"
p1 = barplot(co.rept$eig)
barplot(co.rept$eig,xlab="Rang",ylab="Valeur propre",col="black")
axis(side=1,at=p1,labels=1:23)

# vecteurs propres - ici les deux premiers
s.label(co.rept$li)

# l'éboulis des variances avec les deux premières valeurs propres colorées
barplot(co.rept$eig,xlab="Rang",ylab="Valeur propre",col=rep(c("black","grey"),times=c(2,21)))
axis(side=1,at=p1,labels=1:23)

# même chose : 
screeplot(co.rept)

# les vecteurs propres 2 et 3 : 
co.rept=dudi.coa(msp,scannf=F,nf=3)
s.label(co.rept$li,xax=2,yax=3)

# l'éboulis des variances avec les axes correspondants colorés
barplot(co.rept$eig,xlab="Rang",ylab="Valeur propre",col=rep(c("grey","black","grey"),times=c(1,2,20)))
axis(side=1,at=p1,labels=1:23)

# interprétation du plan 1,2 de l'AFC
s.label(co.rept$li)

# on veut savoir dans quels sites se trouve PODMUR (lézard des murailles)
podmur.loc = factor(msp$PODMUR)
s.class(co.rept$li,fac=podmur.loc,col=c("darkred","steelblue"),
        axesell=F,cellipse=0,cstar=0,cpoint=2)

malmon.loc = factor(msp$MALMON)
s.class(co.rept$li,fac=malmon.loc,col=c("darkred","steelblue"),
        axesell=F,cellipse=0,cstar=0,cpoint=2)

vipber.loc = factor(msp$VIPBER)
s.class(co.rept$li,fac=vipber.loc,col=c("darkred","steelblue"),
        axesell=F,cellipse=0,cstar=0,cpoint=2)

scatter(co.rept)
scatter(co.rept,xax=2,yax=3)

#--------------------------------------------#
#### L'Analyse en Composantes Principales ####
#--------------------------------------------#


# charger les données
clim.lr=read.table("climat_LR_ACP.txt",header=T,sep="\t",row.names=1)

## Le principe de l'ACP
#----------------------

# centrer-réduire les variables
clim.lr.sc=apply(clim.lr,2,scale)

# matrice des corrélations
cor.clim=cor(scale(clim.lr))
library(corrplot)
corrplot(cor.clim)
round(cor.clim,2)

## L'ACP sous R
#--------------

# faire l'acp
library(ade4)
library(adegraphics)
pc=dudi.pca(clim.lr) # NB: sélectionner 2 composantes + entrée
coords=pc$li
head(coords)

# éboulis des variances (voir aussi powerpoint pour figure modifiée)
barplot(pc$eig,col="gray70",xlab="Composantes principales",ylab="Valeurs propres",names.arg= c(1:9))

# éboulis des variances avec % variance expliquée
in.expl=pc$eig/sum(pc$eig)*100
in.expl

cum.in.expl=100*cumsum(pc$eig)/sum(pc$eig)
cum.in.expl

barplot(pc$eig,col="gray70",xlab="Composantes principales",ylab="Valeurs propres",names.arg= c(1:9))
text(x=c(1:9),y=6,round(in.expl,1))
text(x=c(1:9),y=5,round(cum.in.expl,1),font=3)

# quelques exemples fictifs d'éboulis de variance plus ou moins faciles (modifié dans ppt)

par(mfrow=c(3,2))
x3=c(0.50,0.30,0.10,0.05,0.025,0.015,0.005,0.004,0.001)
barplot(x3,xlab="Composantes principales",ylab="Valeurs propres",names.arg= c(1:length(x3)))

x1=c(0.90,0.08,0.009,0.008,0.002,0.0005,0.00025,0.00025)
barplot(x1,xlab="Composantes principales",ylab="Valeurs propres",names.arg= c(1:length(x1)))

x6=c(0.30,0.25,0.20,0.10,0.09,0.05,0.005,0.004,0.001)
barplot(x6,xlab="Composantes principales",ylab="Valeurs propres",names.arg= c(1:length(x6)))

x2=c(0.45,0.30,0.15,0.05,0.025,0.0125,0.0124,0.0001)
barplot(x2,xlab="Composantes principales",ylab="Valeurs propres",names.arg= c(1:length(x2)))

x4=c(0.30,0.25,0.24,0.16,0.05)
barplot(x4,xlab="Composantes principales",ylab="Valeurs propres",names.arg= c(1:length(x4)))

x5=c(0.30,0.25,0.20,0.15,0.10)
barplot(x5,xlab="Composantes principales",ylab="Valeurs propres",names.arg= c(1:length(x5)))


# cercle des corrélations
s.corcircle(pc$co,clabel=1.5,box=F)

# inertie
inertia.dudi(pc, col.inertia=TRUE)


## Projection des individus dans une ACP
#---------------------------------------

s.label(pc$li)

# sortie numérique
coord.var=pc$co
coord.var

#------------------------------------------#
#### Analyse en coordonnées principales ####
#------------------------------------------#


## Les données : analyser un échantillonnage de communautés 
# par ADN environnemental
#-----------------------------------------------------------

library(vegan)

# données
amphi=read.table("amphibiens_edna.txt",header=T,sep="\t",row.names=1)
head(amphi)

## L'analyse en coordonnées principales : 
# l'ordination d'une matrice de dissimilarité
#--------------------------------------------

# calcul de distance
amphi.dist=vegdist(amphi,method="bray")

# ordination
amphi.pco=dudi.pco(amphi.dist) # sélectionner 2 et entrée 

# Le problème des valeurs propres négatives (on ne garde que 2 composantes pour simplifier mais il en faudrait plutôt 3-4 )
amphi.pco2=dudi.pco(sqrt(amphi.dist)) # 2 et entrée

## Représenter les communautés dans l'espace de la PCoA
#------------------------------------------------------

# représentation des individus
s.label(amphi.pco2$li)

## Projeter les espèces sur la PCoA : 
# l'idée de variables supplémentaires
#-----------------------------------

# projeter les variables supplémentaires
sp=supcol(amphi.pco2,amphi)
s.arrow(sp$cosup)

#--------------------------------#
##### L'analyse discriminante ####
#--------------------------------#

# données
petrel.acou=read.table("petrel_acoustique.txt",header=T,sep="\t")
head(petrel.acou)

rownames(petrel.acou)=petrel.acou$PhCode

## Première étape : 
# construire une ordination sur notre table de variables
#-------------------------------------------------------

pc.petrel=dudi.pca(petrel.acou[,-c(1:3)]) # sélectionner 3 axes et entrée

# cercles de corrélation 
x11()
s.corcircle(pc.petrel$co,xax=1,yax=2,sub="PC1 vs PC2")
x11()
s.corcircle(pc.petrel$co,xax=1,yax=3,sub="PC1 vs PC3")
x11()
s.corcircle(pc.petrel$co,xax=2,yax=3,sub="PC2 vs PC3")

## Deuxième étape : différencier graphiquement des
# groupes dans l'espace d'ordination de l'ACP
#-------------------------------------------------

s.label(pc.petrel$li,clabel=0,pch=21)
s.class(pc.petrel$li,petrel.acou$Species,col=c("darkblue","darkred"))

## Troisième étape : étudier les différences 
# entre les groupes : l'analyse discriminante
#--------------------------------------------

discr.petrel=discrimin(pc.petrel,petrel.acou[,c("Species")],scannf=F) # 2 groupes donc 1 seul axe discriminant

# contribution des variables
sum(discr.petrel$eig)/sum(pc.petrel$eig)

## Tester la robustesse de l'analyse 
# discriminante au moyen d'un test par permutations
#---------------------------------------------------

test=randtest(discr.petrel)
plot(test)

## Etudier l'impact du facteur de groupement 
# sur les variables de l'ordination
#-------------------------------------------

plot(discr.petrel)

#----------------------------#
#### Analyse de coinertie ####
#----------------------------#

# analyse des correspondances sur une matrice sites x espéces (reptiles)
msp=read.table("LR_reptiles_sites_especes.txt",header=T,sep="\t")
msp[msp>0]=1
co.rept=dudi.coa(msp,scannf=F,nf=5)

# réduire les données climatiques aux cellules avec données de reptiles
clim.xy=read.table("LR_coordonnees_cellules.txt",header=T,sep="\t")
clim.lr0=read.table("LR_climat_topo.txt",header=T,sep="\t")
# données climatiques sur le Languedoc Roussillon (rasters de Météo France)
clim.lr0=read.table("LR_climat_topo.txt",header=T,sep="\t")
rownames(clim.lr0)=paste("ID",clim.lr0[,1],sep="")

pc.clim2=merge(clim.lr0,clim.xy,by.x=0,by.y="ID",all=F)
rownames(pc.clim2)=pc.clim2$Row.names
pc.clim2=pc.clim2[,-c(1,2)]


# analyse en composantes principales avec pondération des lignes par les poids des lignes
#  dans l'analyse en composantes principales
pc.pcclim2=dudi.pca(pc.clim2[,1:9],row.w = co.rept$lw,scannf=F,nf=3)

# analyse de coinertie: co-structure entre l'acp sur les var climatiques et la cca sur les reptiles
acoi=coinertia(pc.pcclim2,co.rept,scannf=F,nf=2)

# quelques sorties graphiques
s.corcircle(acoi$aX)
s.corcircle(pc.pcclim2$co)
s.corcircle(pc.pcclim2$co,xax=2,yax=3)
s.corcircle(acoi$aY)

# projections des individus
s.label(co.rept$co)

# plot complet
plot(acoi)

# test de robustesse (par permutations). On indique que la table "reptiles" est la table dont
# les poids sont fixés
rd <- randtest(acoi,fixed=2) 

# quelques graphiques supplémentaires
library(adegraphics)
g1 <- s.arrow(acoi$l1, plab.cex = 0.7)
g2 <- s.arrow(acoi$c1, plab.cex = 0.7)
g3 <- s.corcircle(acoi$aX, plot = FALSE)
g4 <- s.corcircle(acoi$aY, plot = FALSE)
cbindADEg(g3, g4, plot = TRUE)
g5 <- plot(acoi)


#------------------------------------------------#
#### Analyse de redondance (données Poissons) ####
#------------------------------------------------#

# ACP sur les variables morphologiques des poissons
pois=read.table("poisson_RDA.txt",sep="\t",header=T)
pois$Age=factor(pois$Age)
pois$Pop=factor(pois$Pop)

pc.pois=dudi.pca(pois[,c(1:3,7)],scannf=F,nf=3)
# sorties graphiques
s.corcircle(pc.pois$co)
s.label(pc.pois$li)
# différencier les sexes et les pops
s.class(pc.pois$li,fac=pois$sex,col=c("red","blue"))
s.class(pc.pois$li,fac=pois$Pop,col=1:8)

# analyse de redondances
rda.pois=pcaiv(pc.pois,pois[,c("Pop","sex")])
summary(rda.pois)

# quelques sorties graphiques
s.label(rda.pois$fa)
s.arrow(rda.pois$co)
s.label(rda.pois$cor)
s.corcircle(rda.pois$as)

s.label(rda.pois$l1, 2, 1, clab = 0, cpoi = 1.5)
s.label(rda.pois$co, 2, 1, add.plot = TRUE)
