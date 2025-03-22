#---------------------------------------------------------------------------#
### Formation "Analyse de données pour les écologues renforcement" - EPHE ###
### Jean-Yves Barnagaud : jean-yves.barnagaud@ephe.psl.eu ###
# TD du lundi : exploration de données par les méthodes multivariées
#---------------------------------------------------------------------------#

# Note : seuls les deux premiers cas d'étude sont entièrement détaillés. Pour les autres, 
# seules les lignes de code sont données car les méthodes employées sont les mêmes.

library(ggplot2)
library(cowplot)
library(ade4)
library(factoextra)
library(sf)
library(mapview)
webshot::install_phantomjs()  # pour les cartes

rm(list=ls())
setwd('D:/AAAAAAAAAAAAA/certificat/certificat_2023/niveau2/supports/donnees')

#-------------------------------------------------------------------#
#### Cas d'étude n°1 : communautés d'oiseaux en Nouvelle Zélande ####
#-------------------------------------------------------------------#

## données

sp <- read.csv2("especes_NZ.csv")
com <- read.csv2("communautes_oiseaux_NZ.csv")
habitats <- read.csv2("habitats_NZ.csv")

## analyse de la matrice sites*espèces

co.sp <- dudi.coa(com) # je choisis de garder 3 axes
co.sp <- dudi.coa(com,nf=3,scannf=F) # nb : les axes étant orthogonaux, le choix du nombre d'axes n'a d'impact que sur la quantité de données à stocker dans la mémoire R, pas sur les  coordonnées de l'analyse.

## plot des sites

s.label(co.sp$li) # on voit d'emblée que le site 1241 est très éloigné des autres et va fortement tirer l'analyse, il faut l'explorer
s.arrow(co.sp$co)
com["1241",] # il s'agit du seul site avec le Pluvier roux - il semble plus pertinent d'éliminer
# cette espèce afin d'éviter qu'elle ne tire exagérément l'analyse (mais !! être transparent sur cette
# manipulation et sa justification dans les méthodes). 

com.nodot <- com[,-which(colnames(com)=="DOT")]
co.sp2 <- dudi.coa(com.nodot,nf=2,scannf=F)
fviz_ca(co.sp2) # on n'y voit pas grand chose, il vaut mieux séparer les graphiques

s.label(co.sp2$li) # on voit un étalement de quelques sites sur l'axe 2, et un étalement de l'ensemble sur l'axe 1
s.arrow(co.sp2$co) # cet étalement s'explique par MAL (Canard colvert), ce qui nous oriente vers des sites aquatiques
# pour l'axe 1, on ne voit toujours rien de bien net car le nuage de point est écrasé. On oppose ROB (Miro rubisole, un endémique forestier)
# à BCT (Guifette des galets, oiseau d'eau) et globalement, l'ensemble nous suggère un gradient de la forêt vers les milieux ouverts. 

## positionnement des espèces introduites

sp.nodot <- sp[-which(rownames(sp)=="DOT"),]
s.class(co.sp2$co,fac=factor(sp.nodot$Status))

## espèces introduites seules
sp.intro  <- subset(sp.nodot,Status=="Introduced")
com.intro <- com.nodot[,rownames(sp.intro)]
co.sp3 <- dudi.coa(com.intro,scannf,nf=2)

s.label(co.sp3$li)
s.label(co.sp3$co) # on explore l'AFC des espèces introduites

s.class(co.sp3$co,fac=factor(sp.intro$Reason)) # raison de l'introduction
s.class(co.sp3$co,fac=factor(sp.intro$Realm_biome)) # raison de l'introduction

# Conclusion de cette première analyse : il n'y a pas de structuration liée aux espèces introduites dans les communautés étudiées

## relation à l'habitat

summary(apply(habitats[,-c(8,9,10)],1,"sum")) # on vérifie que les différentes classes d'habitat ne somment pas à 100% (sinon il faudrait en tenir compte via d'autres méthodes)

pc.hab <- dudi.pca(habitats) # 3 axes clairs
s.corcircle(pc.hab$co) # axe 1 sépare les milieux natifs en montagne, dominés par de la forêt, des milieux ouverts de plaine, dominés par des prairies exotiques

# l'axe 2 sépare des paysages homogènes et des paysages hétérogènes dominés par des plantations 

s.corcircle(pc.hab$co,xax=1,yax=3)
s.corcircle(pc.hab$co,xax=2,yax=3) # l'axe 3 sépare des milieux de buissons de plaine / moyenne altitude vs la forêt

s.label(pc.hab$li)
s.label(pc.hab$li,xax=1,yax=3) # on voit un léger effet d'arche
s.label(pc.hab$li,xax=2,yax=3)

## coinertie

coi.nz  <- coinertia(co.sp2,pc.hab) # on ne peut pas faire la coinertie parce que les poids des lignes ne sont pas homogènes : il faut donc recalculer l'ACP en lui imposant les poids de lignes de l'AFC
pc.habw <- dudi.pca(habitats,row.w=co.sp2$lw)
coi.nz <- coinertia(co.sp2,pc.habw) # classique de la coinertie, un axe nettement majoritaire - ici on en garde 2
plot(coi.nz)

# test par permutations de la robustesse du lien entre les tables

test.coi=RV.rtest(pc.habw$tab, co.sp2$tab, 99)
plot(test.coi) # le coefficient RV évalue la similitude entre deux tables. Ici, il est bien plus élevé que l'attendu des tests de permutation

## interprétation des sorties

s.arrow(coi.nz$l1, clab = 0.7) # on distingue sur l'axe 1 d'une part les espèces associées aux prairies exotiques et de l'autre, un mélange d'espèces associées aux milieux natifs, aux paysages fragmentés, à de la végétation haute. Il faudra la sortie des espèces pour l'interpréter finement
# le deuxième axe sépare les espèces de plantations des espèces d'altitude, il s'agit probablement d'un gradient altitudinal (les plantations sont  à basse et moyenne altitude)

s.label(coi.nz$c1, clab = 0.7) # en bas (négatif axe 2, positif axe 1), on a un groupe d'espèces de forêts natives. En haut (positif axe 2), on trouve la chevêche (LGO)
# et le cacatoès (SUL), isolés du côté des plantations : ces deux espèces ont donc des profils d'habitat particuliers qui tirent cet axe. De l'autre côté, le groupe des 
# espèces natives est logiquement associé aux forêts natives. Au niveau du centroïde, les espèces décalées en bas à droite (YEL, SKY...) sont toutes des espèces de milieux
# ouverts, pour la plupart exotiques - elles sont donc logiquement positionnées près des prairies exotiques. En négatif de l'axe 1 (proche du centroide : MOR, TOM), on a surtout
# des espèces forestières.

s.class(coi.nz$c1,fac=factor(sp.nodot$Status)) # les associations espèces-habitats sont bien distinctes entre les 3 groupes (alors que nous avons précédemment montré
# que les communautés elles-mêmes n'étaient pas distinctes). Les espèces endémiques sont plutôt forestières et d'altitude, les espèces introduites s'étalent plutôt entre 
# des milieux ouverts et de plantation, et les espèces natives sont entre les deux.

# on peut voir la relation entre l'AFC, l'ACP et la coinertie avec des cercles de corrélation : 
s.corcircle(coi.nz$aX)
s.corcircle(coi.nz$aY)

## Conclusion : Les communautés sont mélangées : des espèces natives et exotiques sont simultanément présentes dans les communautés. Néanmoins, il y a 
# un clivage assez net lié à l'habitat : les espèces natives dominent en milieux natifs, donc en altitude et en forêt native, alors que les espèces exotiques 
# sont plus présentes en plantations, en plaine et dans une moindre mesure en milieu ouvert - avec des groupes relativement homogène d'espèces natives ou 
# exotiques partageant des affinités d'habitat proches. Il y a quelques espèces à part : chevêche, cacatoès par exemple. 

#----------------------------------------------------------------------------#
#### Cas d'étude n°2 : traits écologiques des oiseaux de Nouvelle Zélande ####
#----------------------------------------------------------------------------#

## on reprend les données du cas pratique 1 et on ajoute les suivantes :

traits <- read.csv2("traits_NZ.csv")

# deux possibilités pour explorer un jeu de données à multiples types de variables : 
# l'analyse en coordonnées principales : elle se basera alors sur une matrice de distance de Gower
# c'est la méthode la plus flexible, mais avec deux difficultés : il faut pondérer correctement les 
# calculs de dissimilarité (fonction dist.ktab dans ade4, que nous n'utiliserons pas ici) et l'interprétation
# est indirecte
# l'analyse Hill & Smith, qui est un hybride entre une AFC et une ACP : elle est moins flexible et ne permet
# pas de résumer certains types de traits (comme les traits exprimés en proportions, dits "fuzzy"), mais elle 
# a le mérite de la simplicité : nous allons l'utiliser dans cet exemple.

hs.nz <- dudi.hillsmith(traits) # deux décrochements : nous choisissons le premier

## projection des traits

s.arrow(hs.nz$co) 

# l'axe 1 sépare de grosses espèces sociales, granivores, nichant au sol (= probablement espèces de milieux ouverts)
# d'espèces nichant dans les arbres, plutôt insectivores ou nectarivores. Noter que contrairement à une ACP, l'espace n'est pas normé
# (donc la longueur absolue de la flèche ne veut rien dire). L'axe 2 sépare d'un côté les frugivores et les carnivores, de l'autre côté des espèces sédentaires
# nichant dans des cavités d'arbres.

## projection des espèces

s.label(hs.nz$li) 

# référez-vous à l'objet "sp" pour les noms d'espèces complets. On retrouve bien les espèces de milieux ouverts (dont
# oiseaux d'eau) dans le négatif de l'axe 1, les espèces forestières insectivores vers le bas et dans le positif de l'axe 1, un frugivore
# (PIG, pigeon) dans le haut de l'axe 2. On a donc un agencement des traits qui reflète assez bien les préférences d'habitats des espèces.
# On appelle ces axes principaux construits sur des traits des "syndromes de traits", c'est à dire des faisceaux de traits corrélés qui 
# portent une information sur une stratégie écologique - bien sûr, pour s'assurer qu'il y a bien une relation évolutive entre ces syndromes
# et l'habitat, il faudrait aller bien plus loin, mais c'est une première exploration intéressante. 

## projection des guildes native / endemic / exotic

s.class(hs.nz$li,fac=factor(sp$Status))

# on ne voit pas une grosse différence entre les 3 guildes, qui sont légèrement décalées dans l'espace de traits mais se recoupent beaucoup.
# notez que les espèces natives sont entièrement englobées dans les traits des espèces endémiques : la diversité en traits est (paradoxalement) 
# plus grande chez les endémiques et couvre toute la gamme de traits des natives. Les espèces introduites sont légèrement décalées vers les traits
# liés aux milieux ouverts. ATTENTION : ne pas surinterpréter les sorties de s.class ; tenir compte de la manière dont sont calculées les ellipses. On 
# peut d'ailleurs la changer : 

s.class(hs.nz$li,fac=factor(sp$Status),cellipse=0.5)
s.class(hs.nz$li,fac=factor(sp$Status),cellipse=0.7)
s.class(hs.nz$li,fac=factor(sp$Status),cellipse=0.3)
s.class(hs.nz$li,fac=factor(sp$Status),cellipse=0.1)

## analyse discriminante

disc.nz <- discrimin(hs.nz,fac=factor(sp$Status),scannf=F)
disc.nz$eig # % de variance expliqué par les deux axes de l'analyse
disc.nz

plot(randtest(disc.nz)) 

# ici, on teste la relation entre l'ordination et le facteur à partir de permutations. On calcule une statistique, 
# dite "statistique de Pillail", qui résume l'intensité de cette relation à partir d'un calcul sur la somme des valeurs propres. Ici, on 
# voit que la statistique observée est en marge supérieure des permutations, ce qui indique qu'il y a bien une structure de traits entre les
# trois guildes, mais elle reste faible. 

# sortie graphique de l'analyse

plot(disc.nz)

## Pour aller plus loin, il faudrait coupler la table de traits et la table de variables d'habitat : cela permettrait de savoir si certains 
# traits sont associés à certains types d'habitats, reflétant ainsi un filtrage environnemental - c'est bien possible vu que le cas pratique 1
# nous a montré que les espèces étaient ségrégées entre habitats. C'est assez attendu pour une espèce, mais l'impact sur les traits est beaucoup
# moins prévisible. Nous avons fait cette étude avec une matrice de traits un peu augmentée, et nous montrons qu'en effet, il y a bien une ségrégation
# des traits entre habitats, et que cette ségrégation est elle-même associée à l'origine des espèces (native / endémique / exotique) et des habitats
# (natifs ou non) (Barnagaud et al, soumis ; Mossion, 2019). Nous verrons comment faire cette analyse au niveau 3. 

## Le choix des traits a une grosse influence sur l'analyse - et ce, même si un trait ne domine pas l'ensemble. Par exemple, dans notre étude, nous
# avons ajouté quelques traits relatifs à la reproduction et à l'alimentation - cela suffit à faire apparaître des structures claires de relations
# traits-communautés. Cela ne veut pour autant pas simplement dire que ce sont les traits que nous avons ajoutés qui sont responsables de ces 
# structures - c'est l'ensemble des corrélations de l'analyse qui est affecté. Ayez cela en mémoire quand vous analysez des traits.

#-----------------------------------------------------------------#
#### Cas d'étude n°3 : paysages sonores de Polynésie française ####
#-----------------------------------------------------------------#

# données

polynesie.acou <- read.csv2("polynesie_acoustique.csv")

df.pca = polynesie.acou[,c("ACI","BIOAC","H","ADI","NP")]
rownames(df.pca) = polynesie.acou$FILE

# ACP sur les indices acoustiques

pca.acou <- dudi.pca(df.pca,scannf=F,nf=3)

s.corcircle(pca.acou$co)
s.label(pca.acou$li)

# Dans ce script, on va tenir compte du fait qu'il y a plusieurs enregistrements par site en calculant une analyse inter-classe
# (between class analysis). Cette ordination ressemble à une analyse discriminante dans laquelle le site serait le facteur de 
# groupement. Elle maximise la variabilité inter-sites à partir des multiples lignes de chaque site. Il en résulte une ordination
# qui s'interprète comme une ACP, mais avec un seul point par site (et non un point par fichier son, comme dans notre analyse initiale)
# La BCA est une manière de tenir compte des données répétées dans une analyse multivariée. NB : ce n'est pas nécessaire - contrairement 
# aux approches par modèles linéaires qui nécessitent que les individus soient indépendants, les ordinations ne mettent pas de telle condition. 
# On peut donc très bien, si on préfère, s'en tenir à l'ACP précédente. L'intérêt de la BCA est de faciliter l'exploration visuelle (graphique)
# des différences entre sites. A noter qu'il existe une version qui maximise la variabilité intra-facteur de groupement : la within class analysis
# (WCA). Utilisez l'une ou l'autre selon vos besoins. 

bca.acou <- bca(pca.acou,fac=factor(polynesie.acou$sites),scannf=F,nf=2)

s.corcircle(bca.acou$co)
s.label(bca.acou$li)
s.class(bca.acou$ls,fac=factor(polynesie.acou$MODALITE))

# On fait ensuite l'analyse discriminante sur la BCA pour différencier la vallée test et la vallée contrôle

s.class(bca.acou$ls,fac=factor(polynesie.acou$MODALITE))

z <-  unique(polynesie.acou[,c("sites","MODALITE")])
discrim <- factor(z$MODALITE)
disc.acou <- discrimin(bca.acou,fac=discrim,scannf=F,nf=1)
plot(disc.acou)

# NB : contrairement à l'analyse discriminante du cas pratique 2, on ne peut pas faire de test de permutations (non encore implémenté dans ade4 pour une analyse
# discriminante à un seul axe)

# Quelques éléments de conclusion. On détecte effectivement une différence entre vallée test et vallée contrôle qu'on serait tenté d'attribuer aux mesures de conservation.
# En réalité, cette différence est partiellement liée à la présence d'une rivière dans la vallée test, qui modifie fortement le paysage sonore voire le sature
# complètement sur certains enregistreurs. Cet effet confondant était inévitable vu la configuration du terrain, mais il pose un gros problème d'interprétabilité du résultat.
# En 03/2023, la solution (s'il y en a une) est en cours d'investigation.

#---------------------------------------------------------------#
#### Cas d'étude n°4 : communautés de rhopalocères landaises ####
#---------------------------------------------------------------#

# données (nb : noms des espèces dans aquitaine_rhopaloceres_index.csv)

rhopa <- read.csv2("aquitaine_rhopaloceres.csv",row.names=1)
rhopa.hab <- read.csv2("aquitaine_rhopaloceres_habitat.csv",row.names=1)


# On peut  étudier les distributions des espèces sur chacune des variables paysagères: c'est tout simplement le barycentre des comptages de l'espèce sur le gradient voulu, 
# et son écart-type. Prenons par exemple la distribution des espèces sur le gradient de % de strate arborée:

sco.distri(rhopa.hab$cover_strate_arboree,rhopa)

# AFC sur la matrice sites x espèces

coa.rhopa <- dudi.coa(rhopa,scannf=F,nf=2)

scatter(coa.rhopa)
s.label(coa.rhopa$li)
s.arrow(coa.rhopa$co)

# ACP sur la matrice d'habitat

pca.rhopa <- dudi.pca(rhopa.hab,scannf=F,nf=3)

s.corcircle(pca.rhopa$co)
s.label(pca.rhopa$li)
scatter(pca.rhopa)

# On recalcule l'ACP afin de la pondérer par les lignes de l'AFC (indispensable pour l'analyse de coinertie)

pca.rhopa <- dudi.pca(rhopa.hab,scannf=F,nf=3,row.w=coa.rhopa$lw)

# l'analyse de coinertie

coi.rhopa <- coinertia(coa.rhopa,pca.rhopa,scannf=F,nf=2)
scatter(coi.rhopa)
s.label(coi.rhopa$co)
s.label(coi.rhopa$li)

# On voit assez nettement se dessiner des gradients et des groupes d'espèces dans cet espace de coinertie, qui donne l'impression qu'il y a  une structure liée au paysage.
# Testons-la avec un test par permutations : 

plot(randtest(coi.rhopa))

# On voit qu'en réalité, l'association espèces - habitats n'est pas robuste. Cela ne veut pas dire qu'elle n'existe pas : l'analyse de coinertie ne fait que décrire les données
# et par conséquent, le patron observé est bien réel. Cependant, on aurait assez facilement pu le retrouver si les espèces avaient été distribuées aléatoirement sur les gradients
#  paysagers, ce qui suggère que l'association observée n'est peut-être pas assez forte pour en tirer la conclusion que les espèces sont distribuées entre communautés par le paysage.

# On peut aussi passer par une analyse de niche (dite analyse "Outlying Mean Index"), qui par rapport à la coinertie a l'avantage d'être plus flexible quant aux 
# relations espèces habitats (la coinertie implique des relations linéaires somme toute peu probables). Par ailleurs, l'OMI s'accommode mieux des espèces rares.

niche.rhopa <- niche(pca.rhopa,rhopa,scannf=F,nf=2)

scatter(niche.rhopa)

# NB : pour cette analyse, aucun test de permutations n'est implémenté. Il faut donc se contenter d'une interprétation précautionneuse du patron observé. 
# Un intérêt de l'OMI est qu'elle donne accès à quelques quantités intéressantes : 

niche.param(niche.rhopa)

# inertia : la quantité d'inertie portée par l'espèce - c'est à dire à quel point elle s'écarte du centroïde (de l'espèce moyenne fictive)
# OMI : un indice qui exprime la marginalité de l'espèce dans l'espace de niche, c'est à dire à quel point elle se situe à la marge du nuage de points dessiné par l'ensemble des espèces
# Tol : un indice qui exprime la variance de la position de l'espèce dans l'espace de niche, qu'on pourrait envisager comme une mesure de spécialisation
# Rtol : la tolérance résiduelle, l'inertie de l'espèce (= la proportion de sa position dans l'espace de niche) qui n'est ni liée à sa marginalité, ni à sa tolérance

#-------------------------------------------------------------------------------------------#
#### Cas d'étude n°5 : caractéristiques acoustiques et morphologiques des oiseaux de mer ####
#-------------------------------------------------------------------------------------------#

# Cet exemple est traité en anglais pour les participants à la formation non francophones d'origine. 

# data

petrels0 <- read.csv2("petrel_acoustique.csv",dec=".")

## You can notice that there are multiple recordings per individuals. Now let's keep simple. Below is just a small trick to keep one line per biological individual. This operation is done in order to ease the
# analysis during the workshop. In the real life, you should NEVER ignore intra-individual variability. There are ways to account
# for it, we might discuss it later if we have time or if you want to learn more about it. Some parts of the script below deal for this specific point.

# first option : chose randomly one recording per biological individual

petrels <- do.call(rbind, lapply(split(petrels0, petrels0$Indiv), function(x){x[sample(seq_along(x$Indiv), 1),]}))

# second option (which we will use here) : aggregate records per biological individual (here we take the median values across recordings)

petrels <- aggregate(petrels0[,c("Ph.NbSy","Ph.Rythme","Ph.Du","Ph.F0","Ph.Q50","NbPh")],by=list(petrels0$Indiv,petrels0$Species),FUN="median")
colnames(petrels)[1:2] <- c("Indiv","Species")

# with ADE4

petrel1 <- petrels[,-c(1:2,9:10)]

pc.petrel <- dudi.pca(petrel1) # in the R console, you select the number of axes you wish, then type "Enter"
pc.petrel <- dudi.pca(petrel1,scannf=F,nf=3) # same thing, if you already know how many axes you want to keep

# plotting the individuals in the two first dimensions

s.label(pc.petrel$li,xax=1,yax=2)

# in dimensions 2 vs 3

s.label(pc.petrel$li,xax=1,yax=3)

# and dimensions 1 vs 3

s.label(pc.petrel$li,xax=2,yax=3)

# correlation circle (shows how variables are correlated together and are correlated to the axes)

s.corcircle(pc.petrel$co)

fviz_pca_var(pc.petrel) # same plot with another design, maybe a bit easier to read
fviz_pca_var(pc.petrel,axes=c(1,3))
fviz_pca_var(pc.petrel,axes=c(2,3))

## so now, how are the species separated? 

petrels$Species <- factor(petrels$Species)
s.class(pc.petrel$li,petrels$Species,col=c("darkblue","darkred"))
s.class(pc.petrel$li,petrels$Species,col=c("darkblue","darkred"),xax=1,yax=3)

# you may wish to test this formally with a discriminant analysis

discr.petrel <- discrimin(pc.petrel,petrels$Species,scannf=F,nf=2) 

# permutation testing

test <- randtest(discr.petrel)
plot(test)

# projections

plot(discr.petrel)

## morphometric data

morpho <- read.csv2("petrel_biometrie.csv")
rownames(morpho) <- morpho$Indiv

# just in case : check that the individuals are the same

setdiff(morpho$Indiv,petrels$Indiv)

# this new table has to be ordered like the acoustic table 

morpho <- morpho[order(rownames(petrels)),]

## Principal Component Analysis on morphometric measurements

pc.morpho <- dudi.pca(morpho[,2:7],scannf=F,nf=2)

s.corcircle(pc.morpho$co)
fviz_pca_var(pc.morpho) 
s.label(pc.morpho$li)

# We may look at how species are separated in the morphospace (!! individuals need 
# to be in the same order as in 'petrels')

s.class(pc.morpho$li,petrels$Species)

## Coinertia analysis

coi.petrels <- coinertia(pc.morpho,pc.petrel,scannf=F,nf=2) # we may keep just one axis here, but it is more practical to keep 2 for graphical display

# permutation testing --> are covariance structures stronger than expected by chance?

test.coi <- randtest(coi.petrels)
plot(test.coi) 

# graphical summary

plot(coi.petrels)

# projection of the acoustic and morphological PCA axes on the new axes of the coinertia

s.corcircle(coi.petrels$aX, sub="acoustic PCA axes")
s.corcircle(coi.petrels$aY, sub="morphometric PCA axes")

# projection of the variables on the coinertia space

s.arrow(coi.petrels$l1)
s.arrow(coi.petrels$c1)
s.class(coi.petrels$lX,petrels$Species)

# Next part : accounting for intra-individual heterogeneity in sounds.

# This part won't probably be developed during the workshop but it is important that you consider it if you are to work on similar data.
# The main point is that in the original data, some individuals are better sampled than others, with huge variation from an individual to another : 

tapply(petrels0$Indiv,INDEX=petrels0$Indiv,FUN="length")

# In the previous analysis, we ignored this. However, intra-individual variation in sound characteristics could well exceed interspecific variance. 
# If we ignore this possibility, we may over-assess the difference in the songs of the two species (=interspecific variance in acoustic variables)
# Although in this precise example the species are really well separated, there will be many instances where there is overlap in the acoustic parameters
# of two species (or among several populations, or two sexes, or any other classification). In these situations, you want to account for intra-individual 
# variance in order to avoid type 2 error (incorrectly assessing differences among categories that don't actually exist). 

# For this, you can first examine between-individual variation in acoustic parameters, accounting for the fact that there are several recordings per individual.
# This comes to find a representation of the songs that focuses on inter-individual variance. This is done by a modification of the principal component analysis 
# called the "between class analysis" : 

pc.petrel.allind <- dudi.pca(petrels0[,4:ncol(petrels0)],scannf=F,nf=5) # first we do a "standard" PCA on all the data (ie not the aggregated data used so far : here we have all recordings for all individuals). We've set nf=5 but it does not matter
bca.petrel <- bca(pc.petrel.allind,fac=factor(petrels0$Indiv),scannf=F,nf=2) # the between class analysis specifies a grouping factor (birds' ring codes)
pc.petrel # this is the PCA on the aggregated data (averages per individuals) used in the workshop

# compare the ordinations for the pca and the bca in order to grasp what was done : 

s.label(pc.petrel.allind$li)
x11()
s.label(bca.petrel$li) # the ordination plot only represents the birds, not each recording as done in the PCA

s.class(pc.petrel.allind$li,fac=factor(petrels0$Indiv)) # here we've plotted the birds on the PCA. Compare it to the BCA : the positions of birds in the BCA correspond to the centroids of the ellipses.

s.corcircle(pc.petrel.allind$co)
x11()
s.corcircle(bca.petrel$co) # in this example, differences between the correlation circles of the PCA with all data and the BCA are slight. They will be higher in instances when intra-individual variance is higher
x11()
s.corcircle(pc.petrel$co) # on the aggregated PCA, ignoring inter-individual variance, the role of Ph.Q50 (50% energy quantile) on the first axis was over-assessed.
s.label(pc.petrel$li)  # and as expected, the separation between the two species in the ordination space was more marked. 

# Now let's do our discriminant analysis on the BCA, so that to distinguish the two species once intra-individual variance has been accounted for

discrim.petrel.intra <- discrimin(bca.petrel,fac=factor(petrels$Species),scannf=F,nf=2)
plot(discrim.petrel.intra)

# We can move on with the comparison with morphological data. We may perform the coinertia analysis with the bca : 

coi.petrels.bca <- coinertia(bca.petrel,pc.morpho) 

# This does not work. The most sampled individuals (birds with more recordings) have more weight in the BCA, because their intra-individual variance is better estimated. 
# However, there's only one morphological measurement per bird, so there is an inconsistency in the two matrices we want to compare. If we want it to work, we need to state 
# which individuals will have more weight in the morphological PCA also. There are several ways of doing. Here is one : 

pc.morpho.reweight <- dudi.pca(morpho[,2:7],scannf=F,nf=2,row.w=bca.petrel$lw) # we chose to apply the weights of the BCA to the morphological PCA. Thus, the birds that have more recordings will have more weight in both table


# we could also decide that we do NOT want the individuals sampled more times to weight more, in which case we would do a PCA on the BCA in order to equalize weights among individuals : 
pca.bca <- dudi.pca(bca.petrel$tab,scannf=F,nf=2) # note : this comes to doing the PCA on the individual-averaged data that we were using during the workshop, so the interest of that much complexity is limited if you choose this option

# Now that the weights of individuals are equal on the two tables, we can perform the coinertia : 
coi.petrels.bca <- coinertia(bca.petrel,pc.morpho.reweight,scannf=F,nf=2)

# you won't be able to do the permutation testing with this coinertia, since it is not implemented in ade4 yet. It may come later.
# projection of the variables on the coinertia space
s.arrow(coi.petrels.bca$l1)
s.arrow(coi.petrels.bca$c1)
s.class(coi.petrels.bca$lX,fac=petrels$Species) # this coinertia accounts for the intra-individual variance through the BCA
s.class(coi.petrels$lX,fac=petrels$Species)  # and this coinertia does not (it is the one used in the workshop)

# You can see how the weighing changes the interpretation. When you account for intra-individual variance, you learn that variations
# among individuals in prions is explained by bird size (wing length, surface, mass) and some song properties (rythm, frequency, nb of syllabs)
# while variation among individuals in blue petrels is explained by the 50% energy quantile and number of phrases and bill length

# when we didn't account for intra-individual variance (example of the workshop), the individuals of both species differed along the same axis defined 
# by another combination of variables (see above, in section 2). This should warn you about how important intraindividual (or intraspecific and in general within-group
#) variation is in statistical inference. 


#-----------------------------------------------------------------------------------------#
#### Cas d'étude n°6 : peuplements malacologiques de la vallée de Chaudefour, Auvergne ####
#-----------------------------------------------------------------------------------------#

# Cet exemple ressemble beaucoup au cas d'étude 1 et 4. La grosse différence est qu'on ne 
# connait rien à l'écologie des espèces, il s'agit donc d'une étude réellement exploratoire.

# Le fichier malaco_chaudefour_index.csv contient les noms complets des espèces et des variables environnementales

# si vous terminez suffisamment tôt, utilisez les coordonnées des diverses analyses multivariées pour situer vos compositions
# de communautés dans l'espace par des cartographies. Les coordonnées spatiales des points sont dans le fichier "malaco_chaudefour.csv"

# communautés de mollusques
pplt_foret0 <- read.csv2("malaco_chaudefour.csv",row.names=1)

# on commence par explorer les communautés avec une AFC
afc_moll_foret <- dudi.coa(pplt_foret0,scannf=F,nf=2)

# représentations graphiques
scatter(afc_moll_foret)
s.label(afc_moll_foret$li)
s.arrow(afc_moll_foret$co)

# NB : il peut être utile de supprimer Lehm avant de continuer, car elle écrase complètement l'AFC.
pplt_foret <- pplt_foret0[,-which(colnames(pplt_foret0)=="Lehm")]
afc_moll_foret <- dudi.coa(pplt_foret,scannf=F,nf=2)


# variables environnementales
var_env_foret0 <- read.csv2("environnement_chaudefour.csv",row.names=1)

# note : pour pouvoir continuer avec une ACP, il va falloir enlever la contrainte sur les % de recouvrement qui somment à 100% à chaque ligne
# pour cela, on enlève une des classes de recouvrement

var_env_foret <- var_env_foret0[,-which(colnames(var_env_foret0)=="pourc_musc_m2")]

# ACP sur les variables environnementales
acp_var_env_foret <- dudi.pca(var_env_foret,scannf=F,nf=2)
scatter(acp_var_env_foret)
s.arrow(acp_var_env_foret$co)
s.label(acp_var_env_foret$li)

# Nous allons ensuite tenter une coinertie. Pour cela, il faut pondérer l'ACP par les poids
# de lignes de l'AFC
acp_var_env_foret.rew <- dudi.pca(var_env_foret,scannf=F,nf=2,row.w=afc_moll_foret$lw)
coi.moll=coinertia(afc_moll_foret,acp_var_env_foret.rew,scannf=F,nf=2)

# exploration de la coinertie

scatter(coi.moll)
s.label(coi.moll$co)
s.label(coi.moll$li)

# On peut aussi passer par une analyse OMI (voir commentaires dans l'exemple 4)

OMI <- niche(acp_var_env_foret,pplt_foret,scannf=F,nf=2)
niche.param(OMI)
plot(OMI)
scatter(OMI)
s.label(OMI$li)
s.arrow(OMI$co)

# si vous voulez cartographier les coordonnées des analyses (on prend ici pour exemple l'axe 1 de l'ACP,
# mais vous pouvez aussi le faire avec l'OMI, la coinertie, l'AFC... Un bon exercice serait de reprendre 
# la tolérance de chaque espèce issue de l'OMI, la moyenner en chaque point d'échantillonnage, et ainsi étudier
# la variabilité spatiale de la tolérance/spécialisation des communautés)

xy.malaco <- read.csv2("malaco_chaudefour_coordonnees.csv",row.names=1)
df.malaco <- merge(xy.malaco,acp_var_env_foret$li,by=0)
colnames(df.malaco)[1] = "cellname"
df.malaco$cellname = as.character(df.malaco$cellname)

# carte (passez la souris sur les points pour voir leur code "cellname")
carte_acp1 <- df.malaco%>%
  mapview(zcol = "Axis1",
          xcol="x",ycol="y",
          legend = TRUE,
          color = 'black',
          lwd = 2,
          map.types = 'OpenStreetMap',
          at = c(0,2,4,6,8,10),
          layer.name = 'ACP - axe 1',
          crs=2154,
          label="cellname")

