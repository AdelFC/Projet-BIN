########################################### PROJET : dynamique d'accrochage et de decrochage de parp1


import numpy as np

# Pour les graphiques :
import matplotlib.pyplot as plt
import matplotlib.patches as patches # permet d'ajouter le rectangle noir sur le graphique
import pandas as pd # permet de compiler les donnees en un tableau (necessaire pour faire le graphique avec les types de differentes couleurs)
import seaborn as sns # permet d'attribuer des couleurs selon les types
import scipy.signal # permet de creer et appliquer un filtre sur les donnees pour obtenir une courbe plus propre et visuelle




############################ Fonctions ####################################################################################################################################

### Fonction diffusion :
def diffusion(coef_diff, posX_t, posY_t):
    # np.random.normal(mu, sigma, nombre de tirage) : augmenter sigma si on veut augmenter l'amplitude (avoir plus de valeurs extremes)
    epsilonx = (coef_diff)*(np.random.normal(0, 0.15, 1)) 
    epsilony = (coef_diff)*(np.random.normal(0, 0.15, 1))

    posX_t_plus_1 = posX_t + epsilonx
    posY_t_plus_1 = posY_t + epsilony


## Si la molecule depasse les bords du carre :
#1e version : répeter la commande jusqu'à avoir un mouvement possible (qui ne fait pas sortir du carre)
    while (posX_t_plus_1 < 0 or posX_t_plus_1 > 1) or (posY_t_plus_1 < 0 or posY_t_plus_1 > 1) :
        epsilonx = (coef_diff)*(np.random.normal(0, 0.15, 1))
        epsilony = (coef_diff)*(np.random.normal(0, 0.15, 1))

        posX_t_plus_1 = posX_t + epsilonx
        posY_t_plus_1 = posY_t + epsilony
        
    return((posX_t_plus_1, posY_t_plus_1))

# 2e version : prendre l'oppose d'epsilon (-epsilon au lieu de +epsilon)
#     if (posX_t_plus_1 < 0 or posX_t_plus_1 > 1):
#         posX_t_plus_1 = posX_t - epsilonx
#     if (posY_t_plus_1 < 0 or posY_t_plus_1 > 1):
#         posY_t_plus_1 = posY_t - epsilony
#     # mais il y a un risque pour que la molecule sorte quand meme avec cette option
#     # car le deplacement n'est pas verifie dans ce cas
#     # return((posX_t_plus_1, posY_t_plus_1))
  
# 3e version : la molecule reste a sa position (immobile)
#     if (posX_t_plus_1 < 0 or posX_t_plus_1 > 1):
#         posX_t_plus_1 = posX_t
#     if (posY_t_plus_1 < 0 or posY_t_plus_1 > 1):
#         posY_t_plus_1 = posY_t
#     return((posX_t_plus_1, posY_t_plus_1))



### Fonction courbe_gaussienne :
# affiche la courbe selon les parametres choisis
def courbe_gaussienne(mu, sigma, taille_tirage) :
    gauss=np.random.normal(mu, sigma, taille_tirage)
    count, bins, ignored=plt.hist(gauss, 30, density=True)
    plt.plot(bins, 1/(sigma * np.sqrt(2*np.pi))* np.exp(-(bins-mu) **2 / (2*sigma**2)), linewidth=2, color='r')
    plt.show


### Exemple d'utilisation avec les parametres suivants :
# moyenne = 0
# ecart_type = 0.15
# nb = 1000
# courbe_gaussienne(moyenne, ecart_type, nb)





### Fonction immobile :
def immobile(posX_t, posY_t):
    posX_t_plus_1 = posX_t 
    posY_t_plus_1 = posY_t
    return((posX_t_plus_1, posY_t_plus_1))





############################ Parametres ################################################################################################################################
D = 1.633  # facteur qui multiplie epsilon pour obtenir le pas de deplacement de la molecule
# + D est grand : + le pas est grand

nmol = 3000  # nombre de molecules
pt = 840  # nombre de points de temps
# Dans notre cas, les molecules se deplacent sur 420 secondes (7 minutes) 
# Pour etre plus precis, on prend le double de point de temps
# Donc il y aura 1 point toutes les 0.5 secondes

tpre = 10 # temps ou les detags sont induit par FRAP 
# dans la realite : FRAP 1 image avant donc 5 secondes avant le debut de l'acquisition des donnes
# 5 secondes = 10 points de temps

# Limites de la zone
minX = 0.4
maxX = 0.5
minY = 0.4
maxY = 0.5


######### Probabilite des differentes reactions de PARP1 

#### Type 0 = PARP1 : liaison ou pas liaison ?
# Reaction : PARP1 + c --> PARP1-c
# PARP1-c = type 1
Pbound = 0.35 # probabilite de la liaison
# + Pbound petit : + chance de passer en type 1 (liaison favorisée)
# Dans notre cas, cette probabilite est supposee importante donc nous choisissons un Pbound plutot petit



#### Type 1 = PARP1-c : autoparylation, decrochage ou reste tel quel ?

# Decrochage :  PARP1-c --> PARP1 + c (revient au type 0)
Punbound = 0.2 # probabilite du decrochage
# + Punboud petit : - chance de se decrocher et repasser en type 0
# On considere que cette reaction est peu probable et que c'est la parylation qui est favorisee
# Donc on choisit un Punboud assez petit

# Autoparylation : PARP1-c --> PARP1*-c
# PARP1*-c = type 2
Pautop = 0.8 # probabilite d'autoparylation
# + Pautop grand : + chance de s'autoparyler
# On considere que c'est ce que la molecule a le plus de chance de faire une fois accrochee
# Donc on choisit un Pautop assez grand

# Attention, condition a respecter pour la suite : Punbound + Pautop = ou < 1 

#### Si on veut simuler la presence de'inhibiteur : 
# PARP1 reste fixe : type 1 reste type 1 donc ne se paryle pas et ne se decroche pas
# Donc mettre Pautop = 0 et Punbound = 0



#### Type 2 : PARP1*-c : se deparyle, se decroche ou reste tel quel ?

# Deparylation : PARP1*-c --> PARP1-c (revient au type 1)
Pdep = 0.001 #  probabilite de deparylation
# + Pdep petit : - chance de se deparyler en restant accroche
# On considere cette possibilite quasiment nulle, très improbable, donc on choisit un Pdep tres proche de 0


# Decrochage en restant paryle : PARP1*-c --> PARP1* + c
# PARP1*c = type 3 
Punbound_p = 0.9 # probabilite de decrochage en restant paryle
# + Punbound_p grand : + chance de se decrocher en restant paryle
# On considere que c'est la reaction majoritaire qui est tres favorisee, donc on choisit un Punboud_p tres grand

# Attention, condition a respecter pour la suite : Pdep + Punbound_p = ou < 1






############################ Initialisation #########################################################################################################################
# Matrices remplies de 0 (pt lignes et nmol colonnes) et attribution des positions initiales (aleatoires)
posX = np.zeros((pt, nmol))
posX[0] = np.random.random(nmol)
posY = np.zeros((pt, nmol))
posY[0] = np.random.random(nmol)


# Vecteurs :
# Compter le nombre de molecules de chaque type (dans tout le noyau):
compte_0 = np.zeros(pt)
compte_1 = np.zeros(pt)
compte_2 = np.zeros(pt)
compte_3 = np.zeros(pt)

# Compter le nombre de molecules (toutes ou selon les types) rentrees dans la zone
comptezone = np.zeros(pt) # toutes les molecules
comptezone_type0 = np.zeros(pt)
comptezone_type1 = np.zeros(pt)
comptezone_type2 = np.zeros(pt)
comptezone_type3 = np.zeros(pt)

# Attribution du type 0 (PARP1 libre et non paryle) a toutes les molecules (par defaut)
types=np.zeros(nmol)

# Sauvegarde les types en fonction du temps
state=np.zeros((pt,nmol))




################################## Boucle #############################################################################################################################
# noyau : carre de 0 a 1
# zone de degats : carre avec x et y entre 0.4 et 0.5


# Avant le temps tpre :
# Pas encore de degats donc toutes les molecules sont en type 0, diffusent librement (sans se lier donc sans se paryler)
for i in range(0, tpre):
    cpt_zone = 0
    compte_0 [i] = nmol
    
    for j in range (0,nmol) :
        if (minX < posX[i][j] < maxX and minY < posY[i][j] < maxY):
            cpt_zone = cpt_zone + 1   # compte +1 pour chaque molecule dans la zone
        posX[i+1][j], posY[i+1][j] = diffusion(D, posX[i][j], posY[i][j])
    comptezone[i] = cpt_zone


# A partir du temps tpre (degats) :
# La liaison a la cassure est rendue possible, donc la parylation qui s'en suit egalement
for i in range(tpre, pt-1) :
    cpt_type0 = 0 # compter nombre de molecules de chaque type dans tout le noyau
    cpt_type1 = 0
    cpt_type2 = 0
    cpt_type3 = 0    
    cpt_zone = 0 # compter nombre total dans la zone
    cpt_zone0 = 0 # compter nombre de molecules de chaque type, uniquement dans la zone de degat
    cpt_zone1 = 0
    cpt_zone2 = 0
    cpt_zone3 = 0
            
    
    for j in range (0, nmol):
        if (minX < posX[i][j] < maxX and minY < posY[i][j] < maxY) : # si molecule entre dans la zone
            cpt_zone = cpt_zone + 1 
            nb_aleatoire = np.random.random(1) # tire un nombre aleatoire a utiliser pour la suite 
               
            if (types[j] == 0) : # PARP1 libre
                cpt_type0 = cpt_type0 +1
                cpt_zone0 = cpt_zone0+1
                
                if (nb_aleatoire < Pbound) : # pas de liaison : reste type 0 et continue de diffuser
                    types[j] = 0
                    posX[i+1][j], posY[i+1][j] = diffusion(D, posX[i][j], posY[i][j])
            
                else : # si  nombre aleatoire > Pbound, PARP1 se lie et devient de type 1 (PARP-c : immobile)
                    types[j] = 1
                    posX[i+1][j], posY[i+1][j] = immobile(posX[i][j], posY[i][j])
            
            
            elif (types[j] == 1) : # PARP1-c
                cpt_type1 = cpt_type1 + 1
                cpt_zone1 = cpt_zone1+1
                
                if (nb_aleatoire < Punbound): # se dissocie et redevient de type 0 : PARP1 libre
                    types[j] = 0
                    posX[i+1][j], posY[i+1][j] = diffusion(D, posX[i][j], posY[i][j])
                
                elif (Punbound < nb_aleatoire < (Punbound + Pautop)):  # PARP1 s'autoparyle et devient type 2 : PARP*-c 
                    types[j] = 2
                    posX[i+1][j], posY[i+1][j] = immobile(posX[i][j], posY[i][j])
            
                else : # nombre aleatoire > (Punbound + Pautop) : PARP1 reste associee en type 1 : PARP1-c
                    types[j] = 1
                    posX[i+1][j], posY[i+1][j] = immobile(posX[i][j], posY[i][j])
        
        
            elif (types[j] == 2) : # PARP1*-c
                cpt_type2 = cpt_type2 + 1
                cpt_zone2 = cpt_zone2 +1
                
                if (nb_aleatoire < Pdep) : # PARP1 se deparyle et revient en type 1 : PARP1-c
                    types[j] = 1
                    posX[i+1][j], posY[i+1][j] = immobile(posX[i][j], posY[i][j])
                
                elif (Pdep < nb_aleatoire < (Pdep + Punbound_p)): # PARP1 se decroche et devient type 3 : PARP1*
                    types[j] = 3
                    posX[i+1][j], posY[i+1][j] = diffusion(D, posX[i][j], posY[i][j])
                
                else : # si nb > (Pdep + Punbound_p) : PARP1 reste type 2 : PARP1*-c
                    types[j] = 2
                    posX[i+1][j], posY[i+1][j] = immobile(posX[i][j], posY[i][j])
        
        
            else : # type = 3 : PARP1* (diffuse sans suivre d'autre reaction d'apres notre modele)
                cpt_type3 = cpt_type3 + 1
                cpt_zone3 = cpt_zone3+1
                
                posX[i+1][j], posY[i+1][j] = diffusion(D, posX[i][j], posY[i][j])
        
        
        else: # molecules pas dans la zone : restent en type 0 ou en type 3 et diffusent
            posX[i+1][j], posY[i+1][j] = diffusion(D, posX[i][j], posY[i][j])    
            if (types[i] == 0):
                cpt_type0 = cpt_type0 +1
            else :
                ct_type3 = cpt_type3 +1
           
            
# attribution des valeurs aux conteurs pour chaque temps i
    state[i] = types
    compte_0 [i] = cpt_type0 # nombre de molecules de chaque type dans le noyau
    compte_1 [i] = cpt_type1
    compte_2 [i] = cpt_type2
    compte_3 [i] = cpt_type3    
    comptezone[i] = cpt_zone  # nombre de molecules totales dans la zone
    comptezone_type0 [i] = cpt_zone0 # nombre de molecules de chaque type dans la zone
    comptezone_type1 [i] = cpt_zone1
    comptezone_type2 [i] = cpt_zone2
    comptezone_type3 [i] = cpt_zone3
    





############################ Visualisation ##########################################################################################################################

####### Nombre de molecule de chaque types :
print("le nombre de molécule dans la zone à la fin de la simulation est : ")
print(comptezone[pt-2]) # nombre de molecules dans la zone à la fin

print("le nombre de molécules libres et non parylées (type 0) à la fin de la simulation est de :")
print(cpt_type0)

print("le nombre de molécules liées et non parylées (type 1) à la fin de la simulation est de :")
print(cpt_type1)

print("le nombre de molécules auto-parylées et liées (type 2) à la fin de la simulation est de :")
print(cpt_type2)

print("le nombre de molécules auto-parylées et libre (type 3) à la fin de la simulation est de :")
print(cpt_type3)




###### Trajectoire d'une molecule au cours du temps :
# creer des listes vides pour X et Y (autre maniere que de remplir par des zeros comme on faisait avant) :
trajX = [] 
trajY = []

# Pour la 1e molecule (indice = 0)
# Pour chaque point de temps : ajouter les positions X et Y dans la liste "traj"
for i in range (0, pt-1) : 
    trajX.append(posX[i][0]) # on ajoute dans la liste la position (a la fois x et y) de l'element 0 (premiere molecule)
    trajY.append(posY[i][0])

# Afficher la trajectoire en reliant les points :
plt.plot(trajX, trajY) 
plt.xlabel("X")
plt.ylabel("Y")
plt.title("Trajectoire d'une molecule au cours du temps")
plt.show()
# On peut aussi choisir de n'afficher que les points (sans les relier) avec "plt.scatter" au lieu de "plt.plot"





###### Positions des molecules a plusieurs moments de la simulation :
# Choisir le temps a regarder (nombre entier) :
    
### avant les degats :
pt_avant = int(pt-(pt-2))

# Compiler les donnees a regarder dans un seul tableau pour pouvoir lier les positions au type de chaque molecule
compil_pt_avant = {"X" : posX[pt_avant] , "Y" : posY[pt_avant], "State" : state[pt_avant]}
compil_pt_avant = pd.DataFrame(compil_pt_avant)
# selon le [pt] qu'on met : on choisi quand on regarde

fenetre, figure = plt.subplots()  # on definit la fenetre dans laquelle on va creer la figure (le graphique), pour pouvoir y ajouter la zone apres

# on definit les choses a regarder dans le tableau pour pouvoir separer les types avec les couleurs :
sns.scatterplot(x = "X", 
                y = "Y",
                s = 25, # taille des points
                hue = "State",  # donnees sur lequelles on base la separation des couleurs
                data = compil_pt_avant,  # nom du tableau ou il y a toutes les donnees mentionnees
                palette = "rainbow_r", # palette de couleurs utilisees pour la legende (selon les types)
                alpha = 0.6) # transparence pour plus de lisibilite

# titres des axes et du graph :
plt.xlabel("X")
plt.ylabel("Y")
plt.title("Position des molecules selon leur etat avant les degats")

# Mettre la legende en dehors du graph (sinon elle cache les points) 
plt.legend(bbox_to_anchor = (1.01, 1),
           borderaxespad = 0)

# Rectangle noir : ajouter patch
figure.add_patch(patches.Rectangle((minX, minY), # coordonnees du point en bas a gauche
                                   0.1, # largeur en x
                                   0.1, # hauteur en y
                                   edgecolor = 'black',
                                   facecolor = 'red', 
                                   fill = False))  # false car on remplie pas, on pourrait mettre True avec la couleur de la ligne precedente

## Possibilite d'enregistrer la figure dans en png :
# plt.savefig("NOM.png", format='png', dpi=150)



   
### Au moment des degats :
pt_degat = int(tpre)
    
compil_pt_degat = {"X" : posX[pt_degat] , "Y" : posY[pt_degat], "State" : state[pt_degat]}
compil_pt_degat = pd.DataFrame(compil_pt_degat)

fenetre, figure = plt.subplots()  

sns.scatterplot(x = "X", 
                y = "Y",
                s = 25,
                hue = "State",  
                data = compil_pt_degat,  
                palette = "rainbow_r",
                alpha = 0.6) 

plt.xlabel("X")
plt.ylabel("Y")
plt.title("Position des molecules selon leur etat au moment des degats")

plt.legend(bbox_to_anchor = (1.01, 1),
           borderaxespad = 0)

figure.add_patch(patches.Rectangle((minX, minY),
                                   0.1, 
                                   0.1, 
                                   edgecolor = 'black',
                                   facecolor = 'red', 
                                   fill = False))  

    
    
### A la moitie du temps :
pt_moitie = int(pt/2)

compil_pt_moitie = {"X" : posX[pt_moitie] , "Y" : posY[pt_moitie], "State" : state[pt_moitie]}
compil_pt_moitie = pd.DataFrame(compil_pt_moitie)

fenetre, figure = plt.subplots()  

sns.scatterplot(x = "X", 
                y = "Y",
                s = 25,
                hue = "State",  
                data = compil_pt_moitie,  
                palette = "rainbow_r",
                alpha = 0.6) 

plt.xlabel("X")
plt.ylabel("Y")
plt.title("Position des molecules selon leur etat a la moitie du temps")

plt.legend(bbox_to_anchor = (1.01, 1),
           borderaxespad = 0)

figure.add_patch(patches.Rectangle((minX, minY), 
                                   0.1, 
                                   0.1,
                                   edgecolor = 'black',
                                   facecolor = 'red', 
                                   fill = False ))



### A la fin de la simulation :
pt_fin = int(pt-2)
    
compil_pt_fin = {"X" : posX[pt_fin] , "Y" : posY[pt_fin], "State" : state[pt_fin]}
compil_pt_fin = pd.DataFrame(compil_pt_fin)

fenetre, figure = plt.subplots()  

sns.scatterplot(x = "X", 
                y = "Y",
                s = 25, 
                hue = "State", 
                data = compil_pt_fin, 
                palette = "rainbow_r",
                alpha = 0.6) 

plt.xlabel("X")
plt.ylabel("Y")
plt.title("Position des molecules selon leur etat a la fin du temps")

plt.legend(bbox_to_anchor = (1.01, 1),
           borderaxespad = 0)

figure.add_patch(patches.Rectangle((minX, minY), 
                                   0.1, 
                                   0.1, 
                                   edgecolor = 'black',
                                   facecolor = 'red', 
                                   fill = False ))  

plt.show()





####### Nombre de molecule dans la zone selon leur type au cours du temps :

### Tous les types dans la zone en une meme couleur (noir) :
plt.scatter(list(np.linspace(start=0, stop=(pt/2)-1, num=pt-1)), comptezone[0: len(comptezone)-1], c='#333333') 
# -1 pour les deux variables sinon affiche un temps de trop avec 0 molecules dedans (car dans la boucle on fait des i+1)
# list(np.linspace) : creer la liste "temps" qui va de 0 a pt-1 points avec un pas de 0.5 (on va de 0 a 419 avec 839 points)

# Titres et legendes :
plt.xlabel("Secondes")
plt.ylabel("Nombre de molecules dans la zone")
plt.title("Nombre de molecules (tout types) dans la zone au cours du temps")
plt.show()


# Courbe plus lisible en filtrant le signal :
# filtre passe-bas : enleve hautes frequentes (pour supprimer le bruit)
h = scipy.signal.firwin(numtaps=2*100+1,cutoff=1/30,nyq=0.5,window='hann')
# paramètres du filtre : mettre en forme grace a cutoff = frequence de coupure (coupe les frequences au dessus de ce seuil)

# tracer le signal filtre : passe le compte zone dans le filtre --> y2
y2 = scipy.signal.convolve(comptezone,h,mode='same')

# visualisation du comptezone filtre
plt.plot(list(np.linspace(0, (pt/2)-1, pt)),y2, color="Black")
plt.xlabel("Secondes")
plt.ylabel("Nombre de molecules dans la zone")
plt.title("Nombre de molecules dans la zone au cours du temps")
plt.show()



### Afficher separemment chaque type dans la zone :
# type 0 points :
plt.scatter(list(np.linspace(0, (pt/2)-1, pt-1)), comptezone_type0[0: len(comptezone_type0)-1]) 
plt.xlabel("Secondes")
plt.ylabel("Nombre de types 0 dans la zone")
plt.title("Nombre de molecules de type 0 dans la zone au cours du temps")
plt.show()

# type 0 courbe (reutilise le meme h que precedemment car meme filtre) :
y2_type0 = scipy.signal.convolve(comptezone_type0,h,mode='same')
plt.plot(list(np.linspace(0, (pt/2)-1, pt)),y2_type0)
plt.xlabel("Secondes")
plt.ylabel("Nombre de types 0 dans la zone")
plt.title("Nombre de molecules de type 0 dans la zone au cours du temps")
plt.show()



# type 1 points :
plt.scatter(list(np.linspace(0, (pt/2)-1, pt-1)), comptezone_type1[0: len(comptezone_type1)-1], c='darkorange')
plt.xlabel("Secondes")
plt.ylabel("Nombre de types 1 dans la zone")
plt.title("Nombre de molecules de type 1 dans la zone au cours du temps")
plt.show()

# type 1 courbe : 
y2_type1 = scipy.signal.convolve(comptezone_type1,h,mode='same')
plt.plot(list(np.linspace(0, (pt/2)-1, pt)),y2_type1, color="darkorange")
plt.xlabel("Secondes")
plt.ylabel("Nombre de types 1 dans la zone")
plt.title("Nombre de molecules de type 1 dans la zone au cours du temps")
plt.show()


# type 2 points :
plt.scatter(list(np.linspace(0, (pt/2)-1, pt-1)), comptezone_type2[0: len(comptezone_type2)-1], c='Green') 
plt.xlabel("Secondes")
plt.ylabel("Nombre de types 2 dans la zone")
plt.title("Nombre de molecules de type 2 dans la zone au cours du temps")
plt.show()

# # type 2 courbe :
y2_type2 = scipy.signal.convolve(comptezone_type2,h,mode='same')
plt.plot(list(np.linspace(0, (pt/2)-1, pt)),y2_type2, color="Green")
plt.xlabel("Secondes")
plt.ylabel("Nombre de types 2 dans la zone")
plt.title("Nombre de molecules de type 2 dans la zone au cours du temps")
plt.show()    


# type 3 points :
plt.scatter(list(np.linspace(0, (pt/2)-1, pt-1)), comptezone_type3[0: len(comptezone_type3)-1], c='Brown') 
plt.xlabel("Secondes")
plt.ylabel("Nombre de types 3 dans la zone")
plt.title("Nombre de molecules de type 3 dans la zone au cours du temps")
plt.show()

# # type 3 courbe :
y2_type3 = scipy.signal.convolve(comptezone_type3,h,mode='same')
plt.plot(list(np.linspace(0, (pt/2)-1, pt)),y2_type3, color="Brown")
plt.xlabel("Secondes")
plt.ylabel("Nombre de types 3 dans la zone")
plt.title("Nombre de molecules de type 3 dans la zone au cours du temps")
plt.show()




### Tous les types dans la zone, sur le meme graphe avec differentes couleurs selon le type :
plt.scatter(list(np.linspace(0, (pt/2)-1, pt-1)), comptezone_type0[0: len(comptezone_type0)-1], s=8, alpha=0.7) 
plt.scatter(list(np.linspace(0, (pt/2)-1, pt-1)), comptezone_type1[0: len(comptezone_type1)-1], s=8, alpha=0.7)
plt.scatter(list(np.linspace(0, (pt/2)-1, pt-1)), comptezone_type2[0: len(comptezone_type2)-1], s=8, alpha=0.7) 
plt.scatter(list(np.linspace(0, (pt/2)-1, pt-1)), comptezone_type3[0: len(comptezone_type3)-1], s=8, alpha=0.7) 
# s = taille des points, alpha = transparance des points (1=opaque)
plt.legend(("Type0", "Type1", "Type2", "Type3"))
plt.xlabel("Secondes")
plt.ylabel("Nombre de molecules dans la zone")
plt.title("Nombre de molecules dans la zone au cours du temps")
plt.show()


## Courbe : 
y2_type0 = scipy.signal.convolve(comptezone_type0,h,mode='same')
plt.plot(list(np.linspace(0, (pt/2)-1, pt)),y2_type0)

y2_type1 = scipy.signal.convolve(comptezone_type1,h,mode='same')
plt.plot(list(np.linspace(0, (pt/2)-1, pt)),y2_type1)

y2_type2 = scipy.signal.convolve(comptezone_type2,h,mode='same')
plt.plot(list(np.linspace(0, (pt/2)-1, pt)),y2_type2)

y2_type3 = scipy.signal.convolve(comptezone_type3,h,mode='same')
plt.plot(list(np.linspace(0, (pt/2)-1, pt)),y2_type3)

plt.legend(labels=["Type0", "Type1", "Type2", "Type3"])
plt.xlabel("Secondes")
plt.ylabel("Nombre de types 3 dans la zone")
plt.title("Nombre de molecules de type 3 dans la zone au cours du temps")
plt.show()

