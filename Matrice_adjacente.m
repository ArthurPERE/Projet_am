% Pour importé un fichier : uiimport
% a regader : unique
% uiimport

%% Introduction

A = EntrezGeneInteractorA;
B = EntrezGeneInteractorB;

A1 = unique(A); %creation d'un vecteur d'unique pour les genes A
B1 = unique(B); %creation d'un vecteur d'unique pour les genes B

ones = unique([A1;B1]);  %creation d'un vecteur d'unique pour les 2 gene

M_A= zeros(length(ones)-2); %création d'une matrice M_A remplie de 0 de coté = a la longueur de ones -2 car les deux derniers termes de ones sont NA 


ones([length(ones),length(ones)-1])=[];






%% construction de la matrice d'adjacence A%%

for i = 1 : 1 : length(ones)
    z1 = find(A == ones(i)); %chercher les indices dans A des valeurs uniques de ones
    z2 = find(B == ones(i)); %chercher les indices dans B des valeurs uniques de ones
    
    if ~isempty(z1) && ~isempty(z2)  %Si il y au moins une valeur qu'on cherche dans A et dans B
        for t = 1 : 1 : length(z1) 
            f =find(ones == B( z1(t) )); %chercher la valeur de associée a B(z1(t)) dans le vecteur ones pour la coordonée j dans la matrice final
            M_A( i, f(1) ) = 1 + M_A( i, f(1) ); %pour mettre les valeur dans la matrice
            M_A( f(1),i ) = 1 + M_A( f(1),i ); %pour mettre les valeurs dans la partie haute de la matrice pour quelle soit symetrique
        end
    end
    
end;











%% le degré de centralité
%==========================================================================================================

%% degré de centralité
n = sum(M_A);
Cd = zeros(length(ones),1);
max_interaction = sum(sum(M_A))-1;

for i = 1 : 1 : length(ones)
    Cd(i) = n(i) /max_interaction;
end; % trouvé toutes les centralités de degré

ones(Cd==max(Cd)) %pour trouver la proteine correspodante

%% centralité par valeur propre


%enlever par ce que cela prend pas mal de temps
%[vecM_A,valM_A]=eig(M_A);  %vecM_A est la matrice de passage avec les vecteurs propres associé aux valeurs propres dans l'ordre croissant,
                            %valM_A est la matrice diagonale
Ce = vecM_A(:,length(vecM_A));  %vecteur propre associé a la plus grande valeur propre

ones(Ce==max(Ce))

%% regresion lineaire avec les max
p1 = polyfit(Ce,Cd,1);  %retour matrice avec Cd = p(1)*Ce+p(2)
plot(Ce,Cd,'o',Ce,p1(1)*Ce+p1(2));

r1 = abs(Cd-p1(1)*Ce-p1(2));  % residue avec les max 

max(r1)


%% regresion lineaire sans les max

Ce1=Ce;
Cd1=Cd;
Ce1(Ce == max(Ce))=[];
Cd1(Cd == max(Cd))=[];

p2 = polyfit(Ce1,Cd1,1);  %retour matrice avec Cd = p(1)*Ce+p(2)
plot(Ce1,Cd1,'o',Ce1,p2(1)*Ce1+p2(2));

r2=abs(Cd1-(p2(1)*Ce1-p2(2))); %residue sans les max


max(r2)

abs(max(r2)-max(r1))









%% ranking
%===========================================================================================================================
%%

[rCd,rCd1] = sort(Cd,'descend'); %rCd est le vecteur dont valeurs de Cd sont triée
                                 %rCd1 est le vecteur dont les rangs des valeurs de Cd sont trié
 
[rCe,rCe1] = sort(Ce,'descend');


%% 


compare=zeros(40,1);
%compare les 40 1ere valeurs des vecteurs
for i = 1 : 1 : 40
    compare(i,1) = abs(rCd1(i)-rCe1(i));
end;






proteine=zeros(40,2);
%Pour voir les proteine correspondant
for i = 1 : 1 : 40
        proteine(i,1) = ones(rCd1(i));
end

for i = 1 : 1 : 40
        proteine(i,2) = ones(rCe1(i));
end







































































