% Pour importé un fichier : uiimport
% a regader : unique
% uiimport

%% Introduction

A = EntrezGeneInteractorA;
B = EntrezGeneInteractorB;

A1 = unique(A); %creation d'un vecteur d'unique pour les genes A
B1 = unique(B); %creation d'un vecteur d'unique pour les genes B

ones1 = unique([A1;B1]);  %creation d'un vecteur d'unique pour les 2 gene

M_A= zeros(length(ones1)-2); %création d'une matrice M_A remplie de 0 de coté = a la longueur de ones -2 car les deux derniers termes de ones sont NA 


ones1([length(ones1),length(ones1)-1])=[]; %Supprimer les deux derniere valeurs (qui sont NA et NA)





%% construction de la matrice d'adjacence A

for i = 1 : 1 : length(ones1)
    z1 = find(A == ones1(i)); %chercher les indices dans A des valeurs uniques de ones
    z2 = find(B == ones1(i)); %chercher les indices dans B des valeurs uniques de ones
    
    if ~isempty(z1) && ~isempty(z2)  %Si il y au moins une valeur qu'on cherche dans A et dans B
        for t = 1 : 1 : length(z1) 
            f =find(ones1 == B( z1(t) )); %chercher la valeur de associée a B(z1(t)) dans le vecteur ones pour la coordonée j dans la matrice final
            M_A( i, f(1) )= 1 + M_A( i, f(1) ); %pour mettre les valeurs dans la matrice
            M_A( f(1),i ) = 1 + M_A( f(1),i ); %pour mettre les valeurs dans la partie haute de la matrice pour quelle soit symetrique
        end
    end
    
end;







%% le degré de centralité
%==========================================================================================================

%% degré de centralité
n = sum(M_A);
Cd = zeros(length(ones1),1);
max_interaction = sum(sum(M_A))-1;

for i = 1 : 1 : length(ones1)
    Cd(i) = n(i) /max_interaction;
end; % trouvé toutes les centralités de degré

max_int=ones1(Cd==max(Cd)) %pour trouver la proteine correspodante

%% centralité par valeur propre
[vecM_A,valM_A]=eigs(M_A);  %vecM_A est la matrice de passage avec les vecteurs propres associé aux valeurs propres dans l'ordre croissant,
                            %valM_A est la matrice diagonale
Ce = vecM_A(:,1);  %vecteur propre associé a la plus grande valeur propre

max_int = ones1(Ce==max(Ce))

%% regresion lineaire avec les max
p1 = polyfit(Ce,Cd,1);  %retour matrice avec Cd = p1(1)*Ce+p1(2)
plot(Ce,Cd,'o',Ce,p1(1)*Ce+p1(2));
text(0.04,0.06,'0.0827957087531107*Ce+0.000780913497058389')

r1 = abs(Cd-p1(1)*Ce-p1(2));  % residue avec les max 

max_r1=max(r1)

m=mean(r1)

svdCe=ones(length(Ce),2);
svdCe(:,1)=Ce;
[U,S,V]=svd(svdCe);

ab=V*pinv(S)*U'*Cd

r10 = abs(Cd-ab(1)*Ce-ab(2));

max_r10=max(r10)

%% regresion lineaire sans les max

Ce1=Ce;
Cd1=Cd;
Ce1(Ce == max(Ce))=[];
Cd1(Cd == max(Cd))=[];

p2 = polyfit(Ce1,Cd1,1);  %retour matrice avec Cd = p(1)*Ce+p(2)
plot(Ce1,Cd1,'o',Ce1,p2(1)*Ce1+p2(2));
text(-0.15,0.02,'-0.0118808609457656*Ce+0.000223211576972437')

r2=abs(Cd1-(p2(1)*Ce1-p2(2))); %residue sans les max
m=mean(r2)

r2(length(r2)+1)=0;

max_r2=max(r2)

diff_r = abs(r2-r1);

max_diff_r=max(diff_r)

svdCe=ones(length(Ce1),2);
svdCe(:,1)=Ce1;
[U,S,V]=svd(svdCe);

ab=V*pinv(S)*U'*Cd

r20 = abs(Cd1-ab(1)*Ce1-ab(2));

max_r20=max(r20)


%% ranking
%===========================================================================================================================
%%

[~,rCd] = sort(Cd,'descend'); %rCd est le rang des valeurs de Cd sont triée

 
[~,rCe] = sort(Ce,'descend');


%% 


proteine=zeros(40,2);
%Pour voir les proteine correspondant
for i = 1 : 1 : 40
        proteine(i,1) = ones1(rCd(i));
end;

for i = 1 : 1 : 40
        proteine(i,2) = ones1(rCe(i));
end;


cent=zeros(40,1);
for i= 1:1:40
    z1=find(proteine(:,1)==proteine(i,2));
    
    if ~isempty(z1)
        cent(i)=proteine(z1);
    end
end;

cent

































































