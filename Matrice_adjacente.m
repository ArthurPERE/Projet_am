% Pour importé un fichier : uiimport
% a regader : unique
% uiimport

A = EntrezGeneInteractorA;
B = EntrezGeneInteractorB;

A1 = unique(A); %creation d'un vecteur d'unique pour les genes A
B1 = unique(B); %creation d'un vecteur d'unique pour les genes B

ones = unique([A1;B1]);  %creation d'un vecteur d'unique pour les 2 genes

M_A= zeros(length(ones)-2); %création d'une matrice M_A remplie de 0 de coté = a la longueur de ones -2 car les deux derniers termes de ones sont NA 


%===========================================================================
% construction de la matrice d'adjacence A%%
%===========================================================================
for i = 1 : 1 : length(ones)-2
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


%% 


%===========================================================================
% le degré de centralité
%===========================================================================
d = find(M_A>0);
n = length(d);
Cd = zeros(length(ones)-2,1);


for i = 1 : 1 : n/2
    Cd(i) = M_A( d(i) ) /(n-1);
end; % trouvé toutes les centralités de degré

%% 


M_A(d(Cd==max(Cd))); % nombre d'interaction = 9
ones(mod(d(Cd==max(Cd)),length(ones)-2)); %pour trouver la proteine correspodante


%enlever par ce que cela prend pas mal de temps
%[vecM_A,valM_A]=eig(M_A);  %vecM_A est la matrice de passage avec les vecteurs propres associé aux valeurs propres dans l'ordre croissant,
                            %valM_A est la matrice diagonale
Ce = vecM_A(:,length(vecM_A));  %vecteur propre associé a la plus grande valeur propre


p1 = polyfit(Ce,Cd,1);  %retour matrice avec Cd = p(1)*Ce+p(2)
%plot(Ce,Cd,'o',Ce,p1(1)*Ce+p1(2));

r1=abs(Cd-p1(1)*Ce-p1(2));  % residue avec les max 




%% 


%======================================
%enlevé les max
%======================================
Ce1=Ce;
Cd1=Cd;
Ce1(Ce == max(Ce))=0;
Cd1(Cd == max(Cd))=0;

p2 = polyfit(Ce1,Cd1,1);  %retour matrice avec Cd = p(1)*Ce+p(2)
plot(Ce1,Cd1,'o',Ce1,p2(1)*Ce1+p2(2));

r2=abs(Cd1-p2(1)*Ce1-p2(2)); %residue sans les max


abs(r2-r1);  %difference des residues de l'ordre du 10^-7










%% 


%===========================================================================
%question ranking
%===========================================================================
[rCd,rCd1] = sort(Cd,'descend'); %rCd est le vecteur dont valeurs de Cd sont triée
                       %rCd1 est le vecteur dont les rangs des valeurs de Cd sont trié
[~,rank_Cd]= sort(rCd1,'descend');
 
[rCe,rCe1] = sort(Ce,'descend');
[~,rank_Ce]= sort(rCe1,'descend');
%% 

rank_Cd([1:9])=[];

%% 


compare=zeros(40,1);
%compare les 40 1ere valeurs des vecteurs
for i = 1 : 1 : 40
    compare(i,1) = abs(rank_Cd(i)-rank_Ce(i));
end;






proteine=zeros(40,1);
%Pour voir les proteine correspondant
for i = 1 : 1 : 40
        proteine(i)=ones(mod(d(Cd==rCd(rank_Cd(i))),length(ones)-2));
end;









































