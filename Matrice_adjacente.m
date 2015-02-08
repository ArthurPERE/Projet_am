%% Introduction

A = EntrezGeneInteractorA;
B = EntrezGeneInteractorB;


ones1 = unique([A;B]);  %creation d'un vecteur d'unique pour les 2 genes
ones1([length(ones1),length(ones1)-1])=[]; % Supprimer les deux derniere valeurs (qui sont NA et NA)





%% construction de la matrice d'adjacence A

M_A= zeros(length(ones1));

indA = zeros(length(A),1);
indB = zeros(length(B),1);

for i = 1 : 1 : length(A)
    indA(i) = find(ones1 == A(i));
end;

for i = 1 : 1 : length(B)
    indB(i) = find(ones1 == B(i));
end;

for i = 1 : 1 : max(length(A),length(B))
    if indA(i) ~= 0 && indB(i) ~= 0
        M_A(indA(i),indB(i)) = M_A(indA(i),indB(i)) + 1;
        M_A(indB(i),indA(i)) = M_A(indB(i),indA(i)) + 1;
    end;
end;


%% le degré de centralité
%==========================================================================================================

%% degré de centralité
n = sum(M_A);
Cd = zeros(length(ones1),1);
max_interaction = sum(sum(M_A))-1;

for i = 1 : 1 : length(ones1)
      Cd(i) = n(i) /max_interaction;
end; % pour trouver toutes les centralités de degré

max_int=ones1(Cd==max(Cd)) %pour trouver la proteine correspondante

%% centralité par valeur propre
[vecM_A,valM_A]=eigs(M_A);  %vecM_A est la matrice de passage avec les vecteurs propres associé aux valeurs propres dans l'ordre décroissant.
                            %valM_A est la matrice diagonale
Ce = vecM_A(:,1);  %vecteur propre associé a la plus grande valeur propre

max_int = ones1(Ce==max(Ce))

%% regresion lineaire avec les max
p1 = polyfit(Ce,Cd,1);  %retour matrice avec Cd = p1(1)*Ce+p1(2)
plot(Ce,Cd,'o',Ce,p1(1)*Ce+p1(2));
text(0.04,0.06,'0.0833757332095258*Ce-0.000278082785080759')

r1 = abs(Cd-p1(1)*Ce-p1(2));  % residus avec les max

max_r1=max(r1)

m=mean(r1)




% Pour calculer les coefficients de la regresion lineaire sans utilisé la fonction polyfit (ce qui donne la même chose)
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
text(0.02,0.02,'0.0256757352049048*Ce+0.000117646007445322')

r2=abs(Cd1-(p2(1)*Ce1-p2(2))); %residus sans les max
m=mean(r2)

r2(length(r2)+1)=0;

max_r2=max(r2)

diff_r = abs(r2-r1);

max_diff_r=max(diff_r)



% Pour calculer les coefficients de la regresion lineaire sans utilisé la fonction polyfit (ce qui donne la même chose)
svdCe1=ones(length(Ce1),2);
svdCe1(:,1)=Ce1;
[U1,S1,V1]=svd(svdCe1);
ab=V1*pinv(S1)*U1'*Cd1

r20 = abs(Cd1-ab(1)*Ce1-ab(2));

max_r20=max(r20)


%% ranking
%===========================================================================================================================
%% Creation des vecteurs de rang

[~,rCd] = sort(Cd,'descend'); %rCd est le rang des valeurs de Cd triés
[~,rCe] = sort(Ce,'descend');

 
%% Voir pour quel proteines cela corespond

proteine=zeros(40,2);
for i = 1 : 1 : 40
    proteine(i,1) = ones1(rCd(i));
    proteine(i,2) = ones1(rCe(i));
end;


%% Voir les proteines en double

cent=zeros(40,1);
for i= 1:1:40
    z1=find(proteine(:,1)==proteine(i,2));

    if ~isempty(z1)
        cent(i)=proteine(z1);
    end;
end;


cent