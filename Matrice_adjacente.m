% Pour importé un fichier : uiimport
% a regader : unique
% uiimport

A = EntrezGeneInteractorA;
B = EntrezGeneInteractorB;

A1 = unique(A); %creation d'un vecteur d'unique pour les genes A
B1 = unique(B); %creation d'un vecteur d'unique pour les genes B

ones = unique([A1;B1]);  %creation d'un vecteur d'unique pour les 2 genes

M_A= zeros(length(ones)-2); %création d'une matrice M_A remplie de 0 de coté = a la longueur de ones -2 car les deux derniers termes de ones sont NA 

z = find(A == ones(1));

length(z);


for i = 1 : 1 : length(ones)-2
    z1 = find(A == ones(i)); %chercher les indices dans A des valeurs uniques de ones
    z2 = find(B == ones(i)); %chercher les indices dans B des valeurs uniques de ones
    
    if ~isempty(z1) && ~isempty(z2)  %Si il y au moins une valeur qu'on cherche dans A et dans B
        for t = 1 : 1 : length(z1) 
            f =find(ones == B( z1(t) )); %chercher la valeur de associée a B(z1(t)) dans le vecteur ones pour la coordonée j dans la matrice final
            M_A( i, f(1) ) = 1 + M_A( i, f(1) ); %pour mettre les valeur dans la matrice
            M_A( f(1),i ) = 1 + M_A( f(1),i ); %pour mettre les valeurs dans la partie haute de la matrice
        end
    end
    
end;
