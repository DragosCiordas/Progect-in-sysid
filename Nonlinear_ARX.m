clearvars;
close all;
load("iddata-01.mat");
%id-array este pt. cei care nu au instalat identification toolbox
%utilizez variabilele id si val
clear id_array;
clear val_array; %nu le utilizez

%vizualizarea datelor de identificare:
figure(1);
plot(id);

y_id = id.y;
u_id = id.u;

na = 2;
nb = 2;
m = 4;

%%
vector_puteri = combinare_unica(na,nb,m);
DKid = generare_PHI(id,na,nb);  %vectorul de semnale precedente
                                %echivalent cu PHI din ARX (luat de la
                                %lab6)
PHIid = phi_narx(vector_puteri,DKid,length(y_id));
THETA = PHIid\y_id; %s-au gasit coeficientii theta al polinomului NARX

DKval_predictie = generare_PHI(val,na,nb);
PHIval_predictie = phi_narx(vector_puteri,DKval_predictie,length(val.y));
yhat_predictie = PHIval_predictie*THETA; %echivalent cu predictia cu un pas inainte

val_sim = iddata(yhat_predictie,val.u,val.Ts); 
%generez structura iddata pt. functiile create anterior, continand iesirea deja prezisa
Dk_simulare = generare_PHI(val_sim, na,nb); %vectorii de semnale intarziate, luate din simularea anterioara
PHIsim = phi_narx(vector_puteri,Dk_simulare,length(yhat_predictie)); %generare matrice de regresori
yhat_simulare = PHIsim*THETA;

%close all;
figure(2);
plot(val.y), hold on;
plot(yhat_predictie);
%plot(yhat_simulare);
hold off;
legend("y validare","predictie","simulare");
grid
N = length(yhat_predictie);
MSE_pred = 1/N * sum((val.y-yhat_predictie).^2);
MSE_sim = 1/N * sum((val.y - yhat_simulare).^2);
title(strcat("MSE pred = ",strcat(num2str(MSE_pred),strcat(", MSE sim = ",num2str(MSE_sim)))));



function vector_puteri = combinare_unica(na,nb,m) %returneaza toate combinatiile de puteri posibile ai parametrilor

    v = zeros(1,(na+nb)*(m+1));
    for i = 0:m %popularea vectorului cu numere de la 0 la m
        M = 0;
        while(M <= (na+nb)*m)
            v(1,i+1+M)=i;
%             if(i==0) 
%                 M = 1; %pe prima iterare trebuie ca popularea sa fie corecta
%             end
            M = M + m+1;
        end
%         v(1,i+1+m+1)=i;
%         v(1,i+1+2*(m+1))=i;
    end
    vector_puteri = nchoosek(v, na+nb); %returneaza toate combinatiile posibile de lungime na+nb din vectorul v
    vector_puteri = unique(vector_puteri,'rows'); %returneaza doar liniile unice care nu sunt dublate
    N = length(vector_puteri);  %returneaza numarul de linii
    lungime_noua = N;
    i = 0;
    while(i <= N) %liniile in care suma este mai mare decat m, trebuie sa fie sterse
        i = i + 1;
        suma_puterilor = 0;
        for j = 1:na+nb
            suma_puterilor = suma_puterilor + vector_puteri(i,j); %se insumeaza numerele din coloane
        end
        if(suma_puterilor > m)
            vector_puteri(i,:) = [];
            i = i-1;
            lungime_noua = length(vector_puteri);
        end
        if(i == lungime_noua)
            break;
        end
    end
end

function PHI = generare_PHI(data, na, nb)

    y = data.y;
    u = data.u;

    N = length(y);
    PHI = zeros(N, na+nb);
    for i = 1:N
        for j = 1:na
            if((i-j)>0)
                PHI(i,j) = y(i-j);
            end
        end

        for j = 1:nb
            if((i-j)>0)
                PHI(i, na+j) = u(i-j);
            end
        end
    end

end

function PHI = phi_narx(vector_puteri, d_k, N)
%genereaza matricea regresorilor pt. ARX neliniar
    PHI = ones(N,length(vector_puteri));
    for i = 1:N
        for j = 1:length(d_k(1,:)) %=na+nb
            for k = 1:length(vector_puteri)
                PHI(i,k) = PHI(i,k)*d_k(i,j)^vector_puteri(k,j);
            end
        end
    end

end