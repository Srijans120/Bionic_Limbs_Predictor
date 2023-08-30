
clc
clear all
close all
files={'101.txt','102.txt','103.txt','104.txt','105.txt','106.txt','107.txt','108.txt','109.txt','110.txt','111.txt','112.txt','113.txt','114.txt','115.txt','116.txt','117.txt','118.txt','119.txt','120.txt','201.txt','202.txt','203.txt','204.txt','205.txt','206.txt','207.txt','208.txt','209.txt','210.txt','211.txt','212.txt','213.txt','214.txt','215.txt','216.txt','217.txt','218.txt','219.txt','220.txt','301.txt','302.txt','303.txt','304.txt','305.txt','306.txt','307.txt','308.txt','309.txt','310.txt','311.txt','312.txt','313.txt','314.txt','315.txt','316.txt','317.txt','318.txt','319.txt','320.txt','401.txt','402.txt','403.txt','404.txt','405.txt','406.txt','407.txt','408.txt','409.txt','410.txt','411.txt','412.txt','413.txt','414.txt','415.txt','416.txt','417.txt','418.txt','419.txt','420.txt'};
A=cell(1,80);
for i=1:80
    file=fopen(files{i});
    A{i}=fscanf(file,'%d',[8,Inf]);
end
N = 655;
TDS = zeros(80,15,4,8);
for i = 1:80
    for j = 1:15
        for k = 1:4
            for l = 1:8
                temp = 0;
                for m = 1:N
                   temp= temp + (s(A,l,j,i,m)-TDS(i,j,1,l))^k;
                end
                TDS(i,j,k,l) = temp/N;
            end
        end
    end
end
ICS = zeros(80,15,12);
for i = 1:80
    for j = 1:15
        index = 0;
        for k = 1:4
            for l = 1:k-1
                temp = 0;
                for m = 1:N
                    temp = s(A,l,j,i,m)*s(A,k,j,i,m);
                end
                index = index+1;
                ICS(i,j,index) = temp/N;
            end
        end

        for k = 5:8
            for l = 5:k-1
                temp = 0;
                for m = 1:N
                    temp = s(A,l,j,i,m)*s(A,k,j,i,m);
                end
                index = index+1;
                ICS(i,j,index) = temp/N;
            end
        end
    end
end

psi = zeros(80,15,8,N);
for i = 1:80
    for j = 1:15
        for k = 1:8
            psi(i,j,k,:) = abs(fft(s(A,k,j,i,1:N))).^2;
        end
    end
end

psiprime = zeros(80,15,8,7);
for i = 1:80
    for j = 1:15
        for k = 1:8
            for l = 1:7
                temp = 0;
                for m = 1:N
                    temp = temp + psi(i,j,k,m)*(m^(l-1));
                end
                psiprime(i,j,k,l) = sqrt(temp);
            end
        end
    end
end

psipro = zeros(80,15,8,17);
for i = 1:80
    for j = 1:15
        for k = 1:8
            % ho = 1;
            psipro(i,j,k,1) = log(psiprime(i,j,k,1));
            psipro(i,j,k,2) = log(psiprime(i,j,k,3));
            psipro(i,j,k,3) = log(psiprime(i,j,k,5));
            psipro(i,j,k,4) = log(psiprime(i,j,k,1)) - .5*log((psiprime(i,j,k,1) - psiprime(i,j,k,3))/(psiprime(i,j,k,1) - psiprime(i,j,k,5)));
            psipro(i,j,k,5) = psipro(i,j,k,3) -.5*(psipro(i,j,k,1) + psipro(i,j,k,5));
            psipro(i,j,k,6) = psipro(i,j,k,1) - .25*log(psiprime(i,j,k,2)*psiprime(i,j,k,4));
            psipro(i,j,k,7) = psipro(i,j,k,1) - .25*log(psiprime(i,j,k,3)*psiprime(i,j,k,7));    
            iter = 8;
            for l = 2:5
                for m = 1:l-1
                    psipro(i,j,k,iter) = .5*log(psiprime(i,j,k,l)*psiprime(i,j,k,m));
                    iter = iter + 1;
                end
            end
        end
    end
end
% function hola = H(4,8)
spb=zeros(80,15,8,10);
for i=1:80
    for j = 1:15
        for k = 1:8
            x=zeros(1,600);
            for l=1:600
                x(l)=s(A,k,j,i,l);
            end
            p=pburg(x,10);
            for o =1:10
                spb(i,j,k,o)=p(o);
            end
        end
    end
end 

lbp = zeros(80,15,8,2);
for i =1:80
    for j = 1:15
        for k = 1:8
            ho = 0;
            ind = 0;
            for l = 1:5
                ho = ho + s(A,k,j,i,l);
            end
            ho = .2*ho;
            for l = 1:N
                ind = ind+ heaviside(s(A,k,j,i,l) - ho);
            end
            lbp(i,j,k,1) = ind;
            lbp(i,j,k,2) = N - ind;
        end
    end
end
features = zeros(1200,276);
for i = 0:79
    for j = 1:15
        for k = 0:3
            for l = 1:8
                features(i*15 + j,k*8 + l) = TDS(i+1,j,k+1,l);
            end
        end
    end
end

for i = 0:79
    for j = 1:15
        for k = 1:12
            features(i*15 + j,32 + k) = ICS(i+1,j,k);
        end
    end
end

for i = 0:79
    for j = 1:15
        for k = 0:7
            for l = 1:17
                features(i*15 + j,k*17 + 44 + l) = psipro(i+1,j,k+1,l);
            end
        end
    end
end

for i = 0:79
    for j = 1:15
        for k = 0:7
            for l = 1:2
                features(i*15 + j,k*2 + 180 + l) = lbp(i+1,j,k+1,l);
            end
        end
    end
end

for i = 0:79
    for j = 1:15
        for k = 0:7
            for l = 1:10
                features(i*15 + j,k*17 + 196 + l) = spb(i+1,j,k+1,l);
            end
        end
    end
end


% svm.train()

function my = s(A,m,w,p,l)
         my=A{p}(m,(w-1)*600 + l);
end

