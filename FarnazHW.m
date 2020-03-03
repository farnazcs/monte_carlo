clc
clear all
close all

%Res=textread('C:\Users\farna\Downloads\MC_result100x100.txt');
Res=textread('C:\Users\farna\Desktop\results.txt');

T=Res(:,1:3:end-1);
N=Res(:,2:3:end-1);
CM=Res(:,3:3:end-1);

[m,n]=size(N);

for i=1:m
    plot(N(i,:),T(i,:));
    hold on
end
xlabel(' X_n');
ylabel('T');
figure
mu=-0.5:0.05:0.5;
for i=1:n
    plot(N(:,i),mu);
    hold on
end
xlabel(' X_n');
ylabel('mu (eV)');
figure
for i=1:m
    plot(T(i,:),CM(i,:),'*');
    hold on
end
xlabel(' T');
ylabel('Heat capacity');
ylim([0,10])


figure
 i=(m+1)/2;
 plot(T(i,:),CM(i,:),'-*','linewidth',2);
xlabel(' T');
ylabel('Heat capacity');