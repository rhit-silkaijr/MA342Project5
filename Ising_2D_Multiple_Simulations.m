clc; clear all; close all;

T = linspace(1,4,100);
Et = zeros(length(T),1);
Mt = zeros(length(T),1);
Ct = zeros(length(T),1);
Chit = zeros(length(T),1);
i = 1;
k = 1;
N = 100;

for t=T
    [E, M] = Ising_2D_Func(t, N);
    Et(i) = mean(E);
    Ct(i) = (mean(E.^2) - mean(E)^2)/(k*t^2);
    Mt(i) = mean(abs(M));
    Chit(i) = (mean(abs(M).^2) - mean(abs(M))^2)/(k*t);
    i = i + 1;
end

scatter(T, Et/N);
ylabel('Sample Energy per N');
xlabel('Temperature (T)');

figure

scatter(T, Ct/N);
ylabel('Sample Heat Capacity per N');
xlabel('Temperature (T)');

figure

scatter(T, Mt/N);
ylabel('Sample Abs. Mag. per N');
xlabel('Temperature (T)');

figure

scatter(T, Chit/N);
ylabel('Sample Susceptibility per N');
xlabel('Temperature (T)');