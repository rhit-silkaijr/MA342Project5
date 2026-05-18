T = linspace(2.15,2.3,10);
i = 1;
k = 1;
N = linspace(4,20,17);

M = zeros(17,1);
C = zeros(17,1);
Chi = zeros(17,1);

for n=N

    Et = zeros(length(T),1);
    Mt = zeros(length(T),1);
    Ct = zeros(length(T),1);
    Chit = zeros(length(T),1);

    for t=T
        [E, m] = Ising_2D_Func(t, n);
        Et(i) = mean(E);
        Ct(i) = (mean(E.^2) - mean(E)^2)/(k*t^2);
        Mt(i) = mean(abs(m));
        Chit(i) = (mean(abs(m).^2) - mean(abs(m))^2)/(k*t);
        i = i + 1;
    end

    M(N==n) = mean(Mt);
    C(N==n) = mean(Ct);
    Chi(N==n) = mean(Chit);

end

x = log(N);
yAlpha = log(C');
yBeta = log(M');
yGamma = log(Chi');

[pAlpha,~] = polyfit(x,yAlpha,1);
[pBeta,~] = polyfit(x,yBeta,1);
[pGamma,~] = polyfit(x,yGamma,1);

xTemp = 1.2:.05:3.2;
alphaLinear = polyval(pAlpha, xTemp);
betaLinear = polyval(pBeta, xTemp);
gammaLinear = polyval(pGamma, xTemp);


scatter(x, yAlpha);
hold on
plot(xTemp, alphaLinear,'-');
ylabel('ln(Heat Cap.)');
xlabel('ln(N)');

figure

scatter(x, yBeta);
hold on
plot(xTemp, betaLinear,'-');
ylabel('ln(Magnetism)');
xlabel('ln(N)');


figure

scatter(x, yGamma);
hold on
plot(xTemp, gammaLinear,'-');
ylabel('ln(Susceptibility)');
xlabel('ln(N)');