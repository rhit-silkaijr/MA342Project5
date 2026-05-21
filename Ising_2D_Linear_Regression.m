T = linspace(2.15,2.3,15);
k = 1;
N = linspace(4,20,17);

M = zeros(17,1);
C = zeros(17,1);
Chi = zeros(17,1);

for n=N

    i = 1;

    Et = zeros(length(T),1);
    Mt = zeros(length(T),1);
    Ct = zeros(length(T),1);
    Chit = zeros(length(T),1);

    nPadded = n + 2;

    % Set up initial randomized A
    A = rand(nPadded, nPadded);
    A(A>0.5) = 1;
    A(A<=0.5) = -1;

    for t=T
        [E, m, A] = Ising_2D_Func(t, A, i==1, 500, 1);
        Et(i) = mean(E);
        Ct(i) = (mean(E.^2) - mean(E)^2)/(k*t^2);
        Mt(i) = mean(abs(m));
        Chit(i) = (mean(abs(m).^2) - mean(abs(m))^2)/(k*t);
        i = i + 1;
    end

    [~, maxT] = max(Ct);

    M(N==n) = Mt(maxT);
    C(N==n) = Ct(maxT);
    Chi(N==n) = Chit(maxT);

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