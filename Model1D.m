%Set up adjMatrix for A
adjMatrix = zeros(8,8);
for i=1:8
    if i == 1
        adjMatrix(i,8) = 1;
    else
        adjMatrix(i,i-1) = 1;
    end
    if i == 8
        adjMatrix(i,1) = 1;
    else
        adjMatrix(i,i+1) = 1;
    end
end

%1D spin config

sizeSigma = 200;
sigma = ones(sizeSigma,1);
for i=1:sizeSigma
    randVar = rand;
    if randVar > 0.6
        sigma(i) = -1;
    end
end

numIters = 1000;
reversalProb = .25;
T = linspace(1,4,200);
A = adjMatrix;
currEnergy = zeros(numIters+1,length(T));
currMag = zeros(numIters+1,length(T));
iter = linspace(0,numIters,numIters+1);
for t=1:length(T)
    currEnergy(1,t) = ESigma(sigma, A);
    currMag(1,t) = sum(sigma);
    
    for i=1:numIters
        sigmaEnd = sigma;
    
        reversalCandidates = rand(sizeSigma,1);
        reversalIndex = [];
        for k=1:length(sigma)
            if reversalCandidates(k) > reversalProb
                reversalIndex = [reversalIndex k];
            end
        end
    
        for k=1:length(reversalIndex)
            j = reversalIndex(k);
            sigmaFlip = sigma;
            sigmaFlip(j) = sigmaFlip(j)*-1;
            P = ProbReverse(sigma, sigmaFlip, A, t);
            if P > rand
                sigmaEnd(j) = sigma(j)*-1;
            end
        end
        sigma = sigmaEnd;
        currEnergy(i+1,t) = ESigma(sigmaEnd, A);
        currMag(i+1,t) = sum(sigmaEnd);
    end
end

%scatter(iter, currEnergy(:,3)/sizeSigma);
%xlabel('Simulation iteration (k)');
%ylabel('Energy per lattice print (E/N)');

ETN = zeros(length(T),1);
CTN = zeros(length(T),1);
MTN = zeros(length(T),1);
CHITN = zeros(length(T),1);
for i=1:length(T)
    ET = currEnergy(:,i);
    ETN(i) = mean(ET)/numIters;
    CTN(i) = (mean(ET.^2)-mean(ET)^2)/((i^2)*numIters);
    MT = abs(currMag(:,i));
    MTN(i) = mean(MT)/numIters;
    CHITN(i) = (mean(MT.^2)-mean(MT)^2)/(i*numIters);
end
ETN = -ETN(2:length(ETN));

scatter(T(2:length(T)), ETN);
ylabel('Sample Energy per N');
xlabel('Temperature (T)');
%scatter(T, CTN);
%ylabel('Sample Heat Capacity per N');
%xlabel('Temperature (T)');
%scatter(T, MTN);
%ylabel('Sample Abs. Mag. per N');
%xlabel('Temperature (T)');
%scatter(T, CHITN);
%ylabel('Sample Susceptibility per N');
%xlabel('Temperature (T)');


% Now do it for different sizes of sigma
% and use the results for a lin reg
nspace = linspace(4,20,17);
for n=nspace
    sizeSigma = n;
    sigma = ones(sizeSigma,1);
    for i=1:sizeSigma
        randVar = rand;
        if randVar > 0.6
            sigma(i) = -1;
        end
    end
    sigmaStart = sigma;
end
