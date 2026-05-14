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

sizeSigma = 100;
%sigma = [1 1 -1 1 1 1 1 -1];
sigma = ones(sizeSigma,1);
for i=1:sizeSigma
    randVar = rand;
    if randVar > 0.6
        sigma(i) = -1;
    end
end
sigmaStart = sigma;

numIters = 500;
reversalProb = .6;
%T = 1;
T = linspace(1,4,20);
A = adjMatrix;
currEnergy = zeros(numIters+1,length(T));
iter = linspace(0,numIters,numIters+1);
for t=1:length(T)
    sigma = sigmaStart;
    currEnergy(1,t) = ESigma(sigma, A);
    
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
    end
end

%scatter(iter, currEnergy(:,3)/sizeSigma);
%xlabel('Simulation iteration (k)');
%ylabel('Energy per lattice print (E/N)');

ETN = zeros(length(T),1);
for i=1:length(T)
    ETN(i) = mean(currEnergy(:,i))/numIters;
end

scatter(T, ETN);