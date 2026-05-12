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

sizeSigma = 300;
%sigma = [1 1 -1 1 1 1 1 -1];
sigma = ones(sizeSigma);
for i=1:sizeSigma
    if rand > 0.7
        sigma(i) = -1;
    end
end

numIters = 500;
reversalProb = .1;
T = 2;
A = adjMatrix;
currEnergy = zeros(numIters+1);
currEnergy(1) = ESigma(sigma, A);
iter = linspace(0,numIters,numIters+1);

for i=1:numIters
    sigmaEnd = sigma;
    for j=1:length(sigma)
        sigmaFlip = sigma;
        sigmaFlip(j) = sigmaFlip(j)*-1;
        P = ProbReverse(sigma, sigmaFlip, A, T);
        if P > reversalProb
            sigmaEnd(j) = sigma(j)*-1;
        end
    end
    sigma = sigmaEnd;
    currEnergy(i+1) = ESigma(sigmaEnd, A);
end

scatter(iter, currEnergy)