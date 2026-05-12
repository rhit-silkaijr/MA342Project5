function [P] = ProbReverse(sigmaK, sigmaBar, A, T)
    ESigK = ESigma(sigmaK, A);
    ESigBar = ESigma(sigmaBar, A);
    if ESigBar < ESigK
        P = 1;
    else
        P = exp((ESigK-ESigBar)/T);
    end
end