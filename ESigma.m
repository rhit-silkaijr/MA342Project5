function [E] = ESigma(sigma,A)
    E = 0;
    for i = 1:length(A)
        for j=1:i
            if A(i,j) ~= 0
                E = E + sigma(i)*sigma(j);
            end
        end
    end
    E = -E;
end