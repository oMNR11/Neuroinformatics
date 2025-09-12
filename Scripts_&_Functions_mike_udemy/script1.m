%% evaluating code in a script

y=6;
x=5;

%% viewing functions 

z=10;
% edit linspace;

% some functions take multiple inputs
linspace(1,3,5)

%function outputs stored as variables
numsum= sum([1 2 3]);

%some functions take multiple outputs
[maxval,maxindex] = max([1 4 2 2]);

%creating a matrix of random numbers
randmat = randn(10,4,5);

randvect = randn(10); %creates a matrix of 10x10 size with random elements


%%  Section     8   exercises

%create a 4x6 matrix of random integers between 5 and 17
matA = zeros(4,6);
for rowi = 1:4
    for colj = 1:6
            matA(rowi,colj) = randi([5,17]);
    end
end

oddCount = sum(mod(matA(:), 2) == 1);
evenCount = numel(matA) - oddCount;
proportionOdd = oddCount / numel(matA);
proportionEven = evenCount / numel(matA);

%report the proportion of odd and even numbers 
    oddcount= 0;
    evenCount = 0;
for rowi = 1:4
    for colj = 1:6
           if (mod(matA(rowi,colj),2) == 1)
                oddcount = oddcount + 1;
           else
               evenCount = evenCount + 1;
           end
    end
end
