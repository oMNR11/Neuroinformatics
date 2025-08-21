function stats  =   basicstats(invar)

% This function computes the mean, sun and standard deviation of an input 

%   INPUT : invar - input   vector  ,   must    be  numeric atleast 3
%   elements    long.

%   OUPTUT  :   stats - a   3-element   vector  containing  [sum    mean
%   sd]

%% input checks 

%check that the input is numeric 

% input has atleast 3 numbers 

%   must    be  a   vector  ,   not a   matrix  

%% compute  the basic   statistics

%   compute the sum 
    thesum  =   sum(invar);

%   compute the mean
    themean =   mean(invar);
%   compute the standard    deviation
    thestd  =   std(invar);
%%  display the results to  the command window  
    disp(['The sum of the input is ' num2str(thesum) '.'])
    disp(['The mean of the input is ' num2str(themean) '.'])
    disp(['The standard deviation of the input is ' num2str(thestd) '.'])

%%  prepare the output

stats = [thesum themean thestd];