function [OUT]=GetSamplingOp(IN)
   % Function to obtain sampling operator for missing traces.
   % Syntax:
   %        [OUT]=GetSamplingOp(IN)
   % INPUTS:
   %        IN - Seismic section with missing entries
   %
   % OUTPUTS:
   %        OUT - Sampling Operator
   %
   % Breno Bahia
%%

   % Preliminares
   cutoff=1e-10;
   
   OUT=zeros(size(IN));

   % Get sampling operator
   OUT(sqrt(abs(IN))>cutoff)=1.;
  
end   
