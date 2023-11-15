function ReducedSystem = reduceOrder(sys, FreqLO, FreqHI, SepTol, PLOT)
% Reduce LTI model order using mode selection
%

arguments
    sys     (:,:,:) {mustBeA(sys,'tf')}         = false
    FreqLO  (1,1)   {mustBeNumeric}             = 0.01
    FreqHI  (1,1)   {mustBeNumeric}             = 10.0
    SepTol  (1,1)   {mustBeNumeric}             = 10.0
    PLOT            {mustBeNumericOrLogical}    = false
end

System = sys; % Define System to reduce
UpperCutoffFrequency = FreqHI;
LowerCutoffFrequency = FreqLO;
 
% Create option set for freqsep command
Options = freqsepOptions();
% Accuracy loss factor for stable/unstable decomposition
Options.SepTol = SepTol;
 
% Select modes between lower and upper cutoff frequencies
ReducedSystem = freqsep(System,UpperCutoffFrequency,Options);
[~,ReducedSystem] = freqsep(ReducedSystem,LowerCutoffFrequency,Options);
 
% Create comparison plot
if( PLOT )
    figure()
    bode(System,ReducedSystem);
    legend( "Original", "Reduced System" );
end
end