function ReducedSystem = reduceOrder_balred(sys, order, SepTol, PLOT)
% Reduce LTI model order using balanced truncation
%

arguments
    sys     (:,:,:) {mustBeA(sys,'tf')}         = false
    order   (1,1)   {mustBeNumeric}             = -1
    SepTol  (1,1)   {mustBeNumeric}             = 10.0
    PLOT            {mustBeNumericOrLogical}    = false
end

System = sys; % Define System to reduce
Order = order;
 
% Create option set for balred command
Options = balredOptions();
% % Frequency range for computing state contributions
% Options.FreqIntervals = [FreqLO FreqHI];
% Accuracy loss factor for stable/unstable decomposition
Options.SepTol = SepTol;
 
% Compute reduced order approximation on specified frequency range
ReducedSystem = balred(System,Order,Options);
 
% Create comparison plot
if( PLOT )
    figure();
    bode(System,ReducedSystem);
end
end