# QFT
Quantitative Feedback Theory (QFT) based controller for Floating Offshore Wind Turbines (FOWTs)

## About
This repository is a WIP

## Notes
The ```WWTP_MIMO_2x2_Method2*.m``` files have a problem with the sequential method for the diagonal controller $G_{\beta}(s)$ generation.

Will investigate whenever I get the opportunity.

## Additions to the basic QFT package

### Nyquist stability
Based on Dr. Mario Garcia-Sanz's work:

> <cite>Garcia-Sanz, M. (2016) The Nyquist stability criterion in the Nichols chart. Int. J. Robust. Nonlinear Control, 26: 2643– 2651. [doi: 10.1002/rnc.3465](https://doi.org/10.1002/rnc.3465).</cite>

To use, generate your transfer function and pass it to nyquistStability(). Can be used as follows:

```
% Function to calculate closed-loop system stability according to
%   method introduced by Mario Garcia-Sanz (2016)
% 
% Also included in the QFTCT, Controller design window by Dr. Garcia
% -----------------------------------------------------------------
%
%   INPUTS:-
%       - T_OL      : Open-loop transfer function, i.e. L(s) = P(s) C(s)
%       - DEBUG     : If true, prints extra variables to help debug
%
%   RETURN:-
%       Returned variables are stored within the struct (u)
%       - zc        : = N + num_p_RHP
%       - N         : Total number of encirclements, N = Na + Nb + Nc + Nd
%       - num_p_RHP : Number of poles on the RHP in the s-plane
%       - Na        : Number of encirclements uncer case (Na)
%       - Nb        : Number of encirclements uncer case (Nb)
%       - Nc        : Number of encirclements uncer case (Nc)
%       - Nd        : Number of encirclements uncer case (Nd)
%       - zpCancel  : RHP zero-pole cancelations. "0"=no, "1"=yes
%       - k         : = [0, ±1, ±2, ...] depending on phase at w0 and w1
%       - sigma     : Describes sign of gain of T_OL ( = 0 if gain >= 0, 1
%                       otherwise)
%       - alpha     : = [0, +1, +2, ...] depending on phase at w0, w1, and
%                       w_inf
%       - gamma     : = [-2, 0, +2] depending on type and number of poles
%                       and zeros
%

output = nyquistStability( T_OL, DEBUG );

```
Where,

- **T_OL**&emsp;&nbsp;: Open-loop transfer function, L(s)=G(s)P(s)
- **DEBUG** : Optional flag. If true, prints extra information to the command window

### Disturbance rejection with a feedforward element specification
Created genbnd12().m function. Can be called as follows:


```
% Compute QFT bounds for the following closed-loop configuration
%
%      |  M(s) + P(s)G_f(s)  |
%      | ------------------- | <= del_12(w)
%      |     1 + P(s)G(s)    |

bnd12 = genbnds( 12, w, del_12, M, PG_f, 1, P, P_0 );

% OR
% More generically,

bnd = genbnds( ptype, w, del_12, a, b, c, d, Pnom );
```

Where the unknown is the controller G(s), and

- **ptype** &nbsp;: For disturbance rejection with a feedforward element, use 12
- **w**  &emsp;&emsp; : Working frequency range
- **del_12** : Disturbance rejection specification
- **a**   &emsp; &emsp; : Disturbance on plant output dynamics TF, M(s) 
- **b**   &emsp; &emsp; : Product of plant and feedforward controller TF, P(s)*G_f(s)
- **c**   &emsp; &emsp; : Scalar value of 1
- **d**   &emsp; &emsp; : Plant TF, P(s)
- **P_0**&emsp;&nbsp; : Nominal plant TF, P_0(s)
