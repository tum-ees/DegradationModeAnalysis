function OCV_Calc = calculate_OCV_blend(x, myData)
%> Author: Can Korkmaz (can.korkmaz@tum.de)
%> supervised by Mathias Rehm (mathias.rehm@tum.de)
%> Date: 2025-06-18
%
%  FUNCTION: calculate_OCV_blend
%  --------------------------------------------------------------------
%  Builds a full-cell open-circuit-voltage (OCV) curve by combining
%  (i) a parameterised anode potential — which can optionally blend a
%  Blend1/Blend2 sub-curve — and (ii) a parameterised cathode potential.
%
%  INPUTS
%      x      – parameter vector
%               [α_an, β_an, α_cat, β_cat, (γ_Blend2)] where:
%                 α,β  : stretch / shift factors for SOC axis
%                 γ_Blend2 : Blend2 fraction in anode blend (0‒1, optional)
%      myData  – struct with reference half-cell data fields
%               .Q_cell            cell charge grid for final OCV
%               .commonVoltage      anode reference curve  (graphite only)
%               .Q_Blend2_interp    SOC grid for anode reference
%               .normCathode_SOC    SOC grid for cathode reference
%               .normCathode_U      cathode reference voltage
%
%               *If γ_Blend2 is given*  (blend mode):
%               .blendAnodeCurve()  function handle OR local function
%                                   returning [socBlend, uBlend]
%
%  OUTPUT
%      OCV_Calc – vector of full-cell voltage values evaluated on
%                 myData.Q_cell
%
%  ASSUMPTIONS / CAVEATS
%      • If x has only 4 entries, γ_Blend2 is assumed 0 → Blend1-only anode.
%      • Out-of-range SOC points in interp1 are extrapolated with 0 V
%        (linear method + ‘linear,0’ extrapolation). Consider replacing
%        the ‘0’ fill with the nearest boundary value if artefacts appear.
%      • No guards against α,β that warp SOC outside [0,1]; caller should
%        validate feasible parameter bounds.
%
    % -------------------- 1.   Normalise the parameter vector ------------
    % Accept both row and column-vector input
    if size(x,1) > 1
        x = x.';                           % force row orientation
    end

    alpha_an  = x(1);
    beta_an   = x(2);
    alpha_cat = x(3);
    beta_cat  = x(4);

    % ------------------- 2.   Determine blend mode / γ_Blend2 ---------------
    if numel(x) < 5                      % default: no blend
        useBlend = false;
    elseif numel (x) < 6
        gamma_Blend2 = x(5);
        inhomo_mag_an = 0;
        inhomo_mag_ca = 0;
        useBlend = true;
    else
        gamma_Blend2 = x(5);
        inhomo_mag_an = x(6);
        inhomo_mag_ca = x(7);
        useBlend = true;
    end

   

    % ------------------- 3.   Construct anode potential -----------------
    if useBlend
        % Expect helper that returns blended anode curve for given γ_Si
        [blendSOC, blendU] = blend_anode_curve(gamma_Blend2, myData, inhomo_mag_an);

        % SOC axis stretch/shift → interpolate onto cell charge grid
        anodePot = interp1( alpha_an * blendSOC + beta_an, ...
                            blendU, myData.Q_cell, ...
                            'linear', 0 );
    else
        % Pure graphite reference
        anodePot = interp1( alpha_an * myData.Q_Si_interp + beta_an, ...
                            myData.commonVoltage, myData.Q_cell, ...
                            'linear', 0 );
    end

    
    % ------------------- 4.   Construct cathode potential ---------------
    % calculate inhomogeneities
    if ~(inhomo_mag_ca == 0)
        myData.normCathode_U = calculate_inhomogeneity(myData.normCathode_SOC, myData.normCathode_U, inhomo_mag_ca);
    end

    cathPot = interp1( alpha_cat * myData.normCathode_SOC + beta_cat, ...
                       myData.normCathode_U, myData.Q_cell, ...
                       'linear', 0 );

    % ------------------- 5.   Full-cell OCV calculation -----------------
    OCV_Calc = cathPot - anodePot;
end
