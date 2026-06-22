function curvedistance = signalseparation3fmincon(ng1, ng2, ng3, examplecurve, coeffs)

% coeffs(coeffs<0)=0;
% coeffs(coeffs>1)=1;

G1_coeff = coeffs(1);
G2_coeff = coeffs(2);
G3_coeff = coeffs(3);
% bg = coeffs(4);
ng1 = ng1 / max(ng1);
ng2 = ng2 / max(ng2);
ng3 = ng3 / max(ng3);
examplecurve = examplecurve / max(examplecurve);
% examplecurve = max(min(examplecurve,1),0);

% combined_curve = bg + G1_coeff * ng1 + G2_coeff * ng2 + G3_coeff * ng3;
combined_curve = G1_coeff * ng1 + G2_coeff * ng2 + G3_coeff * ng3;
% curve_difference = abs(combined_curve - examplecurve);
% curvedistance = sqrt(mean(curve_difference.^2));
curve_difference = (combined_curve - examplecurve).^2;
mse = abs(mean(curve_difference));
curvedistance = mse;
end