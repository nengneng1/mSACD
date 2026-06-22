function curvedistance = signalseparation2fmincon(ng1, ng2,  examplecurve, coeffs)
G1_coeff = coeffs(1);
G2_coeff = coeffs(2);
ng1 = ng1 / max(ng1);
ng2 = ng2 / max(ng2);
examplecurve = examplecurve / max(examplecurve);
combined_curve = G1_coeff * ng1 + G2_coeff * ng2 ;
curve_difference = (combined_curve - examplecurve).^2;
mse = abs(mean(curve_difference));
curvedistance = mse;
end