function FitResult = poly(k,input,fitParam)
    input_background = input(17:29); % index 17:29 corresponds to 3 ppm to 1 ppm
    k_background = k(17:29);
    input_background(4:10) = []; % index 4:10 corresponds to 2.5 ppm to 1.5 ppm
    k_background(4:10) = [];
    x0_background = [0,0,0,0];
    lb=[-1, -1, -1,-1];
    ub=[100,100,100,100];
    options=optimset('MaxFunEvals',1e6,'TolFun',1e-6,'TolX',1e-6, 'Display',  'off' );
    [FitResult.Coefficents,resnorm]=lsqcurvefit(@CruveFunction_background_poly,x0_background,k_background,input_background,lb,ub,options,fitParam);
    FitResult.Background = CruveFunction_background_poly(FitResult.Coefficents,k,fitParam);
    FitResult.Arex = ((1./input) - (1./FitResult.Background)).*(fitParam.R1).*(1+fitParam.fm_m);
end

