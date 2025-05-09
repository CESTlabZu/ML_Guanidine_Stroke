function FitResult = mfit_poly(k, input, fitParam)
    input_background = input(9:29);
    k_background = k(9:29);
    input_background(10:21) = [];
    k_background(10:21) = [];
    FitResult.xindex = linspace(min(input),max(input),length(input_background));
    x0_background = [0,0,0,0];
    lb=[-1, -1, -1,-0.1];
    ub=[10,10,10,0];
    options=optimset('MaxFunEvals',1e6,'TolFun',1e-6,'TolX',1e-6, 'Display',  'off' );
    [FitResult.Coefficents,resnorm]=lsqcurvefit(@CruveFunction_background_poly,x0_background,k_background,input_background,lb,ub,options,fitParam);
    FitResult.Background = CruveFunction_background_poly(FitResult.Coefficents,k(17:31),fitParam);
end

