 function FitResult = mfitpoly_opt(k, input, fitParam)
    index = [1, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 28, 31, 37, 40, 53, 63]+1;
    input_background = input(index(2:13));
    k_background = k(2:13);
    input_background(5:11) = [];
    k_background(5:11) = [];
    FitResult.xindex = linspace(min(input),max(input),length(input_background));
    x0_background = [0,0,0,0];
    lb=[-1, -1, -1,-0.1];
    ub=[10,10,10,0];
    options=optimset('MaxFunEvals',1e6,'TolFun',1e-6,'TolX',1e-6, 'Display',  'off' );
    [FitResult.Coefficents,resnorm]=lsqcurvefit(@CruveFunction_background_poly,x0_background,k_background,input_background,lb,ub,options,fitParam);
    FitResult.Background = CruveFunction_background_poly(FitResult.Coefficents,k(index(2:13)),fitParam);
end

