function FitResult = poly_opt(k,input,fitParam)
    input_clean = smoothdata(input, "sgolay",6);
    index = [1, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 28, 31, 37, 40, 53, 63]+1;
    input_background = input_clean(index(4:14));
    k_background = k(index(4:14));
    input_background(4:10) = [];
    k_background(4:10) = [];
    x0_background = [0,0,0,0];
    lb=[-1, -1, -1,-1];
    ub=[100,100,100,100];
    options=optimset('MaxFunEvals',1e6,'TolFun',1e-6,'TolX',1e-6, 'Display',  'off' );
    [FitResult.Coefficents,resnorm]=lsqcurvefit(@CruveFunction_background_poly,x0_background,k_background,input_background,lb,ub,options,fitParam);
    FitResult.Background = CruveFunction_background_poly(FitResult.Coefficents,k(index(4:13)),fitParam);
    FitResult.Arex = ((1./input(index(4:13))) - (1./FitResult.Background)).*(fitParam.R1).*(1+fitParam.fm_m);
end

