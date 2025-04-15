function y = spectrum(amp, width, off)
    max=1500;
    step=50;
    offset= -max:step:max;
    rffreq  =[-4000, -3500, -3000, -2500, offset, 2500, 3000,3500,4000];
    y = (amp)./(1+(((rffreq-off).^2)./((0.5*width).^2)));
end
