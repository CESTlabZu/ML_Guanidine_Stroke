function theWidth = fwhm2(y,x)
    % Find the maximum value of the spectrum and its half-maximum
    [maxVal, maxIdx] = max(y);
    halfMax = maxVal / 2;
    
    % Find where the signal crosses the half-maximum
    crossingPoints = find(y >= halfMax);
    
    % Find the left and right boundaries of the FWHM
    leftIdx = crossingPoints(1);
    rightIdx = crossingPoints(end);
    
    % Compute the FWHM
    theWidth = x(rightIdx) - x(leftIdx);
end