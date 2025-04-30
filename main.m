clc;clear;close all;

rawImg = imread('SeqImg.tif');
load('PSF.mat');
runCode(rawImg, PSF);



%% Functions

function runCode(rawImg, PSF)

    rawImg = im2double(rawImg);
    

    iter = 20;
    resImg = deconvlucy(rawImg, PSF, iter);

    figure;
    subplot(1,2,1);imshow(rawImg);title('Raw Image');
    subplot(1,2,2);imshow(resImg);title('Restored Image')

    [rawPowerX, rawPowerY, rawMAE] = Valuation(rawImg);
    [resPowerX, resPowerY, resMAE] = Valuation(resImg);

    fprintf('The MAE between the two power spectrum projection curves for raw image is: %f\n', rawMAE);
    fprintf('The MAE between the two power spectrum projection curves for restored image is: %f\n', resMAE);
    
    ccX = corr(rawPowerX, resPowerX);
    ccY = corr(rawPowerY, resPowerY);
    fprintf('The CC between the power spectrum projection curves of the raw and Restored images in the TDI direction is: %f\n', ccX);
    fprintf('The CC between the power spectrum projection curves of the raw and Restored images in the Cross TDI direction is: %f\n', ccY);

end

function [Power_x, Power_y, MAE] = Valuation(Img)

    [row, ~] = size(Img);
    Coordinate = linspace(-1, 1, row);
    % FFT
    Img_fftshift = fftshift(fft2(Img));
    Img_fftshift_Power = real(Img_fftshift).^2 + imag(Img_fftshift).^2; 
    Img_fftshift_log = log(Img_fftshift_Power);
    Img_fftshift_m = Img_fftshift_log - min(Img_fftshift_log(:));
    Img_fftshift_uint8 = uint8(Img_fftshift_m ./ max(Img_fftshift_m(:)) .* 255);

    % projection
    Power_x = sum(Img_fftshift_uint8, 2);
    Power_x = Power_x'; % TDI direction
    Power_x = Power_x ./ row;
    Power_y = sum(Img_fftshift_uint8, 1); % cross TDI direction
    Power_y = Power_y ./ row;

    MAE = sum(abs(Power_x - Power_y)) / row;

    figure;
    scatter(Coordinate, Power_x, 'o', "MarkerEdgeColor","red");
    hold on
    scatter(Coordinate, Power_y, 'h', "MarkerEdgeColor","black");
    legend('TDI direction', 'cross TDI direction')

    Power_x = Power_x';
    Power_y = Power_y';
    
end

