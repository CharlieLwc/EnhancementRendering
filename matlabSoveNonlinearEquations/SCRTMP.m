 
function [] = start(fileName)


clc;
clear;
size = 16;
halfSize = size/2;


lightTopLeftFront=rand(1,halfSize*halfSize);
lightLeftFront = rand(1, halfSize*size);
lightFrontLeft = rand(1, halfSize*size);
lightDownLeftFront=rand(1,halfSize*halfSize);


PartLightA = [lightTopLeftFront, lightLeftFront, lightFrontLeft, lightDownLeftFront];

lightTopLeftFront=rand(1,halfSize*halfSize);
lightLeftFront = rand(1, halfSize*size);
lightFrontLeft = rand(1, halfSize*size);
lightDownLeftFront=rand(1,halfSize*halfSize);


PartLightB = [lightTopLeftFront, lightLeftFront, lightFrontLeft, lightDownLeftFront];

PartLightAB = [PartLightA, PartLightB];

FullLightAB = symmetryLight(PartLightAB);

FullLightAB = FullLightAB*60;
    T='light4.txt';
        dlmwrite(T,FullLightAB,'delimiter', ' ');

%%
function fullLightAB = symmetryLight(PartLightAB)
 
squreSize = sqrt(size(PartLightAB,2)/3);
halfsqureSize = squreSize/2;

startInd = 1;
endInd = halfsqureSize*halfsqureSize;
lightTopLeftFront= PartLightAB(startInd:endInd);
startInd = startInd+halfsqureSize*halfsqureSize;
endInd = endInd+halfsqureSize*squreSize;
lightLeftFront = PartLightAB(startInd:endInd);
startInd = startInd+halfsqureSize*squreSize;
endInd = endInd+halfsqureSize*squreSize;
lightFrontLeft = PartLightAB(startInd:endInd);
startInd = startInd+halfsqureSize*squreSize;
endInd = endInd+halfsqureSize*halfsqureSize;
lightDownLeftFront = PartLightAB(startInd:endInd);

lightTopLeftFront = reshape(lightTopLeftFront, halfsqureSize, halfsqureSize);
lightLeftFront = reshape(lightLeftFront, halfsqureSize, squreSize);
lightFrontLeft = reshape(lightFrontLeft, halfsqureSize, squreSize);
lightDownLeftFront = reshape(lightDownLeftFront, halfsqureSize, halfsqureSize);

lightTopLeftBack = fliplr(lightTopLeftFront);
%lightTopLeftBack(:) = 0;
lightTopLeft = [lightTopLeftBack, lightTopLeftFront];
lightTopRight = flipud(lightTopLeft);
%lightTopRight(:) = 0;
lightTop = [lightTopLeft; lightTopRight];
lightDownLeftBack = fliplr(lightDownLeftFront);
%lightDownLeftBack(:) = 0;
lightDownLeft = [lightDownLeftFront, lightDownLeftBack];
lightDownRight = flipud(lightDownLeft);
%lightDownRight(:) = 0;
lightDown = [lightDownLeft; lightDownRight];
lightLeftBack = flipud(lightLeftFront);
%lightLeftBack(:) = 0;
lightLeft = [lightLeftBack; lightLeftFront];
lightRight = flipud(lightLeft);
%lightRight(:) = 0;
lightFrontRight = flipud(lightFrontLeft);
%lightFrontRight(:) = 0;
lightFront = [lightFrontLeft; lightFrontRight];
lightBack = fliplr(lightFront);

lightTop = reshape(lightTop, 1,squreSize*squreSize);
lightLeft = reshape(lightLeft, 1,squreSize*squreSize);
lightRight = reshape(lightRight, 1,squreSize*squreSize);
lightFront = reshape(lightFront, 1,squreSize*squreSize);
lightDown = reshape(lightDown, 1,squreSize*squreSize);
lightBack = reshape(lightBack, 1,squreSize*squreSize);

lightA = [lightTop, lightLeft, lightFront, lightRight, lightDown, lightBack];

startInd = startInd+halfsqureSize*halfsqureSize;
endInd = endInd+halfsqureSize*halfsqureSize;
lightTopLeftFront= PartLightAB(startInd:endInd);
startInd = startInd+halfsqureSize*halfsqureSize;
endInd = endInd+halfsqureSize*squreSize;
lightLeftFront = PartLightAB(startInd:endInd);
startInd = startInd+halfsqureSize*squreSize;
endInd = endInd+halfsqureSize*squreSize;
lightFrontLeft = PartLightAB(startInd:endInd);
startInd = startInd+halfsqureSize*squreSize;
endInd = endInd+halfsqureSize*halfsqureSize;
lightDownLeftFront = PartLightAB(startInd:endInd);

lightTopLeftFront = reshape(lightTopLeftFront, halfsqureSize, halfsqureSize);
lightLeftFront = reshape(lightLeftFront, halfsqureSize, squreSize);
lightFrontLeft = reshape(lightFrontLeft, halfsqureSize, squreSize);
lightDownLeftFront = reshape(lightDownLeftFront, halfsqureSize, halfsqureSize);


lightTopLeftBack = fliplr(lightTopLeftFront);
%lightTopLeftBack(:) = 0;
lightTopLeft = [lightTopLeftBack, lightTopLeftFront];
lightTopRight = flipud(lightTopLeft);
%lightTopRight(:) = 0;
lightTop = [lightTopLeft; lightTopRight];
lightDownLeftBack = fliplr(lightDownLeftFront);
%lightDownLeftBack(:) = 0;
lightDownLeft = [lightDownLeftFront, lightDownLeftBack];
lightDownRight = flipud(lightDownLeft);
%lightDownRight(:) = 0;
lightDown = [lightDownLeft; lightDownRight];
lightLeftBack = flipud(lightLeftFront);
%lightLeftBack(:) = 0;
lightLeft = [lightLeftBack; lightLeftFront];
lightRight = flipud(lightLeft);
%lightRight(:) = 0;
lightFrontRight = flipud(lightFrontLeft);
%lightFrontRight(:) = 0;
lightFront = [lightFrontLeft; lightFrontRight];
lightBack = fliplr(lightFront);


lightTop = reshape(lightTop, 1,squreSize*squreSize);
lightLeft = reshape(lightLeft, 1,squreSize*squreSize);
lightRight = reshape(lightRight, 1,squreSize*squreSize);
lightFront = reshape(lightFront, 1,squreSize*squreSize);
lightDown = reshape(lightDown, 1,squreSize*squreSize);
lightBack = reshape(lightBack, 1,squreSize*squreSize);

lightB = [lightTop, lightLeft, lightFront, lightRight, lightDown, lightBack];

fullLightAB = [lightA, lightB];
        
%%