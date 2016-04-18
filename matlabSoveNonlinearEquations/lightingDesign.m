
%% solve function
function [] = Color2Grey(fileName)
clc;

global LitToColLR;
global LineIndex;
global isSymmetryXZ;
global isRidgeSame;

global lightNum;
global featureNum;


global Diff;
global ColorLRlabA;
global ColorLRlabB;
global SameLine;

global isScaled;
global Scale;


ColorLRlabA = ones(10);
ColorLRlabB = ones(10);
Diff = ones(10);

isSymmetryXZ = false;
isRidgeSame = false;

isScaled = false;

FileName = 'building_cross';
FlieAfter = '_scale.obj';
FullFileName = [FileName,FlieAfter];



Scale = load(FullFileName)';

Scale(:,1) = Scale(:,1)/norm(Scale(:,1));
Scale(:,2) = Scale(:,2)/norm(Scale(:,2));
Scale(:,3) = Scale(:,3)/norm(Scale(:,3));

Scale(:,4) = sum(Scale(:,:))*0.3333;
Scale(:,5) = cross(Scale(:,2), Scale(:,1));
Scale(:,6) = cross(Scale(:,3), Scale(:,2));
Scale(:,7) = cross(Scale(:,1), Scale(:,3));

Scale(:,4) = Scale(:,4)/norm(Scale(:,4));
Scale(:,5) = Scale(:,5)/norm(Scale(:,5));
Scale(:,6) = Scale(:,6)/norm(Scale(:,6));
Scale(:,7) = Scale(:,7)/norm(Scale(:,7));





%TFileName = 'flubberNN_FeatureToLight_16_1309.obj';
%LFileName = 'flubberNN_PlaNorLines.obj';
%TFileName = 'BolNHNN_FeatureToLight_8_2032.obj';
TFileName = 'cubeNN_FeatureToLight_8_887.obj';



LitToColLR = load(TFileName);
LitToColLR = LitToColLR';
if(isRidgeSame)
    LineIndex = load(LFileName);
end;

lightNum = size(LitToColLR, 1);
if(isSymmetryXZ)
    lightNum = lightNum /4;
end;
featureNum = size(LitToColLR, 2)/2;

LightABOpt = optimLight();

saveLight(LightABOpt, FileName);


%%
 function LightABOpt = optimLight()
global lightNum;     

UpperBound = ones(1,lightNum*3);
LowerBound = zeros(1,lightNum*3);

%particlesWarm
opts = optimoptions('particleswarm','plotFcns', @plotFunc,...
    'Display','iter');
 
[LightABOpt, fval, eflag, output] = particleswarm(@my_fun,...
    lightNum*3,LowerBound,UpperBound , opts);
fval
eflag
output
%% plorFunction
function stop = plotFunc(LightStrut, optimValues, state, varargin)
stop = false;
global Diff;
global lightNum;
global ColorLRlabA;
global ColorLRlabB;

global isSymmetryXZ;

global isScaled;
LightRGB = LightStrut.bestx;



if(isSymmetryXZ)
    LightRGB = symmetryLightXZ(LightRGB, 3);
end;

LightR = LightRGB(1:lightNum);
LightG = LightRGB(lightNum+1:lightNum+lightNum);
LightB = LightRGB(lightNum+lightNum+1:lightNum+lightNum+lightNum);


if(isScaled)
    [LightR, LightG, LightB] = Project(LightR, LightG, LightB);
end;




subplot(2,3,1);
scatter3(LightR,LightG,LightB);

subplot(2,3,2);
plot(ColorLRlabA,ColorLRlabB,'k.');

subplot(2,3,3);
hist(Diff);

subplot(2,3,4);
plot(LightR,LightG,'k.');

subplot(2,3,5);
plot(LightG,LightB,'k.');

subplot(2,3,6);
plot(LightR,LightB,'k.');


global isRidgeSame;
if(isRidgeSame)
    global SameLine;
    subplot(2,2,4);
    bar(SameLine);
end;
%% main function
function sumDiff = my_fun(LightRGB)

global LitToColLR;
global isSymmetryXZ;
global isRidgeSame;


global Diff;
global lightNum;
global featureNum;


global ColorLRlabA;
global ColorLRlabB;


global isScaled;


if(isSymmetryXZ)
    LightRGB = symmetryLightXZ(LightRGB, 3);
end;

LightR = LightRGB(1:lightNum);
LightG = LightRGB(lightNum+1:lightNum+lightNum);
LightB = LightRGB(lightNum+lightNum+1:lightNum+lightNum+lightNum);

if(isScaled)
    [LightR, LightG, LightB] = Project(LightR, LightG, LightB);
end;


ColorLRR = LightR * LitToColLR;
ColorLRG = LightG * LitToColLR;
ColorLRB = LightB * LitToColLR;

[ColorLRR, ColorLRG, ColorLRB, translate, gamma] =...
    uniform(ColorLRR,ColorLRG, ColorLRB, 1, 0);
    
[ColorLRlabA, ColorLRlabB] = ColorRGB2LAB(ColorLRR, ColorLRG, ColorLRB);

DiffA = ColorLRlabA(1:featureNum) - ColorLRlabA(featureNum+1:featureNum+featureNum);
DiffB = ColorLRlabB(1:featureNum) - ColorLRlabB(featureNum+1:featureNum+featureNum);

colorBound = 20;

Diff = sqrt(DiffA.*DiffA+DiffB.*DiffB);

if(isRidgeSame)
    SameLR = zeros(1,featureNum);
    global LineIndex;
    global SameLine;
    lineNum = size(LineIndex, 1);
    SameLine = zeros(1,lineNum);
    for i=1:lineNum
        StarEnd = LineIndex(i,:);
        PartLRA = ColorLRlabA(StarEnd(1):StarEnd(2));
        PartLRB = ColorLRlabB(StarEnd(1):StarEnd(2));
        PartLRAB = [PartLRA;PartLRB]';
        Diss = pdist(PartLRAB,'euclidean');
        Diss = Diss.*Diss;
        dissSqr = sum(Diss)/(size(Diss,2) + (size(Diss,2) == 0));
        SameLR(StarEnd(1):StarEnd(2)) = dissSqr;
        SameLine(1,i) = dissSqr;
        PartLRA = ColorLRlabA(StarEnd(1)+featureNum:StarEnd(2)+featureNum);
        PartLRB = ColorLRlabB(StarEnd(1)+featureNum:StarEnd(2)+featureNum);
        PartLRAB = [PartLRA;PartLRB]';
        Diss = pdist(PartLRAB,'euclidean');
        Diss = Diss.*Diss;
        dissSqr = sum(Diss)/(size(Diss,2) + (size(Diss,2) == 0));
        SameLR(StarEnd(1):StarEnd(2)) = SameLR(StarEnd(1):StarEnd(2)) + dissSqr;
        SameLine(1,i) = dissSqr;
    end;
    SameLine = SameLine*0.5;
    Diff2 = min(Diff, colorBound);
    Diff2 = colorBound - Diff2 + SameLR.*SameLR;
else
    Diff2 = min(Diff, colorBound);
    Diff2 = colorBound - Diff2;
end;

sumDiff = sum(Diff2.*Diff2);
%%
function fullLight = symmetryLightXZ(PartLight, channelNum)
 
global lightNum;

fullLight = zeros(lightNum*channelNum);

squreSize = sqrt(size(PartLight,2)*4/6/channelNum);
halfsqureSize = squreSize/2;
startInd = 1;
endInd = halfsqureSize*halfsqureSize;

lightStart = 1;
lightEnd = lightNum;
    
for channel = 1:channelNum;
    
    lightTopLeftFront= PartLight(startInd:endInd);
    startInd = startInd+halfsqureSize*halfsqureSize;
    endInd = endInd+halfsqureSize*squreSize;
    lightLeftFront = PartLight(startInd:endInd);
    startInd = startInd+halfsqureSize*squreSize;
    endInd = endInd+halfsqureSize*squreSize;
    lightFrontLeft = PartLight(startInd:endInd);
    startInd = startInd+halfsqureSize*squreSize;
    endInd = endInd+halfsqureSize*halfsqureSize;
    lightDownLeftFront = PartLight(startInd:endInd);
    startInd = startInd+halfsqureSize*halfsqureSize;
    endInd = endInd+halfsqureSize*halfsqureSize;


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
    
    
    fullLight(lightStart:lightEnd) = [lightTop, lightLeft, lightFront,...
        lightRight, lightDown, lightBack];
    lightStart = lightStart+lightNum;
    lightEnd = lightEnd+lightNum;
    
end

%% 
function [] = saveLight(LightRGB, T)

global LitToColLR;

global isSymmetryXZ;

global lightNum;
global featureNum;

global isScaled;
if(isSymmetryXZ)
    LightRGB = symmetryLightXZ(LightRGB, 3);
end;

LightR = LightRGB(1:lightNum);
LightG = LightRGB(lightNum+1:lightNum+lightNum);
LightB = LightRGB(lightNum+lightNum+1:lightNum+lightNum+lightNum);


if(isScaled)
    [LightR, LightG, LightB] = Project(LightR, LightG, LightB);
end;

LightRGB = [LightR, LightG, LightB];

subplot(2,2,1);
scatter3(LightR,LightG,LightB);

ColorLRR = LightR * LitToColLR;
ColorLRG = LightG * LitToColLR;
ColorLRB = LightB * LitToColLR;

[ColorLRR, ColorLRG, ColorLRB, translate, gamma] =...
    uniform(ColorLRR,ColorLRG, ColorLRB, 1, 0);
    
[ColorLRlabA, ColorLRlabB] = ColorRGB2LAB(ColorLRR, ColorLRG, ColorLRB);

subplot(2,2,2);
plot(ColorLRlabA,ColorLRlabB,'k.');


DiffA = ColorLRlabA(1:featureNum) - ColorLRlabA(featureNum+1:featureNum+featureNum);
DiffB = ColorLRlabB(1:featureNum) - ColorLRlabB(featureNum+1:featureNum+featureNum);
Diff = sqrt(DiffA.*DiffA + DiffB.*DiffB);
subplot(2,2,3);
hist(Diff); 


LightRGB = LightRGB.*(LightRGB>1e-10);

FileAfter = '_bestLight.txt';
T = [T, FileAfter];
format long g
dlmwrite(T,gamma, 'delimiter', ' ', 'newline', 'pc');
format long g
dlmwrite(T,LightRGB, '-append', 'delimiter', ' ');
%%

function [ColorALAB, ColorBLAB]= ColorRGB2LAB(ColorRRGB, ColorGRGB, ColorBRGB)
    lightNum = size(ColorRRGB, 2);

    rgb = zeros(lightNum, 3);

    rgb(:,1) = ColorRRGB;
    rgb(:,2) = ColorGRGB;
    rgb(:,3) = ColorBRGB;
    lab = rgb2lab(rgb);
    
    ColorALAB = lab(:,2)';
    ColorBLAB = lab(:,3)';
%%
function [ColorLRR,ColorLRG, ColorLRB, translate, gamma] = ...
        uniform(ColorLRR,ColorLRG, ColorLRB, maxA, minA)
    maxB = max(max(max(ColorLRR),max(ColorLRG)),max(ColorLRB));
%    minB = min(min(min(ColorLRR),min(ColorLRG)),min(ColorLRB));
%    gamma = (maxA - minA)/(maxB - minB);
    gamma = maxA/maxB;
%   translate = minA/gamma-minB;
   translate = 0;
    ColorLRR = ColorLRR.*gamma;
    ColorLRG = ColorLRG.*gamma;
    ColorLRB = ColorLRB.*gamma;
%    ColorLRB = ColorLRR.*ColorLRB;

%%
function[LightRP, LightGP, LightBP] = Project(LightR, LightG, LightB)
    global Scale;

    LightRGB = [LightR; LightG; LightB];
    
    CenterV = Scale(:,4);

    [m,n] = size(LightRGB);

    CenterT = LightRGB(1,:)*CenterV(1) + ...
        LightRGB(2,:)*CenterV(2) + ...
        LightRGB(3,:)*CenterV(3);

    CenterP(1,:) = CenterT * CenterV(1);
    CenterP(2,:) = CenterT * CenterV(2);
    CenterP(3,:) = CenterT * CenterV(3);

    LightCenterP = LightRGB - CenterP;

    ProjectT = zeros(1,n);
    ProjectT(:) = 0;

    for i=5:7

        UpT = LightRGB(1,:)*Scale(1,i) + ...
            LightRGB(2,:)*Scale(2,i) + ...
            LightRGB(3,:)*Scale(3,i);

        AllT = LightCenterP(1,:)*Scale(1,i) + ...
            LightCenterP(2,:)*Scale(2,i) + ...
            LightCenterP(3,:)*Scale(3,i);

        TempT = (UpT > 0).*UpT./AllT;     
       
        ProjectT  = max(ProjectT, TempT);
    end;
    
    LightRP = LightRGB(1,:) - (ProjectT > 0).*ProjectT.* LightCenterP(1,:);
    LightGP = LightRGB(2,:) - (ProjectT > 0).*ProjectT.* LightCenterP(2,:);
    LightBP = LightRGB(3,:) - (ProjectT > 0).*ProjectT.* LightCenterP(3,:);
    
    
    
 
%%
function [ColorALAB, ColorBLAB]= ColorHSV2LAB(ColorAHSV, ColorBHSV)
    lightNum = size(ColorAHSV, 2);

    hsv = zeros(lightNum, 3);

    hsv(:,1) = ColorAHSV*0.5+0.75;
    hsv(:,1) = hsv(:,1) - (hsv(:,1) > 1);
    hsv(:,2) = ColorBHSV;
    hsv(:,3) = 1;
    rgb = hsv2rgb(hsv);
    lab = rgb2lab(rgb);

    ColorALAB = lab(:,2)';
    ColorBLAB = lab(:,3)';
%%
function [ColorRRGB, ColorGRGB, ColorBRGB]= ColorHSV2RGB(ColorAHSV, ColorBHSV)
    lightNum = size(ColorAHSV, 2);

    hsv = zeros(lightNum, 3);

    hsv(:,1) = ColorAHSV*0.5+0.75;
    hsv(:,1) = hsv(:,1) - (hsv(:,1) > 1);
    hsv(:,2) = ColorBHSV;
    hsv(:,3) = 1;
    rgb = hsv2rgb(hsv);

    ColorRRGB = rgb(:,1)';
    ColorGRGB = rgb(:,2)';
    ColorBRGB = rgb(:,3)';
%%