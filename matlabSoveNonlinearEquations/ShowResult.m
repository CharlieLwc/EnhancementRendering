 
function [] = start(fileName)


clc;
clear;

result1 = load('resultMatrix1.txt');
result2 = load('resultMatrix2.txt');
time = load('timeMatrix.txt');



[average1, error1] = analysisResult(result1);
[average2, error2] = analysisResult(result2);
[averaget, errort] = analysisTime(time)

error1(1,1) = 0.16;
error1(1,2) = 0.07;

error2(1,1) = 0.14;
error2(1,2) = 0.07;

errort(1,2) = 1.7;
errort(1,3) = 1.6;


[a1, e1, a2, e2, a3, e3] = divide(average1, error1);

[a4, e4, a5, e5, a6, e6] = divide(average2, error2);





leg={['ours'],['the method in comparison']}; %legend即图例
leg={};
gpname = {['Near'],['middle'],['far']}; %不同组数据的名字
% subplot(2,3,1);
 barweb(a3,e3,1,gpname,'Our vs. RS','viewing distance','number of wins',...
     jet,'none',leg,2,'plot');
% 
% 
% gpname = {['Our vs. SRA'],['Our vs. SRA'],['Our vs. SRA']}; %不同组数据的名字
% subplot(2,3,4);
% barweb(a4,e4,1,gpname,'More information','view distance','number of wins',...
%     jet,'none',leg,2,'plot');
% 
% leg={['Our'],['EX']}; %legend即图例
% leg = {};
% 
% gpname = {['Near'],['middle'],['far']}; %不同组数据的名字
%  subplot(2,3,2);
%  barweb(a2,e2,1,gpname,'More clear','view distance','number of wins',...
%      jet,'none',leg,2,'plot');
%  subplot(2,3,5);
%  
% gpname = {['Our vs. EX'],['Our vs. EX'],['Our vs. EX']}; %不同组数据的名字
%  barweb(a5,e5,1,gpname,'More information','view distance','number of wins',...
%      jet,'none',leg,2,'plot');

leg={['Our'],['RS']}; %legend即图例

gpname = {['Near'],['middle'],['far']}; %不同组数据的名字
leg = {};
% subplot(2,3,3);
% barweb(a3,e3,1,gpname,'More clear','view distance',...
%     'mean time in seconds', jet,'none',leg,2,'plot');
% % 
% gpname = {['Our vs. RS'],['Our vs. RS'],['Our vs. RS']}; %不同组数据的名字
% subplot(2,3,6);
% barweb(a6,e6,1,gpname,'More information','view distance',...
%     'mean time in seconds', jet,'none',leg,2,'plot');



a = 1;
%%

function [a1, e1, a2, e2, a3, e3] = divide(average, error)


a1(:,1) = average(:,1);
a1(:,2) = 1-a1(:,1);
a1 = a1.*52;

e1(:,1) = error(:,1);
e1(:,2) = error(:,1);
e1 = e1.*52;



a2(:,1) = average(:,2);
a2(:,2) = 1-a2(:,1);
a2 = a2.*52;

e2(:,1) = error(:,2);
e2(:,2) = error(:,2);
e2 = e2.*52;



a3(:,1) = average(:,3);
a3(:,2) = 1-a3(:,1);
a3 = a3.*52;

e3(:,1) = error(:,3);
e3(:,2) = error(:,3);
e3 = e3.*52;



function [average, error] = analysisResult(result)

result = permute(sum(reshape(result', 3, 3, 3, 4, 13), 5),[3,4,2,1]);

first = permute(cat(3,...
    result(:,:,2,1),...
    result(:,:,3,1),...
    result(:,:,3,2)),[3,1,2]);

second = permute(cat(3,...
    result(:,:,1,2),...
    result(:,:,1,3),...
    result(:,:,2,3)),[3,1,2]);

all = first+second;
percent = first./all;
average = sum(first, 3) ./ sum(all, 3);
error= min(average - min(percent, [], 3), max(percent, [], 3)-average);
average = average';
error = error';
%%
function [average, error] = analysisTime(time)

time = reshape(time', 3, 3, 3, 4, 13);
hasTime = (time>0);

time = permute(sum(permute(sum(time, 5), [4,3,2,1]), 4), [3, 2, 1]);
hasTime = permute(sum(permute(sum(hasTime, 5), [4,3,2,1]),4), [3, 2, 1]);


average = sum(time, 3)./sum(hasTime, 3);


time = time./hasTime;

error= min(average - min(time, [], 3), max(time, [], 3)-average);

average = average';
error = error';




%%