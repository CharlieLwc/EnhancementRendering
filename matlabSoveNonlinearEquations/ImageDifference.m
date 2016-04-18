[filename, pathname] = uigetfile( ...
    {'*.jpg;*.tif;*.png;*.gif','All Image Files';...
    '*.*','All Files' },...
    '??????????????', ...
    'MultiSelect', 'on');

fileDN = fullfile(pathname, filename);
fileNum = size(fileDN,2);

line = 0;
total = fileNum*(fileNum-1)/2;
for i=1:fileNum
    for j=i+1: fileNum
        first = imread(char(fileDN(i)));
        second = imread(char(fileDN(j)));
        
subplot(total, 4, line*4+1);
imshow(first);
        
subplot(total,4, line*4+2);
imshow(second);

dif1 = first - second;
gamma = 255/max(max(max(dif1)));
dif1 = dif1*gamma;
subplot(total,4, line*4+3);
imshow(dif1);

dif = second - first;
gamma = 255/max(max(max(dif)));
dif = dif*gamma;
subplot(total,4, line*4+4);
imshow(dif);
        
line = line+1;
    end
    
end






