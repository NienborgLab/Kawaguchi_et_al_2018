function mypath = pathfinder
% find the path

mypath = cd;
l = 'Kawaguchi2018_pupil';
c = strfind(mypath, l);
mypath = mypath(1:c+length(l)-1);
addpath(genpath(mypath))