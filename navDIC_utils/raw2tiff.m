

clc

[file,path] = uigetfile('.raw','Select a .RAW file to be converted') ;

readraw;%('clean')

imread([path,file]) ;
