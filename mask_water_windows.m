clear;
close all;

% description: mask water and calculate missing data on land within predefined subsets.
% input:  the image subsets of missing data around city and land/water masks.
% output: masked image subsets and summary statistics for each subset and date.


%-----------------------------------------------------------------------
%set inputs / paths
wd = ('C:\sage\modis\data\cloud_cover\city_windows\');   %working directory
landPath = ('C:\sage\modis\data\cloud_cover\land\');    %directory with water/land mask
S = load('window_size.txt');                                             %text file with sample and line count
outText = fopen('missing_data_summary.txt', 'w');

%get array of city directories
allContents = dir(wd);                  			%get all files/dirs under city_windows directory
allDirs = {allContents([allContents.isdir]).name};  %get the names of all directories

%get rid of parent directories
validIndex = ~ismember(allDirs,{'.','..'});          %get valid directory indices
cityDirs = allDirs(validIndex);                      %get array of city dirs

%loop through city directories
for i=1:length(cityDirs)
   
    city = char(cityDirs(i));                           %get dir name as string
    path = strcat(wd, '\', city, '\*.bsq');             %set search path
        
    mdw = dir(path);                                    %get array of bsq files in city dir
    mdWind = {mdw.name};                                %get mdw as cell
    
    path2 = strcat(landPath, city, '*.bsq');            %get path to water mask directory
    watmsk = dir(path2);                                %get water mask filename
    water = {watmsk.name};                              %get watmsk as cell
    
    s = S(i, 1);                                        %get number of samples for the city window/water mask
    l = S(i, 2);                                        %get number of lines for the city window/water mask
    
    wPath = char(strcat(landPath, water));              %set path to water mask
    fid = fopen(wPath, 'r');                            %open water mask
    W = fread(fid, [s,l]);                              %read water mask with size samples x lines
    fclose(fid);                                        %close image file
    
    
    for j=1:length(mdWind)                              %loop through city windows
        
        mdPath = char(strcat(wd, '\', city, '\', mdWind(j)));        %set path to city window
        mdID = fopen(mdPath, 'r');                                   %open missing data image   
        M = fread(mdID, [s,l]);                                      %read missing data image
        fclose(mdID);                                                %close image file
        
        
        outID = fopen(char(strcat('masked.', mdWind(j))), 'w');      %open output file for writing image
        N = zeros(s, l);                                             %create empty array to write values
        
        for a=1:s                       
            for b=1:l                    %loop through each pixel
                
                if W(a,b) == 0           %if water pixel,     
                    N(a,b) = 0;          %assign NoData
                else
                    N(a,b) = M(a,b);     %otherwise, keep the value from the missing data image
                end
            end    
        end
        
        
        av = round(mean(N(:))*100);
        stdev = round(std(N(:))*100);
        
        fprintf(outText, '%d ', av);
        fprintf(outText, '%d', stdev);
        fprintf(outText, '\n');
        
        
        fwrite(outID, N);                %write N out to binary file         
        fclose(outID);                   %close file
    end
    
    
end


fclose(outText);                        %close text file 






 
 
