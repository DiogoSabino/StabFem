function [] = mydisp(verbositylevel,string)
global ff ffdir ffdatadir sfdir verbosity

if(verbosity>=verbositylevel) 
    disp(string) 
end
