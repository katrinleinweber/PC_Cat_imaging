
clear all

FileList=dir('*_onset_length*');
for file=1:length(FileList)
    clear onsetlist Dvec
    FileName{file}=FileList(file).name;
    
    load(FileName{file})
    
    
    
    for col=1:(length(onsetlist)-1)
        Dvec(col)= onsetlist(col+1)-onsetlist(col);
    end

    
    FileList(file).Mean=mean(Dvec);
    FileList(file).Std=std(Dvec);
end