OnsetPath=[pwd '/optseq/'];
TaskName='Resp';
load([TaskName 'Tables']);
OnsetStruct=load([TaskName 'Tables']);
FileList=fieldnames(OnsetStruct);

for file =1: size(FileList,1)

   Name=FileList{file};
   tbl=eval(Name);
   tbl.Properties.VariableNames = {'onset' 'type' 'duration' 'repeat' 'cond'};
   
   tbl.cond=cell2num(tbl.cond);
   tbl(tbl.cond=='NULL')=[] ;
    
    
end