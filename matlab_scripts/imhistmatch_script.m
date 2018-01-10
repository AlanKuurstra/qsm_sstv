function [result]=imhistmatch_script(inputLoc,refLoc,OutFilename)   
    addpath ./NIfTI_20140122
    inputObj=load_nii(inputLoc);
    input=double(inputObj.img);
    
    refObj=load_nii(refLoc);
    ref=double(refObj.img);
    
    input=input/max(input(:));
    ref=ref/max(ref(:));

    input_indx=input>0;
    ref_indx=ref>0;

    result=input;
    result(input_indx)=imhistmatch(input(input_indx),ref(ref_indx));
    
    tmp=make_nii(result);
    tmp.hdr=inputObj.hdr;
    save_nii(tmp,OutFilename);
return


