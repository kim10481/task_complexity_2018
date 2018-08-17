function [out]=job_del(id,opt)
% (id,opt), opt=0:delete dir&file, opt=1:remove dir only. opt=2: remove
% files only.

name0='/home/swlee/';

for jj=1:1:length(id)
    name1=sprintf('Job%d',id(jj));
    if(opt==1)
        rmdir([name0, name1],'s');
    end
    if(opt==0)
        delete([name0, name1, '*.*']);
        rmdir([name0, name1],'s');
    end
    if(opt==2)
        delete([name0, name1, '*.*']);
    end
    
      
    disp(['- All files associated with ' name1 ' have been deleted.']);
    
end

out=['-done.'];
    
end