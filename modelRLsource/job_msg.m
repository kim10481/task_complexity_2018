function [out]=job_msg(id)

name0='/home/swlee/';
name1=sprintf('Job%d/Task1.out.mat',id);

eval(['load ' name0 name1]);

if(length(errormessage)~=0)
out=['- [Message] ' errormessage];
else
    out='- [No message found]';
end

end
