function [out]=KeycodeChecker()

total_trial=10;

for i=1:1:total_trial
    [secs, keyCode] = KbPressWait;
    id_key=find(keyCode==1);
    disp(sprintf('- key code pressed : %d',id_key));
end


out=1;
end


