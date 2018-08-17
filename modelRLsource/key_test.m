while 1 % temporarily disabled for test APR 21
    [ keyIsDown, timeSecs, keyCode ] = KbCheck;
    if keyIsDown
        [tmp tmp_key_code]=find(keyCode==1);
        tmp_key_code
        if ((tmp_key_code==53)) % trigger

% KbName(keyCode)

%         if (KbName(keyCode) == 5)
            break;
        end
        % If the user holds down a key, KbCheck will report multiple events.
        % To condense multiple 'keyDown' events into a single event, we wait until all
        % keys have been released.
        while KbCheck; end
    end
end
