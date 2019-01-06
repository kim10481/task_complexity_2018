function [out]=SIMUL_cogload_fmri_init(map_name)

switch lower(map_name)


    case {'sangwan2014b'}

        %% relocating all the seed files (from seed_images folder) according to the case #


        seed_path0=pwd;%['D:\0-program\One shot learning']; % laptop
        seed_path=[seed_path0 '\seed\'];
        seed_path_pool=[seed_path '\seed_images\'];

        % [CAUTION] the union must be the entire set. There must be NO the
        % intersection.

        % outcome color - to be randomized
        index_outcome_randomized=[1 2 3];
        index_outcome_fixed=[]; % must be empty

        % non-outcome state indices  - to be randomized
        num_seed_stateimg=30;
        index_state_randomized=[1 2 3 4 10];
        % outcome state indices - to be randomized according to outcome color
        index_state_outcome=[5 6 7 8 9 11];    index_state_outcome_type=[1 1 2 3 2 3];




        %% relocating outcome files
        index_outcome_relo=index_outcome_randomized(randperm(length(index_outcome_randomized)));
%         index_outcome_relo=[1 3 2]; % test only (the same setting as ppt)

        index_outcome_src=[index_outcome_randomized index_outcome_fixed];
        index_outcome_dst=[index_outcome_relo index_outcome_fixed];

        for i=1:1:length(index_outcome_src)

            % relocating outcome state files (except outcome states)
            src_file=[seed_path_pool sprintf('o%03d.png',index_outcome_src(i))];
            dst_file=[seed_path sprintf('o%03d.png',index_outcome_dst(i))];
            if (exist(dst_file, 'file') == 2)
                delete(dst_file); WaitSecs(0.2);
            end
            [status,message,messageId]=copyfile(src_file,dst_file,'f'); WaitSecs(0.3);
            if(status~=1)
                disp('-copy error. try to delete files in the seed folder.');
            end

            file_name_outcome.src{1,i}=src_file;
            file_name_outcome.dst{1,i}=dst_file;

            disp(sprintf('- outcome file copied (%d/%d)',i,length(index_outcome_src)));
        end

        out.file_name_outcome=file_name_outcome;



        %% relocating state files

        total_n_state=length([index_state_randomized index_state_outcome]);

        % 1. only the ones in "index_state_randomized"

        rand_ind=randperm(num_seed_stateimg);
        index_state_relo=rand_ind([1:1:length(index_state_randomized)]);
        
        index_state_src=[index_state_relo];
        index_state_dst=[index_state_randomized];
        
        for i=1:1:length(index_state_dst)
            % relocating outcome state files (except outcome states)
            src_file=[seed_path_pool sprintf('s%03d.png',index_state_src(i))];
            dst_file=[seed_path sprintf('s%03d.png',index_state_dst(i))];
            if (exist(dst_file, 'file') == 2)
                delete(dst_file); WaitSecs(0.2);
            end
            [status,message,messageId]=copyfile(src_file,dst_file,'f'); WaitSecs(0.3);
            if(status~=1)
                disp('-copy error. try to delete files in the seed folder.');
            end
            file_name_state.src{1,i}=src_file;
            file_name_state.dst{1,i}=dst_file;
            disp(sprintf('- stimulus file copied (%d/%d)',i,total_n_state));
        end


        % 2. ones associated with outcomes
        i0=length(index_state_src);
        for i=1:1:max(index_state_outcome_type) % for each O type

            % select state files associated with O type
            ind_file_move=find(index_state_outcome_type==index_outcome_dst(i));
            % move each state file
            for hh=1:1:length(ind_file_move)
                i0=i0+1;
                src_file=[seed_path_pool sprintf('s_o%03d_%d.png',i,hh)];
                dst_file=[seed_path sprintf('s%03d.png',index_state_outcome(ind_file_move(hh)))];
                if (exist(dst_file, 'file') == 2)
                    delete(dst_file); WaitSecs(0.2);
                end
                [status,message,messageId]=copyfile(src_file,dst_file,'f'); WaitSecs(0.3);
                if(status~=1)
                    disp('-copy error. try to delete files in the seed folder.');
                end
                file_name_state.src{1,i0}=src_file;
                file_name_state.dst{1,i0}=dst_file;
                disp(sprintf('- stimulus file copied (%d/%d)',i0,total_n_state));
            end
        end

        out.file_name_state=file_name_state;



        %%
        disp('initialization completed.');
        
        
    otherwise,
        error('- file initialization ERROR. check the map name!!!');

end

end