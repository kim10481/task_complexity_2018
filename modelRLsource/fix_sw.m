rwd=[40 20 10 0];

for session=1:1:5
    
    mat0=HIST_behavior_info{1,session};
    [row tmp]=find(mat0(:,18)==9);
    [row_rest tmp]=find(mat0(:,18)==9);
    for j=1:1:length(row)
        mat0(row(j),18)=-1;
        mat0(row(j),16)=rwd(mat0(row(j),6)-5);
    end
    
    [r_sz c_sz]=size(mat0);    
    mat_part_fixed=mat0(row,:);     mat_part_fixed_ind=1;
    mat1=[];    check_old=10;
    for h=1:1:length(row_rest)        
        check_new=mat0(row_rest(h),18);
        if((check_old==-1)&&(check_new~=-1))
            if(mat_part_fixed_ind<=length(row))
                mat1=[mat1; mat_part_fixed(mat_part_fixed_ind,:)];
                mat_part_fixed_ind=mat_part_fixed_ind+1
            end
            if(mat_part_fixed_ind<=length(row))
                mat1=[mat1; mat_part_fixed(mat_part_fixed_ind,:)];
                mat_part_fixed_ind=mat_part_fixed_ind+1;
            end
        else
            mat1=[mat1; mat0(row_rest(h),:)];
        end
        check_old=check_new;
    end
    
    sum0=0;
    for k=1:1:r_sz
        sum0=sum0+mat1(k,16);
        mat1(k,17)=sum0;
    end
    
    HIST_behavior_info_fix{1,session}=mat1;
end

