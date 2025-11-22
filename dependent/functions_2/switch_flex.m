function trans = switch_flex(state_path)

    for i=1:size(state_path,1)
        count=0;
        for j=1:size(state_path,2)-1
            if state_path(i,j) ~= state_path(i,j+1)
                count=count+1;
            end
        end
        trans(i,1) = count/size(state_path,2);
    end
end
            
            