function [] = non_void_reader_V3_auto_saver_function(geo_file_name, dir_to_read_from, dir_to_save_to, saved_mat_file_name, ms_per_file)
    %% first run this section to initialise
    geometry_file_name = geo_file_name;
    dims = dlmread(geometry_file_name);
    if(size(dims,2) == 1)
        dims = transpose(reshape(dims,[length(dims)/26, 26]));
    end
    Nx = dims(1, 2);
    Ny = dims(1, 3);
    Nz = dims(1, 4);
    directory = dir_to_read_from;
    files = dir(strcat(directory,'/*.vtk')); 
    nfiles = length(files);
    
    save([dir_to_save_to,'/',saved_mat_file_name], 'nfiles', '-v7.3');
    save([dir_to_save_to,'/',saved_mat_file_name], 'ms_per_file', '-append');
    %%
    split_E_step = 10;
    save([dir_to_save_to,'/',saved_mat_file_name], 'split_E_step', '-append');
    split_sizes = split_E_step*ones(1, floor(nfiles/split_E_step));
    if mod(nfiles, split_E_step) ~= 0
        split_sizes = [split_sizes, mod(nfiles, split_E_step)];
    end
    save([dir_to_save_to,'/',saved_mat_file_name], 'split_sizes', '-append');
    
    names_of_split_E_vars = cell(length(split_sizes),1);
    E_combining_string = ['E = cat(4'];
    temp_clear_string = ['clear'];
    for i = 1:length(split_sizes)
        names_of_split_E_vars{i} = ['E_', num2str(1+split_E_step*(i-1)),'_to_', num2str(split_E_step*i)];
        E_combining_string = [E_combining_string,', ' , names_of_split_E_vars{i}];
        temp_clear_string = [temp_clear_string, ' ', names_of_split_E_vars{i}];
    end
    E_combining_string = [E_combining_string,');', temp_clear_string, ';'];
    save([dir_to_save_to,'/',saved_mat_file_name], 'E_combining_string', '-append');
    save([dir_to_save_to,'/',saved_mat_file_name], 'names_of_split_E_vars', '-append');
    %disp(nfiles)
    %disp(split_sizes)
    %disp(' ')
    %% if low memory, non human-readable output, use this section
    V_scale_high = dims(1, 6);
    V_scale_low = dims(1, 7);
    V_scale_gradient = 127.0/(V_scale_high - V_scale_low);
    V_scale_constant = -V_scale_low*V_scale_gradient;
    split_section_index = 1;
    E = zeros([Nx, Ny, Nz, split_sizes(split_section_index)], 'int8');
    for i = 1:1:nfiles
        A = zeros([Nx Ny Nz], 'int8') - 100;
        infile = fullfile(strcat(directory,'/'),files(i).name);
        fid = fopen(infile, 'r');
        aa = fscanf(fid, '%c');
        fclose(fid);
        for j = 2:size(dims,1)-1
            A(dims(j,2)+1, dims(j,3)+1, dims(j,4)+1) = (aa(j-1) - V_scale_constant)/V_scale_gradient;
        end
        
        if ((mod(i-1, split_E_step) == 0) && (i ~= 1))
            var_name = names_of_split_E_vars{split_section_index};
            eval([var_name,' = E;']);
            save([dir_to_save_to,'/',saved_mat_file_name], var_name, '-append');
            clear(var_name);
            
            split_section_index = split_section_index + 1;
            E = zeros([Nx, Ny, Nz, split_sizes(split_section_index)], 'int8');
        end
        
        E(:,:,:, mod(i-1, split_E_step)+1) = A;
        if mod(i,10) == 0
            disp(100*i/nfiles)
        end
    end
    
    var_name = names_of_split_E_vars{split_section_index};
    eval([var_name,' = E;']);
    save([dir_to_save_to,'/',saved_mat_file_name], var_name, '-append');
    clear(var_name);
end