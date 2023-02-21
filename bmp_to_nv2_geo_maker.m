
%% defining geometry bmp
input = imread('SAN_2D_with_hct_centre_shifted.bmp');
geo = zeros(size(input,1), size(input,2), 4);
fib_x = zeros(size(geo));
fib_y = zeros(size(geo));
fib_z = zeros(size(geo));
for i = 1:size(input,1)
    for j = 1:size(input,2)
    if sum(input(i,j,:),'all') == 0 % void
        geo(i,j,3) = 0;
    elseif sum(input(i,j,:), 'all') == 3*255 % wall
        geo(i,j,3) = 100;
    elseif (input(i,j,1) == 255) && (input(i,j,2) == 0) && (input(i,j,3) == 0) % red, atria
        geo(i,j,3) = 1;
    elseif (input(i,j,1) == 0) && (input(i,j,2) == 255) && (input(i,j,3) == 0) % green, pathway
        geo(i,j,3) = 24;
    elseif (input(i,j,1) == 0) && (input(i,j,2) == 0) && (input(i,j,3) > 0) % blue, SAN
        if input(i,j,3) == 155 % SAN head
            geo(i,j,3) = 21; 
        elseif input(i,j,3) == 205 % SAN central
            geo(i,j,3) = 22;
        elseif input(i,j,3) == 255 % SAN tail
            geo(i,j,3) = 23;
        end
    end
    fib_x(i, j, 3) = 0;
    fib_y(i, j, 3) = 0;
    fib_z(i, j, 3) = 0;
    end
end

%% reading in fibres bmp
input = double(imread('hist_fibres_old.bmp'));
fib_x = zeros(size(geo));
fib_y = zeros(size(geo));
fib_z = zeros(size(geo));
for i = 1:size(input,1)
    for j = 1:size(input,2)
        if ~(sum(input(i,j,:), 'all') == 3*255)
            fib_x(i, j, 3) = (input(i,j,1)-127)/127;
            fib_y(i, j, 3) = (input(i,j,2)-127)/127;
            fib_z(i, j, 3) = (input(i,j,3)-127)/127;
        end
    end
end
%% zero-padding
% temp = zeros(size(geo)+3);
% temp(3:end-1, 3:end-1, 3:end-1) = geo;
% geo = temp;
% temp(3:end-1, 3:end-1, 3:end-1) = fib_x;
% fib_x = temp;
% temp(3:end-1, 3:end-1, 3:end-1) = fib_y;
% fib_y = temp;
% temp(3:end-1, 3:end-1, 3:end-1) = fib_z;
% fib_z = temp;
% clear temp

%%
fib_norms = (fib_x.^2 + fib_y.^2 + fib_z.^2).^0.5;
fib_x_n = fib_x./fib_norms;
fib_y_n = fib_y./fib_norms;
fib_z_n = fib_z./fib_norms;
fib_x_n(isnan(fib_x_n)) = 0;
fib_y_n(isnan(fib_y_n)) = 0;
fib_z_n(isnan(fib_z_n)) = 0;
fib_norms_2 = (fib_x_n.^2 + fib_y_n.^2 + fib_z_n.^2).^0.5;
% low-memory output can only output integers between 0 and 127. When the voltage in
% the simulation is greater than or equal to V_scale_high the output is
% 127, when it is less than or equal to V_scale_low the output is 0. The
% voltages in between get scaled evenly to output integers 1 to 126.
V_scale_high = 25;
V_scale_low = -85;
% 
%Remove unconnected regions in data
% temp = logical(geo);
% labels = bwlabeln(temp, 6);
% stats = regionprops(labels, 'Area');
% regions = [stats.Area];
% [~,biggest] = max(regions);
% data = geo;
% data(labels ~= biggest) = 0;
% geo = data;
% clear temp data regions labels stats biggest

%% supplimentary variables
id = geo;
id(id == 0) = -1;
index = 0;
for i = 1:size(geo,1)
    for j = 1:size(geo,2)
        for k = 1:size(geo,3)
            if(geo(i,j,k) > 0)
                id(i,j,k) = index;
                index = index + 1;
            end
        end
    end
end
num_cells = index;



%% create non void geometry as a matlab variable
tic
non_void_geo_file_array = nan(num_cells+2, 1+3+1+3+18);
                               % num_cells, X dimension, Y dimension, Z
non_void_geo_file_array(1,:) = [num_cells, size(geo,1), size(geo,2), size(geo,3), length(unique(geo(1:end-1))), V_scale_high, V_scale_low, (~isempty(nonzeros(cat(2,fib_x,fib_y,fib_z)))), zeros(1,18)];
for i = 1:size(geo, 1)
    for j = 1:size(geo,2)
        for k = 1:size(geo,3)
            if (geo(i,j,k) > 0)
                fibre_vector = [fib_x_n(i,j,k), fib_y_n(i,j,k), fib_z_n(i,j,k)];
                if norm(fibre_vector) ~= 0
                    fibre_vector_normalised = fibre_vector/norm(fibre_vector);
                    if((i>1)&&(i<size(geo,1))&&(j>1)&&(j<size(geo,2))&&(k>1)&&(k<size(geo,3))) % if any connected neghbours of a cell are non-cells, then the cell gets its fibres removed
                        %if(nnz(geo(i-1:i+1, j, k)) < 3 || nnz(geo(i, j-1:j+1, k)) < 3 || nnz(geo(i, j, k-1:k+1)) < 3) % considers 6 neighbours (faces)
                        if(nnz(geo(i-1:i+1, j, k)) < 3 || nnz(geo(i, j-1:j+1, k)) < 3) % considers 4 neighbours (sides) (2D)
                        
                            fibre_vector_normalised = [0, 0, 0];
                        end
                    end
                else
                    fibre_vector_normalised = [0, 0, 0];
                end
                non_void_geo_file_array(id(i,j,k)+1 + 1, :) = [id(i,j,k), i-1, j-1, k-1, geo(i,j,k),...
                    fibre_vector_normalised(1), fibre_vector_normalised(2), fibre_vector_normalised(3),...
                    id(i-1,j-1,k), id(i-1,j,k-1), id(i-1,j,k),   id(i-1,j,k+1), id(i-1,j+1,k),...
                    id(i,j-1,k-1), id(i,j-1,k),   id(i,j-1,k+1), id(i,j,k-1),   id(i,j,k+1), id(i,j+1,k-1), id(i,j+1,k), id(i,j+1,k+1),...
                    id(i+1,j-1,k), id(i+1,j,k-1), id(i+1,j,k),   id(i+1,j,k+1), id(i+1,j+1,k)];
            end
        end
    end
end
geo_name = 'SAN_2D_with_hct_centre_shifted_20p_fib_seed_large_7_small_42';
non_void_geo_file_array(end,:) = 10000+non_void_geo_file_array(1,:);
save([geo_name,'.mat'], 'non_void_geo_file_array')
toc

% create .dat file
tic
fid = fopen([geo_name,'.dat'], 'w');
for i = 1:non_void_geo_file_array(1,1)+2
    for j = 1:size(non_void_geo_file_array,2)-1
        if (j == 6) || (j == 7) || (j == 8) % print float
            fprintf(fid, "%f ", non_void_geo_file_array(i,j));
        else % print int
            fprintf(fid, "%d ", non_void_geo_file_array(i,j));
        end
    end
    fprintf(fid, "%d\n", non_void_geo_file_array(i,end)); 
end
fclose(fid);
toc
