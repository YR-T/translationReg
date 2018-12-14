function [regged_image_series] = translation_reg_log(image_series,dif_y,dif_x)
% Copyright:
% 2014, 2018, Yasuhiro R. Tanaka; 
% License:
% This code is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 2 of the License,
% or any later version. This work is distributed in the hope that it 
% will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See
% version 2 and version 3 of the GNU General Public License for more
% details. You should have received a copy of the GNU General Public
% License along with this program;If not, see http://www.gnu.org/licenses/.
regged_image_series = zeros(size(image_series),'uint16');
[w_row,w_col,num_frame] = size(image_series);

for i=1:num_frame
    image_temp=double(image_series(:,:,i));
    row_disp = dif_y(i);
    col_disp = dif_x(i);
    
    %row direction
    if (row_disp >=0)
        reg_image1=[zeros(floor(row_disp),w_col); image_temp; zeros(101,w_col)];
        reg_image2=[zeros(ceil(row_disp),w_col); image_temp; zeros(100,w_col)];
        ratio = row_disp-floor(row_disp);
    else
        reg_image1 = [image_temp(floor(abs(row_disp))+1:w_row,:); zeros(100,w_col)];
        reg_image2 = [image_temp(ceil(abs(row_disp))+1:w_row,:); zeros(101,w_col)];
        ratio = abs(row_disp)-floor(abs(row_disp));
    end
    reg_image = (1-ratio)*reg_image1 + ratio*reg_image2;
    
    %column direction
    [temp_w_row, temp_w_col]=size(reg_image);
    if (col_disp>=0)
        reg_image1=[zeros(temp_w_row,floor(col_disp)), reg_image, zeros(temp_w_row,101)];
        reg_image2=[zeros(temp_w_row,ceil(col_disp)), reg_image, zeros(temp_w_row,100)];
        ratio = col_disp-floor(col_disp);
    else
        reg_image1 = [reg_image(:,floor(abs(col_disp))+1:temp_w_col), zeros(temp_w_row,100)];
        reg_image2 = [reg_image(:,ceil(abs(col_disp))+1:temp_w_col), zeros(temp_w_row,101)];
        ratio = abs(col_disp)-floor(abs(col_disp));
    end
    reg_image = (1-ratio)*reg_image1 + ratio*reg_image2;
    
    regged_image_series(:,:,i) = uint16(reg_image(1:w_row,1:w_col));
end
