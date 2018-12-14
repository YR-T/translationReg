function [stdData,displacement]=roundtrip_scan_compensate(stdData,ch)
% ch -> channel for correction
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
if nargin==1
    ch=1;
end
for i=[ch,setdiff(1:length(stdData.Image),ch)]
    if i==ch
        image_4d=stdData.Image{ch};
        [stdData.Image{ch},displacement]=roundtrip_scan_comp_z_data(image_4d);
        stdData.Analysis.roundtrip_displacement=displacement;
    else
         image_4d=stdData.Image{i};
         [stdData.Image{i}]=roundtrip_scan_comp_z_data(image_4d,displacement);
    end
end

function [image_4d_reg,displacement]=roundtrip_scan_comp_z_data(image_4d,varargin)
image_4d_reg=image_4d;
if nargin == 1
    ratio=0;
    prev_sumdisp=inf;
    j=0;
    while (ratio<0.95 && j<11 && prev_sumdisp>50)
         % halt at small change
         % 1/10 pixels > compensation
         % cycles>10
        j=j+1;
        [image_4d_reg(:,:,1,:),displacement(j,:)] = roundtrip_scan_compensate_routin(squeeze(image_4d_reg(:,:,1,:)));
        sumdisp=sum(abs(displacement(j,:)),2);
        ratio=sumdisp/prev_sumdisp;
        prev_sumdisp=sumdisp;
    end
    displacement=sum(displacement);
    for i = 1:size(image_4d,3)
        [image_4d_reg(:,:,i,:)] = roundtrip_scan_compensate_routin(squeeze(image_4d(:,:,i,:)),displacement);
    end
else
    displacement = varargin{1};
    for i = 1:size(image_4d,3)
        [image_4d_reg(:,:,i,:)] = roundtrip_scan_compensate_routin(squeeze(image_4d(:,:,i,:)),displacement);
    end
end

function [series_out,displacement]=roundtrip_scan_compensate_routin(series,varargin)

if nargin>3
    disp('error:too many input!')
    return
end

sizeY=size(series,1);

if nargin==1
    average=mean(series,3);
    
    average_1=average(1:2:sizeY,:);
    average_2=average(2:2:sizeY,:);
    
    i_ave_1_t = interp1(1:2:sizeY,average_1,1:sizeY);
    i_ave_2_t = interp1(2:2:sizeY,average_2,1:sizeY);

    strip_size=100;%even number
    win_size=20;%even number
    n_div=ceil((512-strip_size)/win_size);
    %n_div=47;
    for k=1:n_div
        [~,~,dif_x(k)]=poc_reg(i_ave_2_t(2:sizeY-1,(k-1)*win_size+1:(k-1)*win_size+strip_size),1,0,i_ave_1_t(2:sizeY-1,(k-1)*win_size+1:(k-1)*win_size+strip_size));

    end
     p=polyfit(strip_size/2+1:win_size:(512-strip_size/2-1),dif_x,4);
     displacement=polyval(p,1:512);
      %figure;plot(dif_x),hold all;plot(displacement(26:10:486));pause
elseif nargin==2
    displacement=varargin{1};
end
series_out=uint16(series);
for i=3:510
    try
    if displacement(i) >= 0
        int_part=floor(displacement(i)/2);
        temp11=double(series(1:2:sizeY,i+(int_part+1),:));
        temp12=double(series(1:2:sizeY,i+int_part,:));
        temp21=double(series(2:2:sizeY,i-int_part,:));
        temp22=double(series(2:2:sizeY,i-(int_part+1),:));
        ratio=displacement(i)/2-int_part;
    else
        int_part=floor(abs(displacement(i)/2));
        temp11=double(series(1:2:sizeY,i-(int_part+1),:));
        temp12=double(series(1:2:sizeY,i-int_part,:));
        temp21=double(series(2:2:sizeY,i+int_part,:));
        temp22=double(series(2:2:sizeY,i+int_part+1,:));
        ratio=abs(displacement(i)/2)-int_part;
    end
    series_out(1:2:sizeY,i,:)=uint16(temp11*ratio+temp12*(1-ratio));
    series_out(2:2:sizeY,i,:)=uint16(temp21*(1-ratio)+temp22*ratio);
    catch
    end
end
