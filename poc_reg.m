function [regged_image_series,dif_y,dif_x]=poc_reg(image_series,mag,cut,target)

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

debug=0;debug_itr=0;
subset = image_series(cut+1:end-cut,cut+1:end-cut,:);

if nargin==4
    mean_image=target;
else
    mean_image = mean(subset,3);
end
[w_row,w_col,n_im]=size(subset);
[wx,wy]=meshgrid(0.5-0.5*cos(2*pi/size(mean_image,2)*linspace(0,size(mean_image,2),size(mean_image,2))),0.5-0.5*cos(2*pi/size(mean_image,1)*linspace(0,size(mean_image,1),size(mean_image,1))));%2dim hanning window
mean_image=mean_image.*wx.*wy;
ft_mean = fft2(mean_image);

if debug==1
    figure;imagesc(mean_image);title('template')
    figure;imagesc(ifft2(ft_mean./abs(ft_mean)));title('phase')
end

parfor i = 1:n_im
    if ismember(i,debug_itr)==1
        debug=1;
    else
        debug=0;
    end
    temp=double(subset(:,:,i));
    temp=temp.*wx.*wy;
    ft_temp = fft2(temp);
    
    factorial=ft_temp.*conj(ft_mean)./sqrt(abs(ft_mean).*abs(ft_temp));
    if mag~=1
        zeros_ft=zeros(size(factorial));
        factorial = [factorial(1:floor(w_row/2),1:floor(w_col/2)),repmat(zeros(floor(w_row/2),w_col),1,(mag-1)),factorial(1:floor(w_row/2),floor(w_col/2)+1:w_col);
            repmat(zeros_ft,(mag-1),mag);
            factorial(floor(w_row/2)+1:w_row,1:floor(w_col/2)),repmat(zeros(w_row-floor(w_row/2),w_col),1,(mag-1)),factorial(floor(w_row/2)+1:w_row,floor(w_col/2)+1:w_col)];
    end
    c = fftshift(real(ifft2(factorial)));
    [c1,i1] = max(c);
    [~,i2] = max(c1);
    if debug==1
        figure;imagesc(temp);title(['image ',num2str(i)])
        figure;imagesc(real(ifft2(ft_temp./abs(ft_temp))));title('phase')
        figure;imagesc(c);title('correlation');
        hold all;scatter(floor(w_col/2*mag)+1,floor(w_row/2*mag)+1,'.k');
        scatter(i2,i1(i2),'.r');pause
    end
    
 
    val_y=(c(i1(i2)-1,i2)-c(i1(i2)+1,i2))/2/(c(i1(i2)-1,i2)+c(i1(i2)+1,i2)-2*c(i1(i2),i2))+i1(i2);
    val_x=(c(i1(i2),i2-1)-c(i1(i2),i2+1))/2/(c(i1(i2),i2-1)+c(i1(i2),i2+1)-2*c(i1(i2),i2))+i2;
    
    dif_y(i)=-(val_y-(floor(w_row/2*mag)+1))/mag;
    dif_x(i)=-(val_x-(floor(w_col/2*mag)+1))/mag;
    
    if abs(dif_y(i))>100
        dif_y(i)=eps;
        disp(['A displacement > 100 pix was detected in frame ',num2str(i),' and was assumed as eps'])
    end
    if abs(dif_x(i))>100
        dif_x(i)=eps;
        disp(['A displacement > 100 pix was detected in frame ',num2str(i),' and was assumed as eps'])
    end
    if debug==1
        dif_x(i),dif_y(i),pause
    end
end
if debug == 1
    figure;plot(dif_y);hold all;plot(dif_x);
end
[regged_image_series] = translation_reg_log(image_series,dif_y,dif_x);
regged_image_series=uint16(regged_image_series);

end