% Copyright (C) 2013 by The HDF Group.
% All rights reserved.
%
%  This example code illustrates how to access and visualize ICESat/GLAS
% L2 HDF5 file in MATLAB.
%
%  If you have any questions, suggestions, comments on this example, please
% use the HDF-EOS Forum (http://hdfeos.org/forums).
%
%  If you would like to see an  example of any other NASA HDF/HDF-EOS data
% product that is not listed in the HDF-EOS Comprehensive Examples page
% (http://hdfeos.org/zoo), feel free to contact us at eoshelp@hdfgroup.org or
% post it at the HDF-EOS Forum (http://hdfeos.org/forums).

% Tested under: MATLAB R2012a
% Last updated: 2013-1-14

clear;clc;
% Open the HDF5 File.
input_root_path='K:\icesat2\';
res_path="M:\Altimeter\ICESat-2\Icesat2_TP_L3ATL6_slope_geoid.txt";
out_root_path='M:\Altimeter\';
dirOutput1=dir(fullfile(input_root_path,'*.h5'));
fileNames1={dirOutput1.name}';
j=1;
beam={'gt1l','gt1r','gt2l','gt2r','gt3l','gt3r'};
for jj=1:1:6
    for i=1:length(fileNames1)
        Result={};
        try
            FILE_NAME =strcat( input_root_path,fileNames1{i});
            tmp=fileNames1{i};
            file_id = H5F.open (FILE_NAME, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
            
            % Open the datasets.
            LATFIELD_NAME=strcat('',beam{jj},'/land_ice_segments/latitude');
            lat_id=H5D.open(file_id,LATFIELD_NAME);
            
            LONFIELD_NAME=strcat('',beam{jj},'/land_ice_segments/longitude');
            lon_id=H5D.open(file_id, LONFIELD_NAME);
            
            LEVFIELD_NAME=strcat('',beam{jj},'/land_ice_segments/h_li');
            ele_id=H5D.open(file_id, LEVFIELD_NAME);
            
            LEVFIELD_sigma_NAME=strcat('',beam{jj},'/land_ice_segments/h_li_sigma');
            ele_sigma_id=H5D.open(file_id, LEVFIELD_sigma_NAME);
            %
            qf_NAME=strcat('',beam{jj},'/land_ice_segments/atl06_quality_summary');
            qf_ice_id=H5D.open(file_id, qf_NAME);
            
            %             slp_dx_NAME=strcat('',beam{jj},'/land_ice_segments/fit_statistics/dh_fit_dx');
            %             slpdx_ice_id=H5D.open(file_id, slp_dx_NAME);
            %
            %             slp_dx_sigma_NAME=strcat('',beam{jj},'/land_ice_segments/fit_statistics/dh_fit_dx_sigma');
            %             slpdx_ice_sigma_id=H5D.open(file_id, slp_dx_sigma_NAME);
            
            LEVFIELD_geoid_NAME=strcat('',beam{jj},'/land_ice_segments/dem/geoid_h');
            geoid_id=H5D.open(file_id, LEVFIELD_geoid_NAME);
            
            
            lat=H5D.read(lat_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL',...
                'H5P_DEFAULT');
            lon=H5D.read(lon_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', ...
                'H5P_DEFAULT');
            ele=H5D.read(ele_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', ...
                'H5P_DEFAULT');
            ele_sigma=H5D.read(ele_sigma_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', ...
                'H5P_DEFAULT');
            qf_ice=H5D.read(qf_ice_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', ...
                'H5P_DEFAULT');
            %             slp_ice=H5D.read(slpdx_ice_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', ...
            %                 'H5P_DEFAULT');
            geoid=H5D.read(geoid_id,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', ...
                'H5P_DEFAULT');
            
            
            
            % Subset data because there are many fill values in lat/lon.
            index=find(lat>25 & lat<40  & lon >66 & lon<106&qf_ice==0 );
            lat1=lat(index);
            lon1=lon(index);
            ele1=ele(index);
            ele1_sigma=ele_sigma(index);
            %             slp_ice1=slp_ice(index);
            %             slp_sigma_ice1=slp_sigma_ice(index);
            geoid1=geoid(index);
            
            start1=length(Result)+1;
            end1=length(Result)+length(index);
            Result(start1:end1,1)=num2cell(lat1);
            Result(start1:end1,2)=num2cell(lon1);
            Result(start1:end1,3)=num2cell(ele1);
            Result(start1:end1,4)=num2cell(ele1_sigma);
            time=str2double(tmp(7:14))*ones(length(ele1),1);
            Result(start1:end1,5)=num2cell(time);
            %             Result(start1:end1,6)=num2cell(slp_ice1);
            %             Result(start1:end1,7)=num2cell(slp_sigma_ice1);
            Result(start1:end1,6)=num2cell(geoid1);
            
            %           Result(start1:end1,5)=num2cell(wave1);
            % Close and release resources.
            H5D.close(geoid_id)
            %            H5D.close(slpdx_ice_sigma_id)
            %             H5D.close(slpdx_ice_id);
            %             H5D.close(qf_ice_id);
            H5D.close(ele_sigma_id);
            H5D.close(ele_id);
            H5D.close(lon_id);
            H5D.close(lat_id);
            H5F.close(file_id);
            
            dlmwrite(res_path, Result,'-append','precision','%.8f');
        catch
            j;
        end
    end
    
end

