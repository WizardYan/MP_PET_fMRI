#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 10:56:30 2024

@author: yanw4
"""



import os
import numpy as np
from nilearn import plotting, image, input_data
from nilearn.input_data import NiftiMasker
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
from nilearn.plotting import view_img

from scipy import ndimage
from scipy.ndimage import label, find_objects
from scipy.stats import pearsonr

from nilearn.plotting import view_img
import nibabel as nib

# Automatically list all .nii and .nii.gz files under /PET folder
PL_pet_folder = '/Users/yanw4/Documents/LNI_P1_MP_DFNC_202212/data/PET/rac_PL_day/logan_ref/'
PL_pet_images = [os.path.join(PL_pet_folder, f) for f in os.listdir(PL_pet_folder) if f.endswith('.nii') or f.endswith('.nii.gz')]
MP_pet_folder = '/Users/yanw4/Documents/LNI_P1_MP_DFNC_202212/data/PET/rac_MP_day/logan_ref/'
MP_pet_images = [os.path.join(MP_pet_folder, f) for f in os.listdir(MP_pet_folder) if f.endswith('.nii') or f.endswith('.nii.gz')]

# Sort the list for consistent ordering
PL_pet_images.sort()
MP_pet_images.sort()


for dynamic_feature_name in ['fractional_occ_fpn_plus','fractional_occ_som_plus','fractional_occ_vis_minus','dwell_time_fpn_plus','dwell_time_vis_minus','random']:
    if dynamic_feature_name == 'fractional_occ_fpn_plus':
        dynamic_features = np.array([0.0306553749904631,0.00417332731332901,0.0412645157336765,-0.0429395028451632,0.0379269188018588,-0.0183973084569506,0.140883323088296,0.0368003591739000,0.00268771983483712,0.0147607589444203,0.0546524686652075,0.0382503334815474,-0.0438413361169102,0.0403294822421909,-0.0223358350697966,-0.0188830720745614,0.00542986425339365,-0.0398667575138163,-0.0298491075069800,0.0141292293324513,0.0761718750000000,0.0732177263969172,0.00183084678792833,0.00483248225183709,0.0306112778392039,-0.0181409610215054,-0.0153609831029186,-0.0389127495426708,-0.0174947620472912,0.0726376049201253,0.0484400474284867,-0.0398727023495555,0.0978091057878292,0.0167497396220801,0.0376062992125984,0.0721144024514811,0.0282400806605638	])
    if dynamic_feature_name == 'fractional_occ_som_plus':
        dynamic_features = np.array([-0.0351186388952468,-0.0767225144093684,-0.0567454116493550,-0.0120919436957173,-0.0275295208173544,-0.00735204159657441,-0.0765565425610986,-0.0123091888656091,-0.0467808533414895,-0.0423099680485977,0.0834173573027076,-0.0889691714836224,-0.0313152400835073,-0.0257809157038939,0.0555191871507865,-0.0422544677863827,-0.0717020536025061,-0.00264213793625559,0.0364526147378988,-0.112520186003074,-0.00976562500000000,-0.0635838150289017,-0.00876251788268956,-0.0617137875202391,-0.0581385247196335,-0.0283014112903226,-0.0749167946748592,-0.0677125586574405,-0.00604609398383715,0.0558155726608131,-0.00154142581888247,0.0463057227829126,-0.201774289540247,-0.117069632495164,-0.102015748031496,-0.0125297923050732,0.187180478011557])
    if dynamic_feature_name == 'fractional_occ_vis_minus':
        dynamic_features = np.array([-0.0660639353017472,0.105683616232860,-0.0560713300854858,0.00247079964061098,0.0206061252363511,0.0181411530815109,-0.0320488869068425,-0.00894941634241245,0.0623298669521334,0.0523135182044101,-0.115780931067555,0.0959278197717504,0.102296450939457,0.00144415917843388,-0.0675468845082889,0.0427921704517449,0.0637080867850099,0.0461579226285109,-0.0362016500925432,0.0635951514679845,-0.0410156250000000,0.0520231213872832,0.0471421758975407,0.0545356416324158,0.0539132838414153,0.0598496303763441,0.0134408602150538,0.0525172989739919,-0.00313528883567794,-0.0916450702046184,0.0386505113383726,0.0245937778062799,0.152964588602887,0.0784444279125130,0.0982519685039370,0.0761548064918851,-0.0787328829489076])
    if dynamic_feature_name == 'dwell_time_fpn_plus':
        dynamic_features = np.array([0.420750000000000,-2.327082353,1.620771429,0.259875,0.522605769,	0,2.7324,2.512218045,-0.12239011,1.03275,2.685375,1.370769231,-1.55925,2.1384,	-0.22275,1.976793522,-1.782,1.62,1.447875,-0.004895604,	0.537428571,1.896557143,1.188,0.729,0.540964286,-0.3601125,-0.616846154,-0.984016484,-0.794475,4.0392,	0.81675,-1.701635294,1.7982,2.074757143,1.998,3.1878,0.0891 ])  # Make sure this is in the same order as pet_images
    if dynamic_feature_name == 'dwell_time_vis_minus':
        dynamic_features = np.array([-1.50727500000000,0.489789473684211,-0.615599999999999,2.27991176470588,0.247500000000000,0.567000000000000,0.507184615384615,-0.675000000000000,2.33030769230769,-0.102807692307692,-0.594000000000000,1.02446963562753,1.31528571428571,0.648000000000001,-0.447784615384616,0.356400000000001,-0.193923529411764,1.47150000000000,-0.148499999999999,-0.0185624999999998,2.55420000000000,2.70669327731092,2.49480000000000,3.68280000000000,1.50975000000000,1.35135000000000,0.763714285714285,-0.445500000000000,-0.438750000000000,-0.811038461538461,0.483685714285714,2.12355000000000,2.05354285714286,2.15325000000000,2.09727692307692,0.618750000000000,0.0989999999999993])
    if dynamic_feature_name == 'random':
        dynamic_features =np.random.normal(size=37)
    
    # Load the atlas and only keep regions labeled 1 to 32
    atlas_img = image.load_img('./data/rSchaefer2018_200Parcels_7Networks_order_Tian_Subcortex_S2_3T_MNI152NLin2009cAsym_2mm.nii')
    atlas_data = atlas_img.get_fdata()
    #mask_data = np.isin(atlas_data, np.array([5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32]))
    mask_data = np.isin(atlas_data, np.array([9, 10, 13, 14, 15, 16, 25, 26, 29, 30, 31, 32]))
    #mask_data = np.isin(atlas_data, np.array(range(1,233)))

    mask_img = image.new_img_like(atlas_img, mask_data)
    
    
    masker = NiftiMasker(mask_img=mask_img, standardize=True)
    
    # Load images
    PL_pet_imgs = image.concat_imgs(PL_pet_images)
    MP_pet_imgs = image.concat_imgs(MP_pet_images)
    
    PL_pet_data_2d = masker.fit_transform(PL_pet_imgs)
    MP_pet_data_2d = masker.fit_transform(MP_pet_imgs)
    
    delta_pet_data_2d = (PL_pet_data_2d-MP_pet_data_2d)/PL_pet_data_2d
    
    correlation_coefficients = np.zeros(delta_pet_data_2d.shape[1])
    p_values = np.zeros(delta_pet_data_2d.shape[1])
    
    
    # Perform Pearson correlation for each voxel
    for i in range(delta_pet_data_2d.shape[1]):
        correlation_coefficients[i], p_values[i] = pearsonr(delta_pet_data_2d[:, i], dynamic_features)
    
    
    # Threshold coefficients by p-value
    significant_correlation = np.where(p_values < 0.005, correlation_coefficients, 0)
    
    # Convert the significant correlation values back into a Nifti image
    correlation_img = masker.inverse_transform(significant_correlation)
    
    
    # Save this Nifti image
    correlation_file_name = './results/example/rac_dopamine_correlation_'+ dynamic_feature_name +'_stat_map.nii.gz'
    correlation_img.to_filename(correlation_file_name)
    
   
   
    plotting.plot_glass_brain(correlation_img, colorbar=True, cmap='RdBu_r', plot_abs=False, symmetric_cbar=True,vmax=0.5, vmin=-0.5)
    
    # view = view_img(correlation_img)
    # view.open_in_browser()
    
    ######################
    #### filter 20 voxels
    ######################
    img = nib.load(correlation_file_name)
    origin_image_data = img.get_fdata()
    image_data = np.abs(origin_image_data)
        
    labeled_array, num_features = ndimage.label(image_data, structure=np.ones((3, 3, 3)))

    # Step 2: Filter clusters by size
    sizes = ndimage.sum(image_data, labeled_array, range(num_features + 1))
    cluster_mask = sizes > 30
    filtered_clusters_mask = cluster_mask[labeled_array]
    # Step 3: Create a mask for the original image
    # This mask has 'True' for voxels in clusters > 30 voxels, and 'False' elsewhere
    mask = filtered_clusters_mask > 0
    # Step 4: Apply the mask to the original image
    # This will set voxels not in the desired clusters to 0, maintaining the original voxel values otherwise
    filtered_image = np.zeros_like(image_data)
    filtered_image[mask] = origin_image_data[mask]
    
    aa_filtered_image = nib.Nifti1Image(filtered_image, affine=img.affine)

    plotting.plot_glass_brain(aa_filtered_image, colorbar=True, cmap='RdBu_r', plot_abs=False, symmetric_cbar=True,vmax=0.5, vmin=-0.5)
    correlation_file_name = './results/example/rac_dopamine_correlation_'+ dynamic_feature_name +'_stat_map_20voxel_cluster.nii.gz'
    correlation_img.to_filename(correlation_file_name)
