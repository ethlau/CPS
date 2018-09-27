import os

# outputname xl_file input_pk Om w0 h0
exe = './cross_corr/halo_model'
exe_cut = './cross_corr_cut/halo_model'
input_pk = 'input_pk/wmap9_fid_matterpower_z0.dat'
Om = '0.279'
w0 = '-1'
h0 = '0.70'
z_source = '1'

#./halo_model test.cl ../Flender_XSB_profile/xl.dat ../Flender_pressure_profile/yl.dat ../input_pk/wmap9_fid_matterpower_z0.dat 0.279 -1 0.70 0 z_source


param_dict = {'fid':[''], 'agn':['1e-7','1e-5'], 'ac':['0.8','1.2'], 'alpha':['0.3','0'], 'beta':['-1','1'], 'gamma0':['0.01','1.0','1.67'], 'zslope':['0.68','2.72']}
#param_dict = {'fid':[''], 'agn':['1e-7','1e-5'], 'ac':['0.8','1.2'], 'alpha':['0.3','0'], 'beta':['-1','1'], 'gamma0':['0.01','1.0','1.67'], 'zslope':['0.68','2.72']}
#param_dict = {'fid':['']}
for param in param_dict.keys():
    for val in param_dict[param] :
        if param == 'fid':
            name = param
        else: 
            name = param+'_'+val
 
        outputname = './cl/'+name+'.cl'
        xl_file = './xl_'+name+'.dat'
        yl_file = './yl_'+name+'.dat'

        command = exe+' '+outputname+' '+xl_file+' '+yl_file+' '+input_pk+' '+Om+' '+w0+' '+h0+' 0 '+z_source
        os.system(command)

        #outputname_cut = './cl/'+name+'.cl'
        #command = exe+' '+outputname_cut+' '+xl_file+' '+yl_file+' '+input_pk+' '+Om+' '+w0+' '+h0+' 0 '+z_source
        #os.system(command)

 
