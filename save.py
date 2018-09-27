import os

# outputname xl_file input_pk Om w0 h0
xl_exe = './Flender_XSB_profile/save_xl'
yl_exe = './Flender_pressure_profile/save_yl'
input_pk = './input_pk/wmap9_fid_matterpower_z0.dat'
Om = '0.279'
w0 = '-1'
h0 = '0.70'

param_dict = {'fid':[''], 'agn':['1e-7','1e-5'], 'ac':['0.8','1.2'], 'alpha':['0.3','0'], 'beta':['-1','1'], 'gamma0':['0.01','1.0','1.67'], 'zslope':['0.68','2.72']}
#param_dict = {'agn':['1e-7','1e-5'], 'ac':['0.8','1.2'], 'alpha':['0.3','0'], 'beta':['-1','1'], 'gamma0':['0.01','1.0','1.67'], 'zslope':['0.68','2.72']}

for param in param_dict.keys():
    for val in param_dict[param] :
        if param == 'fid':
            name = param
        else: 
            name = param+'_'+val
        inputname = './params_xl/param_'+name
        command = xl_exe+' '+inputname+' 1'
        os.system(command)
        inputname = './params_yl/param_'+name
        command = yl_exe+' '+inputname
        os.system(command)

