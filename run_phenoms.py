import argparse
import pbr
import os
from glob import glob
from pbr.base import _get_output
import json, uuid
from nipype.utils.filemanip import load_json
heuristic = load_json(os.path.join(os.path.split(pbr.__file__)[0], "heuristic.json"))["filetype_mapper"]
from subprocess import Popen, PIPE, check_call, check_output
import pandas as pd
import shutil

base_dir = '/data/henry12/phenoms/IDIBAPS/MulipleMS_Retrospective/'

def rename(final_path):
    print(final_path, final_path.replace(" ","-"))
    os.rename(final_path, final_path.replace(" ","-"))

def get_les(t1_path):
    lesion = glob(t1_path + 'r*les*.nii.gz')
    if len(lesion) >= 1:
        lesion = lesion[0]
    else:
        lesion = ""
    return lesion

def get_t1(t1_path):
    t1_file = glob(t1_path + '/r*MPRAGE*.nii.gz')
    T1 = ""
    if len(t1_file) >= 1: #and not "les" in t1_file and not "brain_mask" in t1_file:
        for t1 in t1_file:
            if not "les" in t1 and not "brain" in t1:
                T1 = t1
    else:
        T1 = ""
    return T1

def get_brain(t1_path):
    BM = glob(t1_path + 'r*brain*.nii.gz')
    if len(BM) >= 1:
        for b in BM:
            if not "mask" in b:
                BM = b
    else:
        BM = ""
    return BM

def get_brain_mask(t1_path):
    BM = glob(t1_path + 'r*brain_mask*.nii.gz')
    if len(BM) >= 1:
        BM = BM[0]
    else:
        BM = ""
    return BM


def run_sienax(t1, les):
    sub = t1.split('/')[6] +'_'+ t1.split('/')[7]
    t1_out =  '/data/henry12/phenoms/IDIBAPS/MulipleMS_Retrospective/sienax_output/' + sub
    if not os.path.exists(t1_out + '/I_brain.nii.gz'):
        if os.path.exists(t1) and os.path.exists(les):
            cmd = ["sienax_optibet", t1, "-lm", les, "-r", "-d", "-o", t1_out ]
            print("sienax_optibet", t1, "-lm", les, "-r", "-d", "-o", t1_out )
        elif os.path.exists(t1):
            cmd = ["sienax_optibet", t1, "-r", "-d", "-o", t1_out ]
        else:
            cmd = ""
        print(cmd)
        try:
            Popen(cmd).wait()
        except:
            pass
    else:
        print(t1_out, "SIENAX FILE EXISTS")


def run_first(t1, in_BET):
    sub = t1.split('/')[6] +'_'+ t1.split('/')[7]
    odir =  '/data/henry12/phenoms/IDIBAPS/MulipleMS_Retrospective/first_output/' + sub
    if not os.path.exists(odir): os.mkdir(odir)

    fname = in_BET.split('/')[-1].split('.nii')[0]
    oname = os.path.join(odir, fname)

    omat = oname+'_to_std_sub'

    flirt_job = ['first_flirt', in_BET, omat]
    mat_name = omat+'.mat'

    if not os.path.exists(mat_name):
        print(flirt_job)
        check_call(flirt_job)

    imsfirst = []
    imscorr = []
    for s in ['L_Accu', 'L_Amyg', 'L_Caud', 'L_Hipp', 'L_Puta', 'L_Thal', 'R_Accu', 'L_Pall',
              'R_Amyg', 'R_Caud', 'R_Hipp', 'R_Pall', 'R_Puta', 'R_Thal', 'BrStem']:

        modelN = '336'
        bcorr = '1'
        intref = '0'
        if 'Accu' in s:
            nmodes = '50'
        elif 'Amyg' in s:
            nmodes = '50'
            intref = '1'
        elif 'Caud' in s:
            nmodes = '30'
            intref= '1'
        elif 'Hipp' in s:
            nmodes = '30'
            intref = '1'
        elif 'Late' in s:
            nmodes = '40'
            intref = '1'
        elif 'Pall' in s:
            nmodes = '40'
            bcorr = '0'
        elif 'Puta' in s:
            nmodes = '40'
            bcorr = '0'
        elif 'Thal' in s:
            nmodes = '40'
            bcorr = '0'
        elif 'Stem' in s:
            nomdes = '40'
        else:
            raise ValueError('The structure {} is not in the structure list'.format(s))

        imfirst = '{}-{}_first'.format(oname, s)
        imsfirst.append(imfirst)

        FSLDIR = '/netopt/fsl5'

        if intref == '0':
            if bcorr == '1':
                cmd = ['{}/bin/run_first'.format(FSLDIR), '-i', in_BET, '-t', mat_name,
                       '-n', nmodes, '-o', imfirst, '-m', '{}/data/first/models_{}_bin/{}_bin.bmv'.format(FSLDIR, modelN, s)]
            else:
                 cmd = ['{}/bin/run_first'.format(FSLDIR), '-i', in_BET, '-t', mat_name,
                       '-n', nmodes, '-o', imfirst, '-m', '{}/data/first/models_{}_bin/05mm/{}_05mm.bmv'.format(FSLDIR, modelN, s)]
        else:
            cmd = ['{}/bin/run_first'.format(FSLDIR), '-i', in_BET, '-t', mat_name,
                   '-n', nmodes, '-o', imfirst, '-m', '{}/data/first/models_{}_bin/intref_thal/{}.bmv'.format(FSLDIR, modelN, s), '-intref', '{}/data/first/models_{}_bin/05mm/{}_Thal_05mm.bmv'.format(FSLDIR, modelN, s.split('_')[0])]
        print('\nRunning Segmentation {} with: {}'.format(s,cmd))
        check_call(cmd)

        imcorr = '{}-{}_corr'.format(oname, s)
        imscorr.append(imcorr)

        btype = 'fast'
        if bcorr != '1': btype='none'
        cmd = ['{}/bin/first_boundary_corr'.format(FSLDIR), '-s', imfirst, '-o', imcorr,
               '-i', in_BET, '-b', btype]
        print(cmd)
        check_call(cmd)

    cmd = ['{}/bin/fslmerge'.format(FSLDIR), '-t', '{}_all_{}_firstsegs'.format(oname, btype)]


    for imCORR in imscorr:
        if os.path.exists(imCORR + '.nii.gz'):
            cmd.append(imCORR+'.nii.gz')
            print(cmd)

    print(cmd)
    check_call(cmd)

    cmd = ['{}/bin/fslmerge'.format(FSLDIR), '-t', '{}_all_{}_origsegs'.format(oname, btype)]
    for imname in imsfirst:
        if os.path.exists(imname + '.nii.gz'):
            cmd.append(imname+'.nii.gz')
    print(cmd)
    check_call(cmd)

    cmd =  ['{}/bin/first_mult_bcorr'.format(FSLDIR),'-i', t1, '-u', '{}_all_{}_origsegs'.format(oname, btype),
            '-c', '{}_all_{}_firstsegs'.format(oname, btype), '-o', '{}_all_{}_firstsegs'.format(oname, btype)]
    print('{}/bin/first_mult_bcorr'.format(FSLDIR),'-i', t1, '-u', '{}_all_{}_origsegs'.format(oname, btype),
            '-c', '{}_all_{}_firstsegs'.format(oname, btype), '-o', '{}_all_{}_firstsegs'.format(oname, btype))
    check_call(cmd)



def brain_Extraction(reorient):
    in_BET = reorient.replace(".nii","_optiBET_brain.nii")
    print("#############", in_BET)
    if not os.path.exists(in_BET):
        cmd = ['/netopt/share/bin/local/optiBET.sh', '-i', reorient]
        print(cmd)
        Popen(cmd).wait()
    return in_BET

def run_all(les, t1, brain, bm, reorient):
    print("****************", reorient)
    #in_BET = brain_Extraction(reorient)
    run_sienax(t1, les)
    #run_first(t1, in_BET)

def reorient_brain(t1):
    reorient = ""
    try:
        odir =  '/data/henry12/phenoms/IDIBAPS/MulipleMS_Retrospective/first_output/'
        sub = t1.split('/')[6] +'_'+ t1.split('/')[7] + '/'
        reorient = odir +sub+ t1.replace(".nii","_reorient.nii").split('/')[-1]
        print(reorient, "^^^^^^^^^^^^^^^^^^REORIENT^^^^^^^^^^^^^^^^^^^^^^^^^^^")
        cmd = ['fslreorient2std', t1, reorient]
        print(cmd)
        Popen(cmd).wait()
    except:
        pass
    return reorient

for N in os.listdir(base_dir):
    if N.startswith('N'):
        new_path = base_dir + N
        for files in os.listdir(new_path):
            final_path = new_path +'/'+ files
            #rename(final_path)
            if os.path.exists(final_path + '/T1MPRAGE/'):
                t1_path = final_path + '/T1MPRAGE/'
                #print(t1_path)
                les = get_les(t1_path)
                t1 = get_t1(t1_path)
                bm = get_brain_mask(t1_path)
                brain = get_brain(t1_path)
                reorient = reorient_brain(t1)
                if os.path.exists(t1):
                    run_all(les, t1, brain, bm, reorient)
                    print(t1)
















"""def brain_extraction(t1):


    in_BET = reorient.replace(".nii","_optiBET_brain.nii")
    print(in_BET, "in BET")
    if not os.path.exists(odir):
        os.mkdir(odir)
    if not os.path.exists(odir + sub):
        os.mkdir(odir+sub)
    print("reorient", reorient)
    #if not os.path.exists(reorient):

    #if not os.path.exists(in_BET):
    cmd = ['/netopt/share/bin/local/optiBET.sh', '-i', reorient]
    print(cmd)
    Popen(cmd).wait()
    return in_BET"""


