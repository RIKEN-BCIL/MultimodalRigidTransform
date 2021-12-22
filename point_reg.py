#!/usr/bin/env python
#Do point registration

#$Id: point_reg.py,v 8bfe38a808d0 2017/05/11 02:07:07 watabe $

import os
import sys
import tempfile
import shutil
import optparse
import string
import numpy as np
from commands import getoutput

def centroid(X):
    """
    Calculate the centroid from a vectorset X
    """
    C = sum(X)/len(X)
    return C

def kabsch(P, Q):
    """
    The optimal rotation matrix U is calculated and then used to rotate matrix
    P unto matrix Q so the minimum root-mean-square deviation (RMSD) can be
    calculated.

    Using the Kabsch algorithm with two sets of paired point P and Q,
    centered around the center-of-mass.
    Each vector set is represented as an NxD matrix, where D is the
    the dimension of the space.

    The algorithm works in three steps:

    - a translation of P and Q
    - the computation of a covariance matrix C
    - computation of the optimal rotation matrix U

    http://en.wikipedia.org/wiki/Kabsch_algorithm

    Parameters:

    P -- (N, number of points)x(D, dimension) matrix
    Q -- (N, number of points)x(D, dimension) matrix

    Returns:

    U -- Rotation matrix
    
    """

    # Computation of the covariance matrix
    C = np.dot(np.transpose(P), Q)

    # Computation of the optimal rotation matrix
    # This can be done using singular value decomposition (SVD)
    # Getting the sign of the det(V)*(W) to decide
    # whether we need to correct our rotation matrix to ensure a
    # right-handed coordinate system.
    # And finally calculating the optimal rotation matrix U
    # see http://en.wikipedia.org/wiki/Kabsch_algorithm

    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    # Create Rotation matrix U
    U = np.dot(V, W)
    
    return U


if __name__ == '__main__':

    option_parser = optparse.OptionParser()
    option_parser.add_option('--pet',help='PET file (assume input)',type='str',default='')
    option_parser.add_option('--mri',help='MRI file (assume reference)',type='str',default='')
    option_parser.add_option('--save',help='Save reoriented file from PET to MRI',type='str',default='')
    option_parser.add_option('--psize',help='point size(mm)',type='float',default=15)
    option_parser.add_option('--mricoord',help='Saved MRI Coordinate file. Default mri.coord',type='str',default='mri.coord')
    option_parser.add_option('--mricfile',help='MRI Coordinate file. Read from file instead of manual type',type='str',default='')
    option_parser.add_option('--petcoord',help='Saved PET Coordinate file. Default pet.coord',type='str',default='pet.coord')
    option_parser.add_option('--petcfile',help='PET Coordinate file. Read from file instead of manual type',type='str',default='')
    option_parser.add_option('--inmat',help='Input Transformation matrix from PET to MRI. If you give it, we apply without anything',type='str',default='')
    option_parser.add_option('--outmat',help='Transformation matrix from PET to MRI. Default out.mat',type='str',default='out.mat')
    option_parser.add_option('--noview',help='No FSLView',action='store_true',default=False)
    option_parser.add_option('--nodel',help='Not delete temp files',action='store_true',default=False)
    option_parser.add_option('--norevx',help='Not Reverse X direction',action='store_true',default=False)
    option_parser.add_option('--tpet',help='Threshold Value for PET. Defalut is no threshold',type='float',default=0)
    option_parser.add_option('--tmri',help='Threshold Value for MRI. Defalut is no threshold',type='float',default=0)
    option_parser.add_option('--kabsch',help='Use Kabsch algorithm instead of pointflirt',action='store_true',default=False)
    options,args = option_parser.parse_args()

    if len(sys.argv)==1:

        option_parser.print_help()
        sys.exit()


    petorig = options.pet
    alldims=getoutput("fslsize "+petorig+" -s")
    listdims=alldims.split()
    petdx=float(listdims[12])
    petdy=float(listdims[14])
    petdz=float(listdims[16])
    petsx=float(listdims[2])
    petsy=float(listdims[4])
    petsz=float(listdims[6])
    mriorig = options.mri
    alldims=getoutput("fslsize "+mriorig+" -s")
    listdims=alldims.split()
    mridx=float(listdims[12])
    mridy=float(listdims[14])
    mridz=float(listdims[16])
    mrisx=float(listdims[2])
    mrisy=float(listdims[4])
    mrisz=float(listdims[6])
    mrisform = getoutput("fslorient -getsform %s" % mriorig)
    saveimg = options.save
    psize = options.psize
    mricoord = options.mricoord
    petcoord = options.petcoord
    outmat = options.outmat
    inmat = options.inmat
    noview_flag = options.noview
    mricfile = options.mricfile
    petcfile = options.petcfile
    nodel_flag = options.nodel
    norevx_flag = options.norevx
    thresh_pet = options.tpet
    thresh_mri = options.tmri
    kabsch_flag = options.kabsch

    if len(saveimg)==0 or len(petorig)==0 or len(mriorig)==0:
        option_parser.print_help()
        sys.exit()

    tmpf = tempfile.mktemp()
    petdelori = "%s_petdelori.nii.gz" % tmpf
    mridelori = "%s_mridelori.nii.gz" % tmpf

    shutil.copy(petorig,petdelori)
    shutil.copy(mriorig,mridelori)

    os.system("fslorient -deleteorient %s" % petdelori)
    os.system("fslorient -deleteorient %s" % mridelori)

    #swap xdim to fix PET and MRI images are radiological and neurological

    #if norevx_flag is False:
        #os.system("fslswapdim %s -x y z %s" % (petdelori,petdelori))

    #determine points

    if len(mricfile)>0 or len(inmat)>0: #already define
        fp = open(mricfile,"rt")
        lines = fp.readlines()
        no_point = len(lines)
        mricoord = mricfile

    else:
        if noview_flag is not True:
            os.system("fslview %s &" % (mridelori))

        no_point = 1
        mripx = []
        mripy = []
        mripz = []
        while 1:

            print "Please select No %d points (type # if you finish) in MRI(reference) image" % no_point
            out = raw_input("center of X,Y,Z (mm): ")

            if out=="#":
                break
            else:
                o2 = out.split()
                mripx.append(float(o2[0]))
                mripy.append(float(o2[1]))
                mripz.append(float(o2[2]))
                no_point = no_point + 1

        no_point = no_point - 1
        mripx = np.array(mripx)
        mripy = np.array(mripy)
        mripz = np.array(mripz)

        for i in range(no_point):
            os.system("fslmaths %s -roi %d %d %d %d %d %d 0 1 %s_mri%d.nii.gz" % (mridelori,int((mripx[i] - psize*0.5)/mridx),int(psize/mridx),int((mripy[i] - psize*0.5)/mridy),int(psize/mridy),int((mripz[i] - psize*0.5)/mridz),int(psize/mridz),tmpf,i+1))
            if thresh_mri>0:
                os.system("fslmaths %s_mri%d.nii.gz -thr %g %s_mri%d.nii.gz" % (tmpf,i+1,thresh_mri,tmpf,i+1))
        #summary all points
        mriout = []

        fpmri = open(mricoord,"wt")
        for ip in range(no_point):
            output = getoutput("fslstats %s_mri%d.nii.gz -c" % (tmpf,ip+1))
            fpmri.write(output+"\n")
        fpmri.close()

    #for pet

    if len(petcfile)>0 or len(inmat)>0: #already define
        petcoord = petcfile
    else:

        if noview_flag is not True:
            os.system("fslview %s &" % (petdelori))
        petpx = []
        petpy = []
        petpz = []

        for ip in range(no_point):
            print "Please select No %d points (Total %d) in PET(target) image" % (ip + 1,no_point)
            out = raw_input("center of X,Y,Z (mm): ")
            o2 = out.split()
            petpx.append(float(o2[0]))
            petpy.append(float(o2[1]))
            petpz.append(float(o2[2]))

        petpx = np.array(petpx)
        petpy = np.array(petpy)
        petpz = np.array(petpz)

        for i in range(no_point):
            os.system("fslmaths %s -roi %d %d %d %d %d %d 0 1 %s_pet%d.nii.gz" % (petdelori,int((petpx[i] - psize*0.5)/petdx),int(psize/petdx),int((petpy[i] - psize*0.5)/petdy),int(psize/petdy),int((petpz[i] - psize*0.5)/petdz),int(psize/petdz),tmpf,i+1))
            if thresh_pet>0:
                os.system("fslmaths %s_pet%d.nii.gz -thr %g %s_pet%d.nii.gz" % (tmpf,i+1,thresh_pet,tmpf,i+1))

        #summary all points

        petout = []
        fppet = open(petcoord,"wt")
        for ip in range(no_point):
            output = getoutput("fslstats %s_pet%d.nii.gz -c" % (tmpf,ip+1))
            fppet.write(output+"\n")
        fppet.close()



    #print "pointflirt -i %s -r %s -o %s" % (petcoord,mricoord,outmat)

    if kabsch_flag is True:
        fppet = open(petcoord,"rt")
        petxyz = fppet.readlines()
        fppet.close()
        fpmri = open(mricoord,"rt")
        mrixyz = fpmri.readlines()
        fpmri.close()

        pet_x = []
        pet_y = []
        pet_z = []

        for ipetxyz in petxyz:
            sp = ipetxyz.split()
            pet_x.append(float(sp[0]))
            pet_y.append(float(sp[1]))
            pet_z.append(float(sp[2]))

        mri_x = []
        mri_y = []
        mri_z = []

        for imrixyz in mrixyz:

            sp = imrixyz.split()
            mri_x.append(float(sp[0]))
            mri_y.append(float(sp[1]))
            mri_z.append(float(sp[2]))



        P = np.zeros((len(pet_x),3))
        Q = np.zeros((len(pet_x),3))
        
        for i in range(len(pet_x)):
            P[i,0] = pet_x[i]
            P[i,1] = pet_y[i]
            P[i,2] = pet_z[i]
            Q[i,0] = mri_x[i]
            Q[i,1] = mri_y[i]
            Q[i,2] = mri_z[i]

        Pc = centroid(P)
        Qc = centroid(Q)

        P -= Pc
        Q -= Qc

        U = kabsch(P, Q)
        transp = np.array([[1,0,0,-Pc[0]],[0,1,0,-Pc[1]],[0,0,1,-Pc[2]],[0,0,0,1]])
        transq = np.array([[1,0,0,Qc[0]],[0,1,0,Qc[1]],[0,0,1,Qc[2]],[0,0,0,1]])

        U4x4 = np.zeros((4,4),'f')
        U4x4[0:3,0:3] = U[:,:]
        U4x4[3,3] = 1.0
        tU4x4 = np.transpose(U4x4)
        transmat = np.dot(transq,np.dot(tU4x4,transp))

        fp = open(outmat,"wt")
        fp.write("%g %g %g %g\n" % (transmat[0,0],transmat[0,1],transmat[0,2],transmat[0,3]))
        fp.write("%g %g %g %g\n" % (transmat[1,0],transmat[1,1],transmat[1,2],transmat[1,3]))
        fp.write("%g %g %g %g\n" % (transmat[2,0],transmat[2,1],transmat[2,2],transmat[2,3]))
        fp.write("%g %g %g %g\n" % (transmat[3,0],transmat[3,1],transmat[3,2],transmat[3,3]))
        fp.close()

    else:
        res = os.system("pointflirt -i %s -r %s -o %s" % (petcoord,mricoord,outmat))

    if len(inmat)>0:
        os.system("flirt -in %s -out %s -ref %s -applyxfm -init %s" % (petdelori,saveimg,mridelori,inmat))
    else:
        os.system("flirt -in %s -out %s -ref %s -applyxfm -init %s" % (petdelori,saveimg,mridelori,outmat))

    #bring back original sform
    os.system("fslorient -setsformcode 1 %s > /dev/null 2>&1" % (saveimg))

    #print("fslorient -setsform %s %s" % (mrisform,saveimg))
    os.system("fslorient -setsform %s %s > /dev/null 2>&1" % (mrisform,saveimg))
    os.system("fslorient -copysform2qform %s > /dev/null 2>&1" % (saveimg))

    #delete related files

    if nodel_flag is not True:
        print "Removing tempfile %s" % tmpf
        
        for ip in range(no_point):
            try:
                os.unlink("%s_pet%d.nii.gz" % (tmpf,ip+1))
            except:
                pass
            try:
                os.unlink("%s" % petdelori)
            except:
                pass
            try:
                os.unlink("%s" % mridelori)
            except:
                pass
