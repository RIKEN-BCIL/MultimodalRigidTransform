# MultimodalRigidTransform

*Point registration*
Usage:
point_reg.py [options]

Options:
  -h, --help           show this help message and exit
  --pet=PET            PET file (assume input)
  --mri=MRI            MRI file (assume reference)
  --save=SAVE          Save reoriented file from PET to MRI
  --psize=PSIZE        point size(mm)
  --mricoord=MRICOORD  Saved MRI Coordinate file. Default mri.coord
  --mricfile=MRICFILE  MRI Coordinate file. Read from file instead of manual
                       type
  --petcoord=PETCOORD  Saved PET Coordinate file. Default pet.coord
  --petcfile=PETCFILE  PET Coordinate file. Read from file instead of manual
                       type
  --inmat=INMAT        Input Transformation matrix from PET to MRI. If you
                       give it, we apply without anything
  --outmat=OUTMAT      Transformation matrix from PET to MRI. Default out.mat
  --noview             No FSLView
  --nodel              Not delete temp files
  --checkorient        Check Orientation Neurological/Radiological
  --norevx             Not Reverse X direction
  --tpet=TPET          Threshold Value for PET. Defalut is no threshold
  --tmri=TMRI          Threshold Value for MRI. Defalut is no threshold
  --kabsch             Use Kabsch algorithm instead of pointflirt
