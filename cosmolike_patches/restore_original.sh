if [ -z ${ROOTDIR} ]; then
    echo "ROOTDIR env variable not defined!"
    exit 1
fi

SOURCE_DIR=${ROOTDIR}/projects/cs2-project/cosmolike_patches/original_files
DES_Y3_DIR=${ROOTDIR}/projects/des_y3
COSMOLIKE_DIR=${ROOTDIR}/external_modules/code/cosmolike

cp ${SOURCE_DIR}/_cosmolike_prototype_base.py ${DES_Y3_DIR}/likelihood/
cp ${SOURCE_DIR}/interface.cpp ${DES_Y3_DIR}/interface/
cp ${SOURCE_DIR}/generic_interface.cpp ${COSMOLIKE_DIR}/
cp ${SOURCE_DIR}/structs.h ${COSMOLIKE_DIR}/
cp ${SOURCE_DIR}/cosmo2D.c ${COSMOLIKE_DIR}/
cp ${SOURCE_DIR}/cosmo3D.h ${COSMOLIKE_DIR}/
cp ${SOURCE_DIR}/cosmo3D.c ${COSMOLIKE_DIR}/