#!/bin/bash
PROG_DIR=/home/calas/Carlos_Membrane_project/BoundaryElementMethod/clovers/precVaried/
OUT_DIR=/home/calas/Carlos_Membrane_project/BoundaryElementMethod/elementary_plots/clover5_GvsEpsilon/
PROG_NAME=clover5eps030
OUT_NAME=${PROG_NAME}.out
export OMP_NUM_THREADS=8
export OMP_WAIT_POLICY=active
export OMP_DYNAMIC=false
##export OMP_PROC_BIND=true
export LD_LIBRARY_PATH="/home/calas/anaconda3/envs/arblib_env/lib"
{ ${PROG_DIR}${PROG_NAME}; } |& tee ${OUT_DIR}${OUT_NAME}
##export OMP_NUM_THREADS=2
##{ ${PROG_DIR}${PROG_NAME}; } |& tee ${OUT_DIR}${OUT_NAME}
##export OMP_NUM_THREADS=3
##{ ${PROG_DIR}${PROG_NAME}; } |& tee ${OUT_DIR}${OUT_NAME}
##export OMP_NUM_THREADS=4
##{ ${PROG_DIR}${PROG_NAME}; } |& tee ${OUT_DIR}${OUT_NAME}
##export OMP_NUM_THREADS=5
##{ ${PROG_DIR}${PROG_NAME}; } |& tee ${OUT_DIR}${OUT_NAME}
##export OMP_NUM_THREADS=6
##{ ${PROG_DIR}${PROG_NAME}; } |& tee ${OUT_DIR}${OUT_NAME}
##export OMP_NUM_THREADS=7
##{ ${PROG_DIR}${PROG_NAME}; } |& tee ${OUT_DIR}${OUT_NAME}
##export OMP_NUM_THREADS=8
##{ ${PROG_DIR}${PROG_NAME}; } |& tee ${OUT_DIR}${OUT_NAME}
##export OMP_NUM_THREADS=9
##{ ${PROG_DIR}${PROG_NAME}; } |& tee ${OUT_DIR}${OUT_NAME}
##export OMP_NUM_THREADS=10
##{ ${PROG_DIR}${PROG_NAME}; } |& tee ${OUT_DIR}${OUT_NAME}
##export OMP_NUM_THREADS=11
##{ ${PROG_DIR}${PROG_NAME}; } |& tee ${OUT_DIR}${OUT_NAME}
##export OMP_NUM_THREADS=12
##{ ${PROG_DIR}${PROG_NAME}; } |& tee ${OUT_DIR}${OUT_NAME}
##export OMP_NUM_THREADS=13
##{ ${PROG_DIR}${PROG_NAME}; } |& tee ${OUT_DIR}${OUT_NAME}
##export OMP_NUM_THREADS=14
##{ ${PROG_DIR}${PROG_NAME}; } |& tee ${OUT_DIR}${OUT_NAME}
##export OMP_NUM_THREADS=15
##{ ${PROG_DIR}${PROG_NAME}; } |& tee ${OUT_DIR}${OUT_NAME}
##export OMP_NUM_THREADS=16
##{ ${PROG_DIR}${PROG_NAME}; } |& tee ${OUT_DIR}${OUT_NAME}

