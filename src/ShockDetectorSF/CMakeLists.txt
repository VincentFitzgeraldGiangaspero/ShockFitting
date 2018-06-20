LIST ( APPEND ShockDetectorSF_files
ComputeGradient.cxx
ComputeGradient.hh
ComputeShockPointsState.cxx
ComputeShockPointsState.hh
DensityGradientMaximaPg2D.cxx
DensityGradientMaximaPg2D.hh
DetectorAlgorithm.cxx
DetectorAlgorithm.hh
DummyShockDetector.cxx
DummyShockDetector.hh
FindShockSurfaceNormals2D.cxx
FindShockSurfaceNormals2D.hh
GnoffoShockSensor.cxx
GnoffoShockSensor.hh
IdentifySpecialPoints.cxx
IdentifySpecialPoints.hh
NormalMachNumberPg2D.cxx
NormalMachNumberPg2D.hh
ShockPointProcessingCGAL.cxx
ShockPointProcessingCGAL.hh
FindShockSurfaceNormals3DCGAL.hh
FindShockSurfaceNormals3DCGAL.cxx
ShockSurfaceRemeshingMMG.hh
ShockSurfaceRemeshingMMG.cxx
ComputeShockPointsState3D.hh
ComputeShockPointsState3D.cxx
LaplacianSmoothing3D.hh
LaplacianSmoothing3D.cxx
AssignShockPointState3D.cxx
AssignShockPointState3D.hh
)

LIST ( APPEND ShockDetectorSF_libs Framework SConfig MathTools )

SF_ADD_PLUGIN_LIBRARY ( ShockDetectorSF )


#example for CGAL (not necessarily to be put here)





#SF_WARN_ORPHAN_FILES()
