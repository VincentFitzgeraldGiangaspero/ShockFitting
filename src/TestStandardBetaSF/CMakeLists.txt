LIST ( APPEND TestStandardBetaSF_files testStandardBetaSF.cxx )

# AL: the following line assumes that the module dirs (SF_MODULES_LIST) share the same name as the corresponding libs
LIST ( APPEND TestStandardBetaSF_libs ${SF_KERNEL_LIBS} ${SF_KERNEL_STATIC_LIBS}
CFDSolverSF
ConverterSF
CopyMakerSF
MeshGeneratorSF
RemeshingSF
ShockDetectorSF
StateUpdaterSF
VariableTransformerSF
WritingMeshSF
)

IF (SF_HAVE_TRILIBRARY)
LIST ( APPEND TestStandardBetaSF_libs triangleLib TriangleSF )
ENDIF()


SF_ADD_PLUGIN_APP ( TestStandardBetaSF )

###############################################################################

#SF_WARN_ORPHAN_FILES()
