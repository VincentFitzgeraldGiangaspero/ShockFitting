LIST ( APPEND WritingMeshSF_files
DummyWritingMesh.cxx
DummyWritingMesh.hh
WriteBackTriangle.cxx
WriteBackTriangle.hh
WriteSdwInfo.cxx
WriteSdwInfo.hh
WriteSdwInfo3D.cxx
WriteSdwInfo3D.hh
WriteTriangle.cxx
WriteTriangleFreez.cxx
WriteTriangleFreezedConnectivity.cxx
WriteTriangleFreezedConnectivity.hh
WriteTriangle.hh
WriteTetgen.hh
WriteTetgen.cxx
WriteTriangleFreez.hh
WriteBackTetgen.cxx
WriteBackTetgen.hh

)

LIST ( APPEND WritingMeshSF_libs Framework SConfig MathTools )

SF_ADD_PLUGIN_LIBRARY ( WritingMeshSF )

#SF_WARN_ORPHAN_FILES()
