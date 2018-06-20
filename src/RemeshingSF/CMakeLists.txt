LIST ( APPEND RemeshingSF_files
BndryFacePtrFreez.cxx
BndryFacePtrFreez.hh
BndryNodePtr.cxx
BndryNodePtrFreez.cxx
BndryNodePtr.hh
BndryNodePtrFreez.hh
ChangeBndryPtr.cxx
ChangeBndryPtr.hh
CoIntrPnt.cxx
CoIntrPnt.hh
CoNorm4Ar.cxx
CoNorm4Ar.hh
CoNorm4Pg.cxx
CoNorm4Pg.hh
CoNorm4Pg3D.cxx
CoNorm4Pg3D.hh
CoNorm4TCneq.cxx
CoNorm4TCneq.hh
CoPntDispl.cxx
CoPntDispl.hh
CoPntDispl3D.cxx
CoPntDispl3D.hh
DummyRemeshing.cxx
DummyRemeshing.hh
FindBEdg.cxx
FindBEdg.hh
FndPhPs.cxx
FndPhPs.hh
FndPhPs3D.cxx
FndPhPs3D.hh
FndPhPs3D_Manual.cxx
FndPhPs3D_Manual.hh
FixMshSps.cxx
FixMshSps.hh
FixMshSps3D.cxx
FixMshSps3D.hh
RdDps.cxx
RdDps.hh
RdDpsEq.cxx
RdDpsEq.hh
RdDpsEq3D.cxx
RdDpsEq3D.hh
Rdshp.cxx
Rdshp.hh
ShpDpndnc.cxx
ShpDpndnc.hh
Smoothing3D.cxx
Smoothing3D.hh
)

LIST ( APPEND RemeshingSF_libs Framework SConfig MathTools )

SF_ADD_PLUGIN_LIBRARY ( RemeshingSF )

#SF_WARN_ORPHAN_FILES()
