MODULE constantes
  INTEGER, PARAMETER :: mi=181, nj=76, itermax=1600000, nsolid = 7
  INTEGER, PARAMETER :: SGL=SELECTED_REAL_KIND(P=6,R=37)
  INTEGER, PARAMETER :: DBL=SELECTED_REAL_KIND(P=15,R=300)
  REAL(kind=DBL),    PARAMETER :: cero=0._DBL, visc_solido = 1.e40_DBL
END MODULE constantes
