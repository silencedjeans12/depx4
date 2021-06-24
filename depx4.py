# Aca van a ir todos tus metodos Andres
def SPELI4(IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,D,N,NBDCND,BDC,BDD,COFX,AN,BN,CN,DN,UN,ZN,AM,BM,CM,DM,UM,ZM,GRHS,USOL,IDMN,W,PERTRB,IERROR):
#    
#     SPELI4 SETS UP VECTORS AND ARRAYS FOR INPUT TO BLKTRI
#     AND COMPUTES A SECOND ORDER SOLUTION IN USOL.  A RETURN JUMP TO
#     SEPX4 OCCURRS IF IORDER=2.  IF IORDER=4 A FOURTH ORDER
#     SOLUTION IS GENERATED IN USOL.
#
    global KSWX,KSWY,K,L
    global AIT,BIT,CIT,DIT
    global MIT,NIT,IS,MS
    global JS,NS,DLX,DLY
    global TDLX3,TDLY3,DLX4,DLY4
#
#
#   SET PARAMETERS INTERNALLY
#
    KSWX=MBDCND+1
    KSWY=NBDCND+1
    K=M+1
    L=N+1
    AIT = A
    BIT = B
    CIT = C
    DIT = D
    DLY=(DIT-CIT)/np.float32(N)
    for I in range(2,M):
        for J in range(2,N):
            USOL[I][J]=DLY**2**GRHS[I][J]
    if (KSWX!=2 or KSWX!=3):
        for J in range(2,N):
            USOL[1][J]=DLY**2**GRHS[1][J]
    if (KSWX!=2 or KSWX!=5):
        for J in range(2,N):
            USOL[K][J]=DLY**2**GRHS[K][J]
    if (KSWY!=2 or KSWY!=3):
        for I in range(2,M):
            USOL[I][1]=DLY**2**GRHS[I][1]
    if (KSWY!=2 or KSWY!=5):
        for I in range(2,M):
            USOL[I][L]=DLY**2**GRHS[I][L]
    if (KSWX!=2 and KSWX!=3 and KSWY!=2 and KSWX!=3):
        USOL[1][1]=DLY**2*GRHS[1][1]
    if (KSWX!=2 and KSWX!=5 and KSWY!=2 and KSWX!=3):
        USOL[K][1]=DLY**2*GRHS[K][1]
    if (KSWX!=2 and KSWX!=3 and KSWY!=2 and KSWX!=5):
        USOL[1][L]=DLY**2*GRHS[1][L]
    if (KSWX!=2 and KSWX!=5 and KSWY!=2 and KSWX!=5):
        USOL[K][L]=DLY**2*GRHS[K][L]
    I1=1
#
#   SET SWITCHES FOR PERIODIC OR NON-PERIODIC BOUNDARIES
#
    MP=1
    if(KSWX==1):
        MP=0
    NP=NBDCND
#
#     SET DLX,DLY AND SIZE OF BLOCK TRI-DIAGONAL SYSTEM GENERATED
#     IN NINT,MINT
#
    DLX=(BIT-AIT)/float(M)
    MIT=K-1
    if(KSWX==2):
        MIT=K-2
    if(KSWX==4):
        MIT=K
    DLY=(DIT-CIT)/float(N)
    NIT=L-1
    if(KSWY==2):
        NIT=L-2
    if(KSWY==4):
        NIT=L
    TDLX3=float(2.000)*DLX**3
    DLX4=DLX**4
    TDLY3=float(2.000)*DLY**3
    DLY4=DLY**4
#
#     SET SUBSCRIPT LIMITS FOR PORTION OF ARRAY TO INPUT TO BLKTRI
#    
    IS=1
    JS=1
    if(KSWX==2 or KSWX==3):
        IS=2
    if(KSWY==2 or KSWY==3):
        JS=2
    NS=NIT+JS-1
    MS=MIT+IS-1
#
#     SET X - DIRECTION
#
    for I in range(1,MIT):
        XI=AIT+float(IS+I-2)*DLX
        COFX(XI,AI,BI,CI)
        AXI=(AI/DLX-float(0.50)*BI)/DLX
        BXI=float(-2.000)*AI/DLX**2+CIT
        CXI=(AI/DLX+float(0.50)*BI)/DLX
        AM[I]=DLY**2*AXI
        BM[I]=DLY**2*BXI
        CM[I]=DLY**2*CXI
#
#     SET Y DIRECTION
#
    for J in range(1,NIT):
        DYJ=float(1.00)
        EYJ=float(-2.00)
        FYJ=float(1.00)
        AN[J]=DYJ
        BN[J]=EYJ
        CN[J]=FYJ
#
#     ADJUST EDGES IN X DIRECTION UNLESS PERIODIC
#
    AX1=AM[1]
    CXM=CM[MIT]
#
#     DIRICHLET-DIRICHLET IN X DIRECTION
#
    AM[1]=float(0.00)
    CM[MIT]=float(0.00)
#
#     MIXED-DIRICHLET IN X DIRECTION
#
    AM[1]=float(0.00)
    BM[1]=BM[1]+float(2.00)*ALPHA*DLX*AX1
    CM[1]=CM[1]+AX1
    CM[MIT]=float(0.00)
#
#     DIRICHLET-MIXED IN X DIRECTION
#
    AM[1]=float(0.00)
    AM[MIT]=AM[MIT]+CXM
    BM[MIT]=BM[MIT]-float(2.00)*BETA*DLX*CXM
    CM[MIT]=float(0.00)
#
#     MIXED - MIXED IN X DIRECTION
#
    AM[1]=float(0.00)
    BM[1] = BM[1]+float(2.00)*DLX*ALPHA*AX1
    CM[1] = CM[1]+AX1
    AM[MIT] = AM[MIT]+CXM
    BM[MIT] = BM[MIT]-float(2.00)*DLX*BETA*CXM
    CM[MIT] = float(0.00)
