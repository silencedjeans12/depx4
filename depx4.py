# Aca van a ir todos tus metodos Andres
from os import mkdir, rename


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
        USOL[I][J] = DLY**2*GRHS[I][J]
        #CONTINUE
    #20
    if (KSWX != 2 or KSWX != 3):#---->40
        for J in range(2,N):
            USOL[1][J] = DLY**2*GRHS[1][J]
    #30
    #40
    if(KSWX != 2 or KSWX != 5):#--->60
        for J in range(2,N):
            USOL[K][J] = DLY**2*GRHS[K][J]
    #50
    #60
    if(KSWY != 2 or KSWY != 3):#--->80
        for I in range(2,M):
            USOL[I][1]=DLY**2*GRHS[I][1]
    #70
    #80
    if(KSWY != 2 or KSWY != 5):#---->100
        for I in range(2,M):
            USOL[I][L] = DLY**2*GRHS[I][L]
    #90
    #100
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
    if(KSWX == 2):
        MIT=K-2
    if(KSWX == 4):
        MIT=K
    DLY=(DIT-CIT)/float(N)
    NIT=L-1
    if(KSWY == 2):
        NIT=L-2
    if(KSWY == 4):
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
    #110
#     SET Y DIRECTION
#
    for J in range(1,NIT):
        DYJ=float(1.00)
        EYJ=float(-2.00)
        FYJ=float(1.00)
        AN[J]=DYJ
        BN[J]=EYJ
        CN[J]=FYJ
    #120
    
#
#     ADJUST EDGES IN X DIRECTION UNLESS PERIODIC
#
    AX1=AM[1]
    CXM=CM[MIT]
    if(KSWX == 1):#---170
        #170
        DY1 = AN[1]
        FYN = CN[NIT]
        GAMA = float(0.00)
        XNU = float(0.00)
    
    if(KSWX == 2):#--->130
        AM[1] = float(0.000)
        CM[MIT] = float(0.000)
        #---->170
        #170
        DY1 = AN[1]
        FYN = CN[NIT]
        GAMA = float(0.00)
        XNU = float(0.00)

    if(KSWX == 3):#--->150
        AM[1] = float(0.00)
        AM[MIT] = AM[MIT]+CXM
        BM[MIT] = BM[MIT]-float(2.00)*BETA*DLX*CXM
        CM[MIT] = float(0.00)
        #--->170
        #170
        DY1 = AN[1]
        FYN = CN[NIT]
        GAMA = float(0.00)
        XNU = float(0.00)
    

    if(KSWX == 4):#--->160
        #160
        AM[1] = float(0.00)
        BM[1] = BM[1]+float(2.00)*DLX*ALPHA*AX1
        CM[1] = CM[1] + AX1
        AM[MIT] =AM[MIT]+CXM
        BM[MIT] = BM[MIT]-float(2.00)*DLX*BETA*CXM
        CM[MIT] = float(0.00)
        #--->170
        #170
        DY1 = AN[1]
        FYN = CN[NIT]
        GAMA = float(0.00)
        XNU = float(0.00)

    if(KSWX == 5):#---->140
        AM[1] = float(0.00)
        BM[1] = BM[1]+float(2.00)*ALPHA*DLX*AX1
        CM[1] = CM[1]+AX1
        CM[MIT] = float(0.00)
        #--->170
        #170
        DY1 = AN[1]
        FYN = CN[NIT]
        GAMA = float(0.00)
        XNU = float(0.00)

    #170
    #GO TO (220,180,200,210,190),KSWY
    if(KSWY != 1):#--->220
        
        if(KSWY == 2):#--->180
            AN[1] = float(0.00)
            CN[NIT] = float(0.00)
            #--->220
        if(KSWY == 3):#---->200
            AN[1] = float(0.00)
            AN[1] = float(0.00)
            BN[NIT] = BN[NIT]-float(2.00)*DLY*XNU*FYN
            CN[NIT] = float(0.00)
            #---->220
        if(KSWY == 4):#--->210
            AN[1] = float(0.00)
            BN[1] = BN[1]+float(2.00)*DLY*GAMA*DY1
            CN[1] = CN[1]+DY1
            AN[NIT] = AN[NIT]*FYN
            BN[NIT] = BN[NIT]-float(2.00)*DLY*XNU*FYN
            CN[NIT] = float(0.00)
            #220
        if(KSWY == 5):#--->190
            AN[1] = float(0.00)
            BN[1] = BN[1]+float(2.00)*DLY*GAMA*DY1
            CN[1] = CN[1] * DY1
            CN[NIT] = float(0.00)
            #--->220
    #220
    if(KSWX != 1):#--->270
        #270
        for J in range(JS,NS):
            if(KSWX == 2 and KSWX == 3):#--->230
                #230
                USOL[IS][J] = USOL[IS][J] - AX1*USOL[1][J]
                #---->240
                
            else:
                #230
                USOL[IS][J] = USOL[IS][J]+float(2.00)*DLX*AX1*BDA[J]
                #--->240
                
            #240
            if(KSWY != 2 and KSWY != 5):
                #250
                USOL[MS][J] = USOL[MS][J]-float(2.00)*DLX*CXM*BDB[J]
                #--->260
            else:
                USOL[MS][J] = USOL[MS][J]-CXM*USOL[K][J]
                #--->260
            #260
            #CONTINUE
    #270
    if(KSWY != 1):#--->320
        for I in range(IS,MS):
            if(KSWY != 2 and KSWY != 3):
                USOL[I][JS] = USOL[I][JS]+float(2.00)*DLY*DY1*BDC[1]
            else:
                USOL[I][JS] = USOL[I][JS]-DY1*USOL[I][1]
            if(KSWY != 2 and KSWY != 5):
                USOL[I][NS] = USOL[I][NS]-float(2.00)*DLY*FYN*BDD[I]
            else:
                USOL[I][NS] = USOL[I][NS]-FYN*USOL[I][L]
    
    if(IORDER == 4):
        for J in range(JS,NS):
            GRHS[IS][J] = USOL[IS][J]
            GRHS[MS,J] = USOL[MS][J]
        for I in range(IS,MS):
            GRHS[I][JS] = USOL[I][JS]
            GRHS[I][NS] = USOL[I][NS]
    IORD = IORDER
    PERTRB = float(0.00)

    CHKSN4(MBDCND,NBDCND,ALPHA,BETA,COFX,SINGLR)

    if(SINGLR):
        TRIS4(MIT,AM,BM,CM,DM,UM,ZM)
    if(SINGLR):
        TRIS4(NIT,AN,BN,CN,DN,UN,ZN)

#   ADJUST RIGHT HAND SIDE IF NECESSARY

    if(SINGLR):
        ORTHO4(USOL,IDMN,ZN,ZM,PERTRB)

#    COMPUTE SOLUTION

#   SAVE ADJUSTED RIGHT HAND SIDE IN GRHS

    for J in range(JS,NS):
        for I in range(IS,MS):
            GRHS[I][J] = USOL[I][J]

    GENBUN(NP,NIT,MP,MIT,AM,BM,CM,IDMN,USOL[IS][JS],IERROR,W)

#    CHECK IF ERROR DETECTED IN POIS  
#   THIS CAN ONLY CORRESPOND TO IERROR=12

    if(IERROR != 0):
#   SET ERROR FLAG IF IMPROPER COEFFICIENTS INPUT TO POIS
        IERROR=12
        return
    if(IERROR == 0):
        return

#   SET PERIODIC BOUNDARIES IF NECESSARY
    if(KSWX == 1):
        for J in range(1,L):
            USOL[K][J] = USOL[1][J]
    if(KSWY == 1):
        for I in range(1,K):
            USOL[I][L] = USOL[I][1]

#   MINIMIZE SOLUTION WITH RESPECT TO WEIGHTED LEAST SQUARES
#   NORM IF OPERATOR IS SINGULAR

    if(SINGLR):
        MINSO4(USOL,IDMN,ZN,ZM,PERTRB)

#   RETURN IF DEFE4RED CORRECTIONS AND A FOURTH ORDER SOLUTION ARE
#   NOT FLAGGED

    if(IORD != 2):
        IORD = 2

    DEFE4(COFX,IDMN,USOL,GRHS)


#DECK TRIS4

def TRIS4 (N, A, B, C, D, U, Z):
#    C***BEGIN PROLOGUE  TRIS4
#C***SUBSIDIARY
#C***PURPOSE  Subsidiary to SEPX4
#C***LIBRARY   SLATEC
#C***TYPE      SINGLE PRECISION (TRIS4-S)
#C***AUTHOR  (UNKNOWN)
#C***DESCRIPTION
#C
#C     This subroutine solves for a non-zero eigenvector corresponding
#C     to the zero eigenvalue of the transpose of the rank
#   C     superdiagonal C , with A(1) in the (1,N) position, with
#   C     C(N) in the (N,1) position, AND all other elements zero.
#   C
#    C***SEE ALSO  SEPX4
#    C***ROUTINES CALLED  (NONE)
#    C***REVISION HISTORY  (YYMMDD)
#    C   801001  DATE WRITTEN
#    C   890831  Modified array declarations.  (WRB)
#    C   891214  Prologue converted to Version 4.0 format.  (BAB)
#    C   900402  Added TYPE section.  (WRB)
#    C***END PROLOGUE  TRIS4
#C
    

    BN= B[N]
    D[1] =  A[2]/B[1]
    V= A[1]
    U[1]= C[N]/B[1]
    NM2 = N - 2

    for J in range (2,NM2):
        DEN =B[J]-C(J-1)*D(J-1)
        D[J] = A(J+1)/DEN
        U[J] = -C(J-1)*U(J-1)/DEN
        BN = BN-V*U(J-1)
        V = -V*D(J-1)

    DEN = B(N-1)-C(N-2)*D(N-2)
    D(N-1) = (A(N)-C(N-2)*U(N-2))/DEN
    AN =  C(N-1)-V*D(N-2)
    BN = BN - V*U(N-2)
    DEN = BN-AN*D(N-1)

    Z[N] = 1.0
    Z(N-1) = -D(N-1)
    NM1 = N-1
    for J in range (2,NM1):
        K = N-J
        Z[K]= -D[K]*Z(K+1)-U[K]+Z[N]
    return
    




    







            







