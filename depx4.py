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
    

    
def TRI3 (M,A,B,C,K,Y1,Y2,Y3,TCOS,D,W1,W2,W3):
    A = np.zeros((1))
    B = np.zeros((1))
    C = np.zeros((1))
    K = np.zeros((4))
    TCOS = np.zeros((1))
    Y1 = np.zeros((1))
    Y2 = np.zeros((1))
    Y3 = np.zeros((1))
    D = np.zeros((1))
    W1 = np.zeros((1))
    W2 = np.zeros((1))
    W3 = np.zeros((1))

    MM1 = M-1
    K1 = K[1]
    K2 = K[2]
    K3 = K[3]
    K4 = K[4]
    F1 = K1+1
    F2 = K2+1
    F3 = K3+1
    F4 = K4+1
    K2K3K4 = K2+K3+K4
    if(K2K3K4 != 0):#--->101
        L1 = F1/F2
        L2 = F1/F3
        L3 = F1/F4
        LINT1 = 1
        LINT2 = 1
        LINT3 = 1
        KINT1 = K1
        KINT2 = KINT1+K2
        KINT3 = KINT2+K3
    
    #101
    for N in range(1,K1):
        X = TCOS[N]
        if(K2K3K4 != 0 ):#---107
            if(N == L1):#---->103
                for I in range(1,M):
                    W1[I] = Y1[I]
                #102
                #CONTINUE
                
            #103
            if(N == L2):#--->105
                for I in range(1,M):
                    W2[I] = Y2[I]
                #CONTINUE
                #104
                
            #105
            if(N == L3):#--->107
                for I in range(1,M):
                    W3[I] = Y3[I]
                #106
                #COTINUE
        #107
        Z = np.double(1.00)/(B[1]-X)
        D[1] = C[1]*Z
        Y1[1] = Y1[1]*Z
        Y2[1] = Y2[1]*Z
        Y3[1] = Y3[1]*Z
        for I in range(2,M):
            Z = np.double(1.00)/(B[I]-X-A[I]*D[I-1])
            D[I] = C[I]*Z
            Y1[I] = (Y1[I]-A[I]*Y1[I-1])*Z
            Y2[I] = (Y2[I]-A[I]*Y2[I-1])*Z
            Y3[I] = (Y3[I]-A[I]*Y3[I-1])*Z
        #108
        for I in range(1,MM1):
            I = M - IP
            Y1[I] = Y1[I]-D[I]*Y1[I+1]
            Y2[I] = Y2[I]-D[I]*Y2[I+1]
            Y3[I] = Y3[I]-D[I]*Y3[I+1]
        #CONTINUE
        #109
        if(K2K3K4 != 0):#--->115
            if(N == L1):#--->111
                I = LINT1+KINT1
                XX = X-TCOS[I]
                for I in range(1,M):
                    Y1[I] = XX*Y1[I]+W1[I]
                #CONTINUE
                #110
                LINT1 = LINT1+1
                L1 = (np.double(LINT1)*F1)/F2

            #111
            if(N == L2):#--->113
                I = LINT2+KINT2
                XX = X-TCOS[I]
                for I in range(1,M):
                    Y2[I] = XX*Y2[I]+W2[I]
                #COTINUE
                #112
                LINT2 = LINT2+1
                L2 = (np.double(LINT2)*F1)/F3
            #113
            if(N == L3):#--->115
                I = LINT3+KINT3
                XX = X-TCOS[I]
                for I in range(1,M):
                    Y3[I] =XX*Y3[I]+W3[I]
                #CONTINUE
                #114
                LINT3 = LINT3+1
                L3 = (np.double(LINT3)*F1)/F4
    #115
    #CONTINUE
    return
                





from cosgenPy import COSGEN
from os import closerange, fdopen


def POISN2(M,N,ISTAG,MIXBND,A,BB,C,Q,IDIMQ,B,B2,B3,W,W2,W3,D,TCOS,P):
    A = np.zeros((1))
    BB = np.zeros((1))
    C = np.zeros((1))
    Q = np.zeros((IDIMQ,1))
    B = np.zeros((1))
    B2 = np.zeros((1))
    B3 = np.zeros((1))
    W = np.zeros((1))
    W2 = np.zeros((1))
    W3 = np.zeros((1))
    D = np.zeros((1))
    TCOS = np.zeros((1))
    K = np.zeros((1))
    P = np.zeros((1))


    FISTAG = 3 - ISTAG
    FNUM = 1.0 / np.double(ISTAG)
    FDEN = 0.5*np.double(ISTAG-1)
    MR = M
    IP = -MR
    IPSTOR = 0
    I2R = 1
    JR = 2
    NR = N
    NLAST = N
    KR = 1
    LR = 0



    if(ISTAG == 1): #-->101
        #101
        for I in range(1,MR+1):
            Q[I][N] = np.double(0.500)*Q[I][N]
        #102
        if (MIXBND == 1):#--->103
            #103
            if(N > 3):#--->155
                while(NR > 1):
                    #104
                    JR = 2*I2R
                    NROD = 0
                    if ((NR/2)*2 == NR):
                        NROD=0
                    if(MIXBND == 1):
                            #105
                            JSTART = 1
                    if (MIXBND == 2):
                        #106
                        JSTART = JR
                        NROD = 1-NROD

                    #107
                    JSTOP = NLAST-JR
                    if (NROD == 0):
                        JSTOP=JSTOP-I2R
                    COSGEN(I2R,1,np.double(0.500),np.double(0.00),TCOS)
                    I2RBY2= I2R/2
                    if(JSTOP < JSTART):#--->108
                        J = JR
                        #--->116
                    #108
                    else:
                        for J in range(JSTART,JSTOP,JR):
                            JP1 = J+I2RBY2
                            JP2 = J+I2R
                            JP3 = JP2+I2RBY2
                            JM1 = J-I2RBY2
                            JM2 = J-I2R
                            JM3 = JM2-I2RBY2
                            if (J == 1):#-->109
                                JM1 = JP1
                                JM2 = JP2
                                JM3 = JP3
                            else:
                                #109
                                if(I2R == 1):#---->111
                                    if(J == 1):
                                        JM2 = JP2
                                    for I in range(1,MR+1):
                                        B[I] = np.double(2.00)*Q[I][J]
                                        Q[I][J] = Q[I][JM2]+Q[I][JP2]
                                    #110
                                    #--->113
                                else:
                                    #111
                                    for I in range(1,MR+1):
                                        FI = Q[I][J]
                                        Q[I][J] = Q[I][J]-Q[I][JM1]-Q[I][JP1]+Q[I][JM2]+Q[I][JP2]
                                        B[I] = FI+Q[I][J]-Q[I][JM3]-Q[I][JP3]
                                    #112
                                    #CONTINUE
                                #113
                                TRIX(I2R,0,MR,A,BB,C,B,TCOS,D,W)
                                for I in range(1,MR+1):
                                    Q[I][J] = Q[I][J]+B[I]
                                #114
                                #COTINUE
                                #115
                                #CONTINUE
                                J = JSTOP+JR
                    #116
                    NLAST = J
                    JM1 = J-I2RBY2
                    JM2 = J-I2R
                    JM3 = JM2-I2RBY2
                    if(NROD != 0):#---->128
                        if(I2R == 1):#--->118
                            for I in range(1,MR+1):
                                B[I] = FISTAG*Q[I][J]
                                Q[I][J] = Q[I][JM2]
                            #117
                            #CONTINUE
                            #---->126
                        else:
                            #118
                            for I in range(1,MR+1):
                                B[I] = Q[I][J]+np.double(0.500)*(Q[I][JM2]-Q[I][JM1]-Q[I][JM3])
                            #119
                            #CONTINUE
                            if(NRODPR == 0):#--->121
                                for I in range(1,MR+1):
                                    II = IP+I
                                    Q[I][J] = Q[I][JM2]+P[II]
                                    #120
                                    #CONTINUE
                                    IP = IP-MR
                                    #--->123
                            else:
                                #121
                                for I in range(1,MR+1):
                                    Q[I][J] = Q[I][J]-Q[I][JM1]+Q[I][JM2]
                                #122
                                #CONTINUE
                            #123
                            if(LR != 0):#---->124
                                COSGEN(LR,1,0,np.double(0.500),FDEN,TCOS(KR+1))
                                #---->126
                            else:
                                #124
                                for I in range(1,MR+1):
                                    B[I] = FISTAG*B[I]
                                #125
                                #CONTINUE
                                #---->126
                        #126
                        COSGEN(KR,1,np.double(0.500),FDEN,TCOS)
                        TRIX(KR,LR,MR,A,BB,C,B,TCOS,D,W)
                        for I in range(1,MR+1):
                            Q[I][J] = Q[I][J]+B[I]
                        #127
                        KR = KR+I2R
                        #---->151
                    else:
                        #128
                        JP1 = J+I2RBY2
                        JP2 = J+I2R
                        if(I2R == 1):#---->135
                            for I in range(1,MR+1):
                                B[I] = Q[I][J]
                            TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                            IP = 0
                            IPSTOR = MR
                            if(ISTAG == 1):#--->133
                                #133
                                for I in range(1,MR+1):
                                    P[I] = B[I]
                                    Q[I][J] = Q[I][JM2]+np.double(2.00)*Q[I][JP2]+np.double(3.00)*B[I]
                                #134
                                #COTINUE
                                #---->150
                            if(ISTAG == 2):#--->130
                                #130
                                for I in range(1,MR+1):
                                    P[I] = B[I]
                                    B[I] = B[I]+Q[I][N]
                                #131
                                #CONTINUE
                                TCOS[1] =  np.double(1.00)
                                TCOS[2] = np.double(0.00)
                                TRIX(1,1,MR,A,BB,C,B,TCOS,D,W)
                                for I in range(1,MR+1):
                                    Q[I][J] = Q[I][JM2]+P[I]+B[I]
                                #132
                                #CONTINUE
                                #---->150
                        else:
                            #135
                            for I in range(1,MR+1):
                                B[I] = Q[I][J]+np.double(0.500)*(Q[I][JM2]-Q[I][JM1]-Q[I][JM3])
                            #136
                            if (NRODPR == 0):#--->138
                                for I in range(1,MR+1):
                                    II = IP+I
                                    B[I] = B[I]+P[II]
                                #137
                                #CONTINUE
                                #--->140
                            else:
                                #138
                                for I in range(1,MR+1):
                                    B[I] = B[I]+Q[I][JP2]-Q[I][JP1]
                                #139
                                #CONTINUE
                            #140
                            TRIX(I2R,0,MR,A,BB,C,B,TCOS,D,W)
                            IP = IP+MR
                            IPSTOR = MAX0(IPSTOR,IP+MR)
                            for I in range(1,MR+1):
                                II = IP+I
                                P[II] = B[I]+np.double(0.500)*(Q[I][J]-Q[I][JM1]-Q[I][JP1])
                                B[I] =P[II]+Q[I][JP2]
                            #CONTINUE
                            #141
                            if (LR != 0):#--->142
                                COSGEN(LR,1,np.double(0.500),FDEN,TCOS(I2R+1))
                                MERGE(TCOS,0,I2R,I2R,LR,KR)
                                #---->144
                            else:
                                #142
                                for I in range(1,MR+1):
                                    II = KR+I
                                    TCOS[II] = TCOS[I]
                                #COTINUE
                                #143
                                #--->144
                            #144
                            COSGEN(KR,1,np.double(0.500),FDEN,TCOS)
                            if(LR == 0):#--->145
                                if(ISTAG == 1):#--->146
                                    #146
                                    for I in range(1,MR+1):
                                        B[I] = FISTAG*B[I]
                                    #CONTINUE
                                    #147
                                    #---->148
                                if(ISTAG == 2): #----->145
                                    #145
                                    TRIX(KR,KR,MR,A,BB,C,B,TCOS,D,W)
                                    #---->148
                            else:
                                #145
                                TRIX(KR,KR,MR,A,BB,C,B,TCOS,D,W)
                                #--->148

                            #148
                            for I in range(1,MR+1):
                                II = IP+I
                                Q[I][J] = Q[I][JM2]+P[II]+B[I]
                            #149
                            #COTINUE
                        #150
                        LR = KR
                        KR = KR+JR
                        #CONTINUE
                    #151
                    if(MIXBND == 1):#--->152
                        #152
                        NR = (NLAST-1)/JR+1
                        if(NR <= 3):#---155
                            #155
                            JM1 = J-I2R
                            JP1 = J+I2R
                            JM2 = NLAST-I2R
                            if(NR == 2):#---->184
                                #184
                                if(N == 2):#--->188
                                    for I in range(1,MR+1):
                                        B[I] = Q[I][1]
                                    #185
                                    #CONTINUE
                                    TCOS[1] = 0
                                    TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                                    for I in range(1,MR+1):
                                        Q[I][1] = B[I]
                                        B[I] = np.double(2.00)*(Q[I][2]+B[I])*FISTAG
                                    #186
                                    #CONTINUE
                                    TCOS[1] = -FISTAG
                                    TCOS[2] = np.double(2.00)
                                    TRIX(2,0,MR,A,BB,C,B,TCOS,D,W)
                                    for I in range(1,MR+1):
                                        Q[I][1] = Q[I][1]+B[I]
                                    #187
                                    #CONTINUE
                                    JR = 1
                                    I2R = 0
                                    #--->194
                                else:
                                    #188
                                    for I in range(1,MR+1):
                                        II = IP+1
                                        B3[I] = np.double(0.00)
                                        B[I] = Q[I][1]+np.double(2.00)*P[II]
                                        Q[I][1] = np.double(0.500)*Q[I][1]-Q[I][JM1]
                                        B2[I] = np.double(2.00)*(Q[I][1]+Q[I][NLAST])
                                    #189
                                    K1 = KR+JR-1
                                    TCOS[K1+1] = np.double(-2.00)
                                    K4 = K1+3-ISTAG
                                    COSGEN(KR+ISTAG-2,1,np.double(0.00),FNUM,TCOS(K4))
                                    K4 = K1+KR+1
                                    COSGEN(JR-1,1,np.double(0.00),np.double(1.00),TCOS(K4))
                                    MERGE(TCOS,K1,KR,K1+KR,JR-1,0)
                                    COSGEN(KR,1,np.double(0.500),FDEN,TCOS(K1+1))
                                    K2 = KR
                                    K4 = K1+K2+1
                                    COSGEN(LR,1,np.double(0.500),FDEN,TCOS(K4))
                                    K3 = LR
                                    K4 = 0
                                    TRI3(MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
                                    for I in range(1,MR+1):
                                        B[I] = B[I]+B2[I]
                                    #COTINUE
                                    #190
                                    TCOS[1] = np.double(2.00)
                                    TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                                    for I in range(1,MR+1):
                                        Q[I][1] = Q[I][NLAST]
                                    #191
                                    # --->194
                                #194
                                J = NLAST - JR
                                for I in range(1,MR+1):
                                    B[I] = Q[I][NLAST]+Q[I][J]
                                #195
                                #CONTINUE
                                #--->196
                            if(LR != 0):#--->170
                                #170
                                for I in range(1,MR+1):
                                    B[I] = np.double(0.500)*Q[I][1]-Q[I][JM1]+Q[I][J]
                                #171
                                if(NROD == 0):#---->173
                                    for I in range(1,MR+1):
                                        II = IP+I
                                        B[I] = B[I]+P[II]
                                        #172
                                        #CONTINUE
                                        #--->175
                                else:
                                    #173
                                    for I in range(1,MR+1):
                                        B[I] = B[I]+Q[I][NLAST]-Q[I][JM2]
                                    #174
                                    #CONTINUE
                                    #--->175
                                #175
                                for I in range(1,MR+1):
                                    T = np.double(0.500)*(Q[I][J]-Q[I][JM1]-Q[I][JP1])
                                    Q[I][J] = T
                                    B2[I] = Q[I][NLAST]+T
                                    B3[I] = Q[I][1]+np.double(2.00)*T
                                #176
                                #177
                                K1 = KR+2*JR-1
                                K2 = KR+JR
                                TCOS[K1+1] = np.double(-2.00)
                                K4 = K1+3-ISTAG
                                COSGEN(K2+ISTAG-2,1,np.double(0.00),FNUM,TCOS(K4))
                                K4 = K1+K2+1
                                COSGEN(JR-1,1,np.double(0.00),np.double(1.00),TCOS(K4))
                                MERGE(TCOS,K1,K2,K1+K2,JR-1,0)
                                K3 = K1+K2+LR
                                COSGEN(JR,1,np.double(0.500),np.double(0.00),TCOS(K3+1))
                                K4 = K3+JR+1
                                COSGEN(KR,1,np.double(0.500),FDEN,TCOS(K4))
                                MERGE(TCOS,K3,JR,K3+JR,KR,K1)
                                if(LR != 0):#---->178
                                    COSGEN(LR,1,np.double(0.500),FDEN,TCOS(K4))
                                    MERGE(TCOS,K3,JR,K3+JR,LR,K3-LR)
                                    COSGEN(KR,1,np.double(0.500),FDEN,TCOS(K4))
                                #178
                                K3 = KR
                                K4 = KR
                                TRI3(MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
                                for I in range(1,MR+1):
                                    B[I] = B[I]+B2[I]+B3[I]
                                #179
                                TCOS[1] = np.double(2.00)
                                TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                                for I in range(1,MR+1):
                                    Q[I][J] = Q[I][J] +B[I]
                                    B[I] = Q[I][1]+np.double(2.00)*Q[I][J]
                                #180
                                COSGEN(JR,1,np.double(0.500),np.double(0.00),TCOS)
                                TRIX(JR,0,MR,A,BB,C,B,TCOS,D,W)
                                if(JR == 1):#--->182
                                    for I in range(1,MR+1):
                                        Q[I][1] = B[I]
                                    #181
                                    #CONTINUE
                                    #---->194
                                else:
                                    #182
                                    for I in range(1,MR+1):
                                        Q[I][1] = np.double(0.500)*Q[I][JM1]+B[I]
                                    #183
                                    #CONTINUE
                                    #--->194
                                #194
                                J = NLAST - JR
                                for I in range(1,MR+1):
                                    B[I] = Q[I][NLAST]+Q[I][J]
                                #195
                                #CONTINUE
                                #--->196
                            if(N != 3):#--->161
                                #161
                                if(ISTAG == 1):#--->162
                                    #162
                                    for I in range(1,MR+1):
                                        B[I] = Q[I][J]+np.double(0.500)*Q[I][1]-Q[I][JM1]+Q[I][NLAST]-Q[I][JM2]
                                    #163
                                    COSGEN(JR,1,np.double(0.500),np.double(0.00),TCOS)
                                    TRIX(JR,0,MR,A,BB,C,B,TCOS,D,W)
                                    for I in range(1,MR+1):
                                        Q[I][J] = np.double(0.500)*(Q[I][J]-Q[I][JM1]-Q[I][JP1])+B[I]
                                        B[I] = Q[I][1]+np.double(2.00)*Q[I][NLAST]+np.double(4.00)*Q[I][J]
                                    #164
                                    JR2 = 2*JR
                                    COSGEN(JR,1,np.double(0.00),np.double(0.00),TCOS)
                                    for I in range(1,JR):
                                        I1 = JR+1
                                        I2 = JR+1-I
                                        TCOS[I1] = -TCOS[I2]
                                    #165
                                    #CONTINUE
                                    TRIX(JR2,0,MR,A,BB,C,B,TCOS,D,W)
                                    for I in range(1,MR+1):
                                        Q[I][J] = Q[I][J]+B[I]
                                        B[I] = Q[I][1]+np.double(2.00)*Q[I][J]
                                    #166
                                    COSGEN(KR,1,np.double(0.500),np.double(0.00),TCOS)
                                    TRIX(KR,0,MR,A,BB,C,B,TCOS,D,W)
                                    for I in range(1,MR+1):
                                        Q[I][1] =np.double(0.500)*Q[I][1]-Q[I][JM1]+B[I]
                                    #COTINUE
                                    #167
                                    #---->194
                                if(ISTAG == 2):#--->170
                                    #170
                                    for I in range(1,MR+1):
                                        B[I] = np.double(0.500)*Q[I][1]-Q[I][JM1]+Q[I][J]
                                    #171
                                    if(NROD == 0):#---->173
                                        for I in range(1,MR+1):
                                            II = IP+I
                                            B[I] = B[I]+P[II]
                                            #172
                                            #CONTINUE
                                            #--->175
                                    else:
                                        #173
                                        for I in range(1,MR+1):
                                            B[I] = B[I]+Q[I][NLAST]-Q[I][JM2]
                                        #174
                                        #CONTINUE
                                        #--->175
                                    #175
                                    for I in range(1,MR+1):
                                        T = np.double(0.500)*(Q[I][J]-Q[I][JM1]-Q[I][JP1])
                                        Q[I][J] = T
                                        B2[I] = Q[I][NLAST]+T
                                        B3[I] = Q[I][1]+np.double(2.00)*T
                                    #176
                                    #177
                                    K1 = KR+2*JR-1
                                    K2 = KR+JR
                                    TCOS[K1+1] = np.double(-2.00)
                                    K4 = K1+3-ISTAG
                                    COSGEN(K2+ISTAG-2,1,np.double(0.00),FNUM,TCOS(K4))
                                    K4 = K1+K2+1
                                    COSGEN(JR-1,1,np.double(0.00),np.double(1.00),TCOS(K4))
                                    MERGE(TCOS,K1,K2,K1+K2,JR-1,0)
                                    K3 = K1+K2+LR
                                    COSGEN(JR,1,np.double(0.500),np.double(0.00),TCOS(K3+1))
                                    K4 = K3+JR+1
                                    COSGEN(KR,1,np.double(0.500),FDEN,TCOS(K4))
                                    MERGE(TCOS,K3,JR,K3+JR,KR,K1)
                                    if(LR != 0):#---->178
                                        COSGEN(LR,1,np.double(0.500),FDEN,TCOS(K4))
                                        MERGE(TCOS,K3,JR,K3+JR,LR,K3-LR)
                                        COSGEN(KR,1,np.double(0.500),FDEN,TCOS(K4))
                                    #178
                                    K3 = KR
                                    K4 = KR
                                    TRI3(MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
                                    for I in range(1,MR+1):
                                        B[I] = B[I]+B2[I]+B3[I]
                                    #179
                                    TCOS[1] = np.double(2.00)
                                    TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                                    for I in range(1,MR+1):
                                        Q[I][J] = Q[I][J] +B[I]
                                        B[I] = Q[I][1]+np.double(2.00)*Q[I][J]
                                    #180
                                    COSGEN(JR,1,np.double(0.500),np.double(0.00),TCOS)
                                    TRIX(JR,0,MR,A,BB,C,B,TCOS,D,W)
                                    if(JR == 1):#--->182
                                        for I in range(1,MR+1):
                                            Q[I][1] = B[I]
                                        #181
                                        #CONTINUE
                                        #---->194
                                    else:
                                        #182
                                        for I in range(1,MR+1):
                                            Q[I][1] = np.double(0.500)*Q[I][JM1]+B[I]
                                        #183
                                        #CONTINUE
                                        #--->194
                                #194
                                J = NLAST - JR
                                for I in range(1,MR+1):
                                    B[I] = Q[I][NLAST]+Q[I][J]
                                #195
                                #CONTINUE
                                #--->196
                            else:
                                #156
                                for I in range(1,MR+1):
                                    B[I] = Q[I][2]
                                #CONTINUE
                                #157
                                TCOS[1] = np.double(0.00)
                                TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                                for I in range(1,MR+1):
                                    Q[I][2] = B[I]
                                    B[I] = np.double(4.00)*B[I][1]+np.double(2.00)*Q[I][3]
                                #158
                                #CONTINUE

                                TCOS[1] = np.double(-2.00)
                                TCOS[2] = np.double(2.00)
                                I1 = 0
                                I2 = 0
                                TRIX(I1,I2,MR,A,BB,C,B,TCOS,D,W)
                                for I in range(1,MR+1):
                                    Q[I][2] = Q[I][2] + B[I]
                                    B[I] = Q[I][1] + np.double(2.00)*Q[I][2]
                                #159
                                TCOS[1] = np.double(0.00)
                                TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                                for I in range(1,MR+1):
                                    Q[I][1] = B[I]
                                #160
                                JR = 1
                                I2R = 0
                                #--->194
                                #194
                                J = NLAST - JR
                                for I in range(1,MR+1):
                                    B[I] = Q[I][NLAST]+Q[I][J]
                                #195
                                #CONTINUE
                                #--->196
                            break
                        else:
                            I2R = JR
                            NRODPR = NROD
                    if(MIXBND == 2):#---->153
                        #153
                        NR = NLAST/JR
                        if(NR <= 1):#--->192
                            #192
                            for I in range(1,MR+1):
                                B[I] = Q[I][NLAST]
                            #--->196
                            break
                        else:
                            I2R = JR
                            NRODPR = NROD
            else:
                #155
                JM1 = J-I2R
                JP1 = J+I2R
                JM2 = NLAST-I2R
                if(NR == 2):#---->184
                    #184
                    if(N == 2):#--->188
                        for I in range(1,MR+1):
                            B[I] = Q[I][1]
                        #185
                        #CONTINUE
                        TCOS[1] = 0
                        TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                        for I in range(1,MR+1):
                            Q[I][1] = B[I]
                            B[I] = np.double(2.00)*(Q[I][2]+B[I])*FISTAG
                        #186
                        #CONTINUE
                        TCOS[1] = -FISTAG
                        TCOS[2] = np.double(2.00)
                        TRIX(2,0,MR,A,BB,C,B,TCOS,D,W)
                        for I in range(1,MR+1):
                            Q[I][1] = Q[I][1]+B[I]
                        #187
                        #CONTINUE
                        JR = 1
                        I2R = 0
                        #--->194
                    else:
                        #188
                        for I in range(1,MR+1):
                            II = IP+1
                            B3[I] = np.double(0.00)
                            B[I] = Q[I][1]+np.double(2.00)*P[II]
                            Q[I][1] = np.double(0.500)*Q[I][1]-Q[I][JM1]
                            B2[I] = np.double(2.00)*(Q[I][1]+Q[I][NLAST])
                        #189
                        K1 = KR+JR-1
                        TCOS[K1+1] = np.double(-2.00)
                        K4 = K1+3-ISTAG
                        COSGEN(KR+ISTAG-2,1,np.double(0.00),FNUM,TCOS(K4))
                        K4 = K1+KR+1
                        COSGEN(JR-1,1,np.double(0.00),np.double(1.00),TCOS(K4))
                        MERGE(TCOS,K1,KR,K1+KR,JR-1,0)
                        COSGEN(KR,1,np.double(0.500),FDEN,TCOS(K1+1))
                        K2 = KR
                        K4 = K1+K2+1
                        COSGEN(LR,1,np.double(0.500),FDEN,TCOS(K4))
                        K3 = LR
                        K4 = 0
                        TRI3(MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
                        for I in range(1,MR+1):
                            B[I] = B[I]+B2[I]
                        #COTINUE
                        #190
                        TCOS[1] = np.double(2.00)
                        TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                        for I in range(1,MR+1):
                            Q[I][1] = Q[I][NLAST]
                        #191
                        # --->194
                    #194
                    J = NLAST - JR
                    for I in range(1,MR+1):
                        B[I] = Q[I][NLAST]+Q[I][J]
                    #195
                    #CONTINUE
                    #--->196
                if(LR != 0):#--->170
                    #170
                    for I in range(1,MR+1):
                        B[I] = np.double(0.500)*Q[I][1]-Q[I][JM1]+Q[I][J]
                    #171
                    if(NROD == 0):#---->173
                        for I in range(1,MR+1):
                            II = IP+I
                            B[I] = B[I]+P[II]
                            #172
                            #CONTINUE
                            #--->175
                    else:
                        #173
                        for I in range(1,MR+1):
                            B[I] = B[I]+Q[I][NLAST]-Q[I][JM2]
                        #174
                        #CONTINUE
                        #--->175
                    #175
                    for I in range(1,MR+1):
                        T = np.double(0.500)*(Q[I][J]-Q[I][JM1]-Q[I][JP1])
                        Q[I][J] = T
                        B2[I] = Q[I][NLAST]+T
                        B3[I] = Q[I][1]+np.double(2.00)*T
                    #176
                    #177
                    K1 = KR+2*JR-1
                    K2 = KR+JR
                    TCOS[K1+1] = np.double(-2.00)
                    K4 = K1+3-ISTAG
                    COSGEN(K2+ISTAG-2,1,np.double(0.00),FNUM,TCOS(K4))
                    K4 = K1+K2+1
                    COSGEN(JR-1,1,np.double(0.00),np.double(1.00),TCOS(K4))
                    MERGE(TCOS,K1,K2,K1+K2,JR-1,0)
                    K3 = K1+K2+LR
                    COSGEN(JR,1,np.double(0.500),np.double(0.00),TCOS(K3+1))
                    K4 = K3+JR+1
                    COSGEN(KR,1,np.double(0.500),FDEN,TCOS(K4))
                    MERGE(TCOS,K3,JR,K3+JR,KR,K1)
                    if(LR != 0):#---->178
                        COSGEN(LR,1,np.double(0.500),FDEN,TCOS(K4))
                        MERGE(TCOS,K3,JR,K3+JR,LR,K3-LR)
                        COSGEN(KR,1,np.double(0.500),FDEN,TCOS(K4))
                    #178
                    K3 = KR
                    K4 = KR
                    TRI3(MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
                    for I in range(1,MR+1):
                        B[I] = B[I]+B2[I]+B3[I]
                    #179
                    TCOS[1] = np.double(2.00)
                    TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                    for I in range(1,MR+1):
                        Q[I][J] = Q[I][J] +B[I]
                        B[I] = Q[I][1]+np.double(2.00)*Q[I][J]
                    #180
                    COSGEN(JR,1,np.double(0.500),np.double(0.00),TCOS)
                    TRIX(JR,0,MR,A,BB,C,B,TCOS,D,W)
                    if(JR == 1):#--->182
                        for I in range(1,MR+1):
                            Q[I][1] = B[I]
                        #181
                        #CONTINUE
                        #---->194
                    else:
                        #182
                        for I in range(1,MR+1):
                            Q[I][1] = np.double(0.500)*Q[I][JM1]+B[I]
                        #183
                        #CONTINUE
                        #--->194
                    #194
                    J = NLAST - JR
                    for I in range(1,MR+1):
                        B[I] = Q[I][NLAST]+Q[I][J]
                    #195
                    #CONTINUE
                    #--->196
                if(N != 3):#--->161
                    #161
                    if(ISTAG == 1):#--->162
                        #162
                        for I in range(1,MR+1):
                            B[I] = Q[I][J]+np.double(0.500)*Q[I][1]-Q[I][JM1]+Q[I][NLAST]-Q[I][JM2]
                        #163
                        COSGEN(JR,1,np.double(0.500),np.double(0.00),TCOS)
                        TRIX(JR,0,MR,A,BB,C,B,TCOS,D,W)
                        for I in range(1,MR+1):
                            Q[I][J] = np.double(0.500)*(Q[I][J]-Q[I][JM1]-Q[I][JP1])+B[I]
                            B[I] = Q[I][1]+np.double(2.00)*Q[I][NLAST]+np.double(4.00)*Q[I][J]
                        #164
                        JR2 = 2*JR
                        COSGEN(JR,1,np.double(0.00),np.double(0.00),TCOS)
                        for I in range(1,JR):
                            I1 = JR+1
                            I2 = JR+1-I
                            TCOS[I1] = -TCOS[I2]
                        #165
                        #CONTINUE
                        TRIX(JR2,0,MR,A,BB,C,B,TCOS,D,W)
                        for I in range(1,MR+1):
                            Q[I][J] = Q[I][J]+B[I]
                            B[I] = Q[I][1]+np.double(2.00)*Q[I][J]
                        #166
                        COSGEN(KR,1,np.double(0.500),np.double(0.00),TCOS)
                        TRIX(KR,0,MR,A,BB,C,B,TCOS,D,W)
                        for I in range(1,MR+1):
                            Q[I][1] =np.double(0.500)*Q[I][1]-Q[I][JM1]+B[I]
                        #COTINUE
                        #167
                        #---->194
                    if(ISTAG == 2):#--->170
                        #170
                        for I in range(1,MR+1):
                            B[I] = np.double(0.500)*Q[I][1]-Q[I][JM1]+Q[I][J]
                        #171
                        if(NROD == 0):#---->173
                            for I in range(1,MR+1):
                                II = IP+I
                                B[I] = B[I]+P[II]
                                #172
                                #CONTINUE
                                #--->175
                        else:
                            #173
                            for I in range(1,MR+1):
                                B[I] = B[I]+Q[I][NLAST]-Q[I][JM2]
                            #174
                            #CONTINUE
                            #--->175
                        #175
                        for I in range(1,MR+1):
                            T = np.double(0.500)*(Q[I][J]-Q[I][JM1]-Q[I][JP1])
                            Q[I][J] = T
                            B2[I] = Q[I][NLAST]+T
                            B3[I] = Q[I][1]+np.double(2.00)*T
                        #176
                        #177
                        K1 = KR+2*JR-1
                        K2 = KR+JR
                        TCOS[K1+1] = np.double(-2.00)
                        K4 = K1+3-ISTAG
                        COSGEN(K2+ISTAG-2,1,np.double(0.00),FNUM,TCOS(K4))
                        K4 = K1+K2+1
                        COSGEN(JR-1,1,np.double(0.00),np.double(1.00),TCOS(K4))
                        MERGE(TCOS,K1,K2,K1+K2,JR-1,0)
                        K3 = K1+K2+LR
                        COSGEN(JR,1,np.double(0.500),np.double(0.00),TCOS(K3+1))
                        K4 = K3+JR+1
                        COSGEN(KR,1,np.double(0.500),FDEN,TCOS(K4))
                        MERGE(TCOS,K3,JR,K3+JR,KR,K1)
                        if(LR != 0):#---->178
                            COSGEN(LR,1,np.double(0.500),FDEN,TCOS(K4))
                            MERGE(TCOS,K3,JR,K3+JR,LR,K3-LR)
                            COSGEN(KR,1,np.double(0.500),FDEN,TCOS(K4))
                        #178
                        K3 = KR
                        K4 = KR
                        TRI3(MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
                        for I in range(1,MR+1):
                            B[I] = B[I]+B2[I]+B3[I]
                        #179
                        TCOS[1] = np.double(2.00)
                        TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                        for I in range(1,MR+1):
                            Q[I][J] = Q[I][J] +B[I]
                            B[I] = Q[I][1]+np.double(2.00)*Q[I][J]
                        #180
                        COSGEN(JR,1,np.double(0.500),np.double(0.00),TCOS)
                        TRIX(JR,0,MR,A,BB,C,B,TCOS,D,W)
                        if(JR == 1):#--->182
                            for I in range(1,MR+1):
                                Q[I][1] = B[I]
                            #181
                            #CONTINUE
                            #---->194
                        else:
                            #182
                            for I in range(1,MR+1):
                                Q[I][1] = np.double(0.500)*Q[I][JM1]+B[I]
                            #183
                            #CONTINUE
                            #--->194
                    #194
                    J = NLAST - JR
                    for I in range(1,MR+1):
                        B[I] = Q[I][NLAST]+Q[I][J]
                    #195
                    #CONTINUE
                    #--->196

                else:
                    #156
                    for I in range(1,MR+1):
                        B[I] = Q[I][2]
                    #CONTINUE
                    #157
                    TCOS[1] = np.double(0.00)
                    TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                    for I in range(1,MR+1):
                        Q[I][2] = B[I]
                        B[I] = np.double(4.00)*B[I][1]+np.double(2.00)*Q[I][3]
                    #158
                    #CONTINUE

                    TCOS[1] = np.double(-2.00)
                    TCOS[2] = np.double(2.00)
                    I1 = 0
                    I2 = 0
                    TRIX(I1,I2,MR,A,BB,C,B,TCOS,D,W)
                    for I in range(1,MR+1):
                        Q[I][2] = Q[I][2] + B[I]
                        B[I] = Q[I][1] + np.double(2.00)*Q[I][2]
                    #159
                    TCOS[1] = np.double(0.00)
                    TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                    for I in range(1,MR+1):
                        Q[I][1] = B[I]
                    #160
                    JR = 1
                    I2R = 0
                    #--->194
                    #194
                    J = NLAST - JR
                    for I in range(1,MR+1):
                        B[I] = Q[I][NLAST]+Q[I][J]
                    #195
                    #CONTINUE
                    #--->196
            #196
            JM2 = NLAST  - I2R
            if(JR == 1):#--->198
                for I in range(1,MR+1):
                    Q[I][NLAST] = np.double(0.000)
                #197
                #CONTINUE
                #--->202
            else:
                #198
                if (NROD == 0):#--->200
                    for I in range(1,MR+1):
                        II = IP+I
                        Q[I][NLAST] = P[II]
                    #199
                    IP = IP+MR
                    #---->202
                else:
                    #200
                    for I in range(1,MR+1):
                        Q[I][NLAST] = Q[I][NLAST]-Q[I][JM2]
                    #201
                    #CONTINUE
                    #--->202
            #202
            COSGEN (KR,1,np.double(0.500),FDEN,TCOS)
            COSGEN (LR,1,np.double(0.500),FDEN,TCOS(KR+1))
            if (LR == 0):#------>204
                for I in range(1,MR+1):
                    B[I] = FISTAG*B[I]
            #203
            #CONTINUE
            #204
            TRIX (KR,LR,MR,A,BB,C,B,TCOS,D,W)
            for I in range(1,MR+1):
                Q[I][NLAST] = Q[I][NLAST]+B
            #CONTIUE
            #205
            NLASTP = NLAST
            #206
            JSTEP = JR
            JR = I2R
            I2R = I2R/2
            while(JR > 0):
                #
                if (JR != 0):#--->222
                    if(MIXBND == 1):#---->207
                        #207
                        JSTART = 1+JR
                        #---->209
                    if(MIXBND == 2):#---->208
                        #208
                        JSTART = JR
                    #209
                    KR = KR-JR
                    if (NLAST+JR < N):#----->210
                        KR = KR-JR
                        NLAST = NLAST+JR
                        JSTOP = NLAST-JSTEP
                        #---->211
                    else:
                        #210
                        JSTOP = NLAST-JR
                        #--->211
                    #211
                    LR = KR-JR
                    COSGEN (JR,1,np.double(0.500),np.double(0.000),TCOS)
                    for J in range(JSTART,JSTOP,JSTEP):
                        JM2 = J-JR
                        JP2 = J+JR
                        if(J == JR):#--->213
                            for I in range(1,MR+1):
                                B[I] = Q[I][J]+Q[I][JP2]
                            #212
                            #CONTINUE
                            #215
                        else:
                            #213
                            for I in range(1,MR+1):
                                B[I] = Q[I][J]+Q[I][JM2]+Q[I][JP2]
                            #214
                            #CONTINUE
                        #215
                        if(JR != 1):#---->217
                            for I in range(1,MR+1):
                                Q[I][J] = np.double(0.00)
                            #216
                            #CONTINUE
                            #---->219
                        else:
                            #217
                            JM1 = J-I2R
                            JP1 = J+I2R
                            for I in range(1,MR+1):
                                Q[I][J] = np.double(0.500)*(Q[I][J]-Q[I][JM1]-Q[I][JP1])
                            #218
                            #CONTINUE
                        #219
                        TRIX (JR,0,MR,A,BB,C,B,TCOS,D,W)
                        for I in range(1,MR+1):
                            Q[I][J] = Q[I][J]+B[I]
                            #220
                            #CONTINUE
                    #221
                    #CONTINUE
                    NROD = 1
                    if(NLAST+I2R > N):
                        NROD = 0
                    if(NLASTP != NLAST):
                        #--->194
                        #194
                        J = NLAST - JR
                        for I in range(1,MR+1):
                            B[I] = Q[I][NLAST]+Q[I][J]
                        #195
                        #CONTINUE
                        #--->196
                        #196
                        JM2 = NLAST  - I2R
                        if(JR == 1):#--->198
                            for I in range(1,MR+1):
                                Q[I][NLAST] = np.double(0.000)
                            #197
                            #CONTINUE
                            #--->202
                        else:
                            #198
                            if (NROD == 0):#--->200
                                for I in range(1,MR+1):
                                    II = IP+I
                                    Q[I][NLAST] = P[II]
                                #199
                                IP = IP+MR
                                #---->202
                            else:
                                #200
                                for I in range(1,MR+1):
                                    Q[I][NLAST] = Q[I][NLAST]-Q[I][JM2]
                                #201
                                #CONTINUE
                                #--->202
                        #202
                        COSGEN (KR,1,np.double(0.500),FDEN,TCOS)
                        COSGEN (LR,1,np.double(0.500),FDEN,TCOS(KR+1))
                        if (LR == 0):#------>204
                            for I in range(1,MR+1):
                                B[I] = FISTAG*B[I]
                        #203
                        #CONTINUE
                        #204
                        TRIX (KR,LR,MR,A,BB,C,B,TCOS,D,W)
                        for I in range(1,MR+1):
                            Q[I][NLAST] = Q[I][NLAST]+B
                        #CONTIUE
                        #205
                        NLASTP = NLAST
                        #206
                        JSTEP = JR
                        JR = I2R
                        I2R = I2R/2
                    else:
                        #206
                        JSTEP = JR
                        JR = I2R
                        I2R = I2R/2
                else:
                    break
            #222
            W[1] = IPSTOR
            return

                        
                        



                            

                





                        



                            
        if(MIXBND == 2):#--->104
            while(NR > 1):

                #104
                JR = 2*I2R
                NROD = 1
                if((NR/2)*2 == NR ):
                    NROD = 0
                if(MIXBND == 1):#--->105
                    #105
                    JSTART = 1
                    #--->107
                if(MIXBND == 2):#--->106
                    #106
                    JSTART = JR
                    NROD = 1-NROD
                    #---->107
                #107
                JSTOP = NLAST-JR
                if(NROD == 0):
                    JSTOP = JSTOP-I2R
                COSGEN(I2R,1,np.double(0.500),np.double(0.00),TCOS)
                I2RBY2= I2R/2
                if(JSTOP < JSTART):#--->108
                    J = JR
                    #--->116
                #108
                else:
                    for J in range(JSTART,JSTOP,JR):
                        JP1 = J+I2RBY2
                        JP2 = J+I2R
                        JP3 = JP2+I2RBY2
                        JM1 = J-I2RBY2
                        JM2 = J-I2R
                        JM3 = JM2-I2RBY2
                        if (J == 1):#-->109
                            JM1 = JP1
                            JM2 = JP2
                            JM3 = JP3
                        else:
                            #109
                            if(I2R == 1):#---->111
                                if(J == 1):
                                    JM2 = JP2
                                for I in range(1,MR+1):
                                    B[I] = np.double(2.00)*Q[I][J]
                                    Q[I][J] = Q[I][JM2]+Q[I][JP2]
                                #110
                                #--->113
                            else:
                                #111
                                for I in range(1,MR+1):
                                    FI = Q[I][J]
                                    Q[I][J] = Q[I][J]-Q[I][JM1]-Q[I][JP1]+Q[I][JM2]+Q[I][JP2]
                                    B[I] = FI+Q[I][J]-Q[I][JM3]-Q[I][JP3]
                                #112
                                #CONTINUE
                            #113
                            TRIX(I2R,0,MR,A,BB,C,B,TCOS,D,W)
                            for I in range(1,MR+1):
                                Q[I][J] = Q[I][J]+B[I]
                            #114
                            #COTINUE
                            #115
                            #CONTINUE
                            J = JSTOP+JR
                #116
                NLAST = J
                JM1 = J-I2RBY2
                JM2 = J-I2R
                JM3 = JM2-I2RBY2
                if(NROD != 0):#---->128
                    if(I2R == 1):#--->118
                        for I in range(1,MR+1):
                            B[I] = FISTAG*Q[I][J]
                            Q[I][J] = Q[I][JM2]
                        #117
                        #CONTINUE
                        #---->126
                    else:
                        #118
                        for I in range(1,MR+1):
                            B[I] = Q[I][J]+np.double(0.500)*(Q[I][JM2]-Q[I][JM1]-Q[I][JM3])
                        #119
                        #CONTINUE
                        if(NRODPR == 0):#--->121
                            for I in range(1,MR+1):
                                II = IP+I
                                Q[I][J] = Q[I][JM2]+P[II]
                                #120
                                #CONTINUE
                                IP = IP-MR
                                #--->123
                        else:
                            #121
                            for I in range(1,MR+1):
                                Q[I][J] = Q[I][J]-Q[I][JM1]+Q[I][JM2]
                            #122
                            #CONTINUE
                        #123
                        if(LR != 0):#---->124
                            COSGEN(LR,1,0,np.double(0.500),FDEN,TCOS(KR+1))
                            #---->126
                        else:
                            #124
                            for I in range(1,MR+1):
                                B[I] = FISTAG*B[I]
                            #125
                            #CONTINUE
                            #---->126
                    #126
                    COSGEN(KR,1,np.double(0.500),FDEN,TCOS)
                    TRIX(KR,LR,MR,A,BB,C,B,TCOS,D,W)
                    for I in range(1,MR+1):
                        Q[I][J] = Q[I][J]+B[I]
                    #127
                    KR = KR+I2R
                    #---->151
                else:
                    #128
                    JP1 = J+I2RBY2
                    JP2 = J+I2R
                    if(I2R == 1):#---->135
                        for I in range(1,MR+1):
                            B[I] = Q[I][J]
                        TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                        IP = 0
                        IPSTOR = MR
                        if(ISTAG == 1):#--->133
                            #133
                            for I in range(1,MR+1):
                                P[I] = B[I]
                                Q[I][J] = Q[I][JM2]+np.double(2.00)*Q[I][JP2]+np.double(3.00)*B[I]
                        #134
                        #COTINUE
                        #---->150
                        if(ISTAG == 2):#--->130
                            #130
                            for I in range(1,MR+1):
                                P[I] = B[I]
                                B[I] = B[I]+Q[I][N]
                            #131
                            #CONTINUE
                            TCOS[1] =  np.double(1.00)
                            TCOS[2] = np.double(0.00)
                            TRIX(1,1,MR,A,BB,C,B,TCOS,D,W)
                            for I in range(1,MR+1):
                                Q[I][J] = Q[I][JM2]+P[I]+B[I]
                            #132
                            #CONTINUE
                            #---->150
                    else:
                        #135
                        for I in range(1,MR+1):
                            B[I] = Q[I][J]+np.double(0.500)*(Q[I][JM2]-Q[I][JM1]-Q[I][JM3])
                         #136
                        if (NRODPR == 0):#--->138
                            for I in range(1,MR+1):
                                II = IP+I
                                B[I] = B[I]+P[II]
                            #137
                            #CONTINUE
                            #--->140
                        else:
                            #138
                            for I in range(1,MR+1):
                                B[I] = B[I]+Q[I][JP2]-Q[I][JP1]
                            #139
                            #CONTINUE
                        #140
                        TRIX(I2R,0,MR,A,BB,C,B,TCOS,D,W)
                        IP = IP+MR
                        IPSTOR = MAX0(IPSTOR,IP+MR)
                        for I in range(1,MR+1):
                            II = IP+I
                            P[II] = B[I]+np.double(0.500)*(Q[I][J]-Q[I][JM1]-Q[I][JP1])
                            B[I] =P[II]+Q[I][JP2]
                        #CONTINUE
                        #141
                        if (LR != 0):#--->142
                            COSGEN(LR,1,np.double(0.500),FDEN,TCOS(I2R+1))
                            MERGE(TCOS,0,I2R,I2R,LR,KR)
                            #---->144
                        else:
                            #142
                            for I in range(1,MR+1):
                                II = KR+I
                                TCOS[II] = TCOS[I]
                            #COTINUE
                            #143
                            #--->144
                        #144
                        COSGEN(KR,1,np.double(0.500),FDEN,TCOS)
                        if(LR == 0):#--->145
                            if(ISTAG == 1):#--->146
                                #146
                                for I in range(1,MR+1):
                                    B[I] = FISTAG*B[I]
                                #CONTINUE
                                #147
                                #---->148
                            if(ISTAG == 2): #----->145
                                #145
                                TRIX(KR,KR,MR,A,BB,C,B,TCOS,D,W)
                                #---->148
                        else:
                            #145
                            TRIX(KR,KR,MR,A,BB,C,B,TCOS,D,W)
                            #--->148

                        #148
                        for I in range(1,MR+1):
                            II = IP+I
                            Q[I][J] = Q[I][JM2]+P[II]+B[I]
                        #149
                        #COTINUE
                    #150
                    LR = KR
                    KR = KR+JR
                    #CONTINUE  
                #151
                if(MIXBND == 1):#--->152
                    #152
                    NR = (NLAST-1)/JR+1
                    if(NR <= 3):#---155
                        #155
                        JM1 = J-I2R
                        JP1 = J+I2R
                        JM2 = NLAST-I2R
                        if(NR == 2):#---->184
                            #184
                            if(N == 2):#--->188
                                for I in range(1,MR+1):
                                    B[I] = Q[I][1]
                                #185
                                #CONTINUE
                                TCOS[1] = 0
                                TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                                for I in range(1,MR+1):
                                    Q[I][1] = B[I]
                                    B[I] = np.double(2.00)*(Q[I][2]+B[I])*FISTAG
                                #186
                                #CONTINUE
                                TCOS[1] = -FISTAG
                                TCOS[2] = np.double(2.00)
                                TRIX(2,0,MR,A,BB,C,B,TCOS,D,W)
                                for I in range(1,MR+1):
                                    Q[I][1] = Q[I][1]+B[I]
                                #187
                                #CONTINUE
                                JR = 1
                                I2R = 0
                                #--->194
                            else:
                                #188
                                for I in range(1,MR+1):
                                    II = IP+1
                                    B3[I] = np.double(0.00)
                                    B[I] = Q[I][1]+np.double(2.00)*P[II]
                                    Q[I][1] = np.double(0.500)*Q[I][1]-Q[I][JM1]
                                    B2[I] = np.double(2.00)*(Q[I][1]+Q[I][NLAST])
                                #189
                                K1 = KR+JR-1
                                TCOS[K1+1] = np.double(-2.00)
                                K4 = K1+3-ISTAG
                                COSGEN(KR+ISTAG-2,1,np.double(0.00),FNUM,TCOS(K4))
                                K4 = K1+KR+1
                                COSGEN(JR-1,1,np.double(0.00),np.double(1.00),TCOS(K4))
                                MERGE(TCOS,K1,KR,K1+KR,JR-1,0)
                                COSGEN(KR,1,np.double(0.500),FDEN,TCOS(K1+1))
                                K2 = KR
                                K4 = K1+K2+1
                                COSGEN(LR,1,np.double(0.500),FDEN,TCOS(K4))
                                K3 = LR
                                K4 = 0
                                TRI3(MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
                                for I in range(1,MR+1):
                                    B[I] = B[I]+B2[I]
                                #COTINUE
                                #190
                                TCOS[1] = np.double(2.00)
                                TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                                for I in range(1,MR+1):
                                    Q[I][1] = Q[I][NLAST]
                                #191
                                # --->194
                            #194
                            J = NLAST - JR
                            for I in range(1,MR+1):
                                B[I] = Q[I][NLAST]+Q[I][J]
                            #195
                            #CONTINUE
                            #--->196
                        if(LR != 0):#--->170
                            #170
                            for I in range(1,MR+1):
                                B[I] = np.double(0.500)*Q[I][1]-Q[I][JM1]+Q[I][J]
                            #171
                            if(NROD == 0):#---->173
                                for I in range(1,MR+1):
                                    II = IP+I
                                    B[I] = B[I]+P[II]
                                    #172
                                    #CONTINUE
                                    #--->175
                            else:
                                #173
                                for I in range(1,MR+1):
                                    B[I] = B[I]+Q[I][NLAST]-Q[I][JM2]
                                #174
                                #CONTINUE
                                #--->175
                            #175
                            for I in range(1,MR+1):
                                T = np.double(0.500)*(Q[I][J]-Q[I][JM1]-Q[I][JP1])
                                Q[I][J] = T
                                B2[I] = Q[I][NLAST]+T
                                B3[I] = Q[I][1]+np.double(2.00)*T
                            #176
                            #177
                            K1 = KR+2*JR-1
                            K2 = KR+JR
                            TCOS[K1+1] = np.double(-2.00)
                            K4 = K1+3-ISTAG
                            COSGEN(K2+ISTAG-2,1,np.double(0.00),FNUM,TCOS(K4))
                            K4 = K1+K2+1
                            COSGEN(JR-1,1,np.double(0.00),np.double(1.00),TCOS(K4))
                            MERGE(TCOS,K1,K2,K1+K2,JR-1,0)
                            K3 = K1+K2+LR
                            COSGEN(JR,1,np.double(0.500),np.double(0.00),TCOS(K3+1))
                            K4 = K3+JR+1
                            COSGEN(KR,1,np.double(0.500),FDEN,TCOS(K4))
                            MERGE(TCOS,K3,JR,K3+JR,KR,K1)
                            if(LR != 0):#---->178
                                COSGEN(LR,1,np.double(0.500),FDEN,TCOS(K4))
                                MERGE(TCOS,K3,JR,K3+JR,LR,K3-LR)
                                COSGEN(KR,1,np.double(0.500),FDEN,TCOS(K4))
                            #178
                            K3 = KR
                            K4 = KR
                            TRI3(MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
                            for I in range(1,MR+1):
                                B[I] = B[I]+B2[I]+B3[I]
                            #179
                            TCOS[1] = np.double(2.00)
                            TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                            for I in range(1,MR+1):
                                Q[I][J] = Q[I][J] +B[I]
                                B[I] = Q[I][1]+np.double(2.00)*Q[I][J]
                            #180
                            COSGEN(JR,1,np.double(0.500),np.double(0.00),TCOS)
                            TRIX(JR,0,MR,A,BB,C,B,TCOS,D,W)
                            if(JR == 1):#--->182
                                for I in range(1,MR+1):
                                    Q[I][1] = B[I]
                                #181
                                #CONTINUE
                                #---->194
                            else:
                                #182
                                for I in range(1,MR+1):
                                    Q[I][1] = np.double(0.500)*Q[I][JM1]+B[I]
                                #183
                                #CONTINUE
                                #--->194
                            #194
                            J = NLAST - JR
                            for I in range(1,MR+1):
                                B[I] = Q[I][NLAST]+Q[I][J]
                            #195
                            #CONTINUE
                            #--->196
                        if(N != 3):#--->161
                            #161
                            if(ISTAG == 1):#--->162
                                #162
                                for I in range(1,MR+1):
                                    B[I] = Q[I][J]+np.double(0.500)*Q[I][1]-Q[I][JM1]+Q[I][NLAST]-Q[I][JM2]
                                #163
                                COSGEN(JR,1,np.double(0.500),np.double(0.00),TCOS)
                                TRIX(JR,0,MR,A,BB,C,B,TCOS,D,W)
                                for I in range(1,MR+1):
                                    Q[I][J] = np.double(0.500)*(Q[I][J]-Q[I][JM1]-Q[I][JP1])+B[I]
                                    B[I] = Q[I][1]+np.double(2.00)*Q[I][NLAST]+np.double(4.00)*Q[I][J]
                                #164
                                JR2 = 2*JR
                                COSGEN(JR,1,np.double(0.00),np.double(0.00),TCOS)
                                for I in range(1,JR):
                                    I1 = JR+1
                                    I2 = JR+1-I
                                    TCOS[I1] = -TCOS[I2]
                                #165
                                #CONTINUE
                                TRIX(JR2,0,MR,A,BB,C,B,TCOS,D,W)
                                for I in range(1,MR+1):
                                    Q[I][J] = Q[I][J]+B[I]
                                    B[I] = Q[I][1]+np.double(2.00)*Q[I][J]
                                #166
                                COSGEN(KR,1,np.double(0.500),np.double(0.00),TCOS)
                                TRIX(KR,0,MR,A,BB,C,B,TCOS,D,W)
                                for I in range(1,MR+1):
                                    Q[I][1] =np.double(0.500)*Q[I][1]-Q[I][JM1]+B[I]
                                #COTINUE
                                #167
                                #---->194
                            if(ISTAG == 2):#--->170
                                #170
                                for I in range(1,MR+1):
                                    B[I] = np.double(0.500)*Q[I][1]-Q[I][JM1]+Q[I][J]
                                #171
                                if(NROD == 0):#---->173
                                    for I in range(1,MR+1):
                                        II = IP+I
                                        B[I] = B[I]+P[II]
                                        #172
                                        #CONTINUE
                                        #--->175
                                else:
                                    #173
                                    for I in range(1,MR+1):
                                        B[I] = B[I]+Q[I][NLAST]-Q[I][JM2]
                                    #174
                                    #CONTINUE
                                    #--->175
                                #175
                                for I in range(1,MR+1):
                                    T = np.double(0.500)*(Q[I][J]-Q[I][JM1]-Q[I][JP1])
                                    Q[I][J] = T
                                    B2[I] = Q[I][NLAST]+T
                                    B3[I] = Q[I][1]+np.double(2.00)*T
                                #176
                                #177
                                K1 = KR+2*JR-1
                                K2 = KR+JR
                                TCOS[K1+1] = np.double(-2.00)
                                K4 = K1+3-ISTAG
                                COSGEN(K2+ISTAG-2,1,np.double(0.00),FNUM,TCOS(K4))
                                K4 = K1+K2+1
                                COSGEN(JR-1,1,np.double(0.00),np.double(1.00),TCOS(K4))
                                MERGE(TCOS,K1,K2,K1+K2,JR-1,0)
                                K3 = K1+K2+LR
                                COSGEN(JR,1,np.double(0.500),np.double(0.00),TCOS(K3+1))
                                K4 = K3+JR+1
                                COSGEN(KR,1,np.double(0.500),FDEN,TCOS(K4))
                                MERGE(TCOS,K3,JR,K3+JR,KR,K1)
                                if(LR != 0):#---->178
                                    COSGEN(LR,1,np.double(0.500),FDEN,TCOS(K4))
                                    MERGE(TCOS,K3,JR,K3+JR,LR,K3-LR)
                                    COSGEN(KR,1,np.double(0.500),FDEN,TCOS(K4))
                                #178
                                K3 = KR
                                K4 = KR
                                TRI3(MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
                                for I in range(1,MR+1):
                                    B[I] = B[I]+B2[I]+B3[I]
                                #179
                                TCOS[1] = np.double(2.00)
                                TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                                for I in range(1,MR+1):
                                    Q[I][J] = Q[I][J] +B[I]
                                    B[I] = Q[I][1]+np.double(2.00)*Q[I][J]
                                #180
                                COSGEN(JR,1,np.double(0.500),np.double(0.00),TCOS)
                                TRIX(JR,0,MR,A,BB,C,B,TCOS,D,W)
                                if(JR == 1):#--->182
                                    for I in range(1,MR+1):
                                        Q[I][1] = B[I]
                                    #181
                                    #CONTINUE
                                    #---->194
                                else:
                                    #182
                                    for I in range(1,MR+1):
                                        Q[I][1] = np.double(0.500)*Q[I][JM1]+B[I]
                                    #183
                                    #CONTINUE
                                    #--->194
                            #194
                            J = NLAST - JR
                            for I in range(1,MR+1):
                                B[I] = Q[I][NLAST]+Q[I][J]
                            #195
                            #CONTINUE
                            #--->196
                        else:
                            #156
                            for I in range(1,MR+1):
                                B[I] = Q[I][2]
                            #CONTINUE
                            #157
                            TCOS[1] = np.double(0.00)
                            TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                            for I in range(1,MR+1):
                                Q[I][2] = B[I]
                                B[I] = np.double(4.00)*B[I][1]+np.double(2.00)*Q[I][3]
                            #158
                            #CONTINUE

                            TCOS[1] = np.double(-2.00)
                            TCOS[2] = np.double(2.00)
                            I1 = 0
                            I2 = 0
                            TRIX(I1,I2,MR,A,BB,C,B,TCOS,D,W)
                            for I in range(1,MR+1):
                                Q[I][2] = Q[I][2] + B[I]
                                B[I] = Q[I][1] + np.double(2.00)*Q[I][2]
                            #159
                            TCOS[1] = np.double(0.00)
                            TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                            for I in range(1,MR+1):
                                Q[I][1] = B[I]
                            #160
                            JR = 1
                            I2R = 0
                            #--->194
                            #194
                            J = NLAST - JR
                            for I in range(1,MR+1):
                                B[I] = Q[I][NLAST]+Q[I][J]
                            #195
                            #CONTINUE
                            #--->196

                        break
                    else:
                        I2R = JR
                        NRODPR = NROD
                if(MIXBND == 2):#---->153
                    #153
                    NR = NLAST/JR
                    if(NR <= 1):#--->192
                        #192
                        for I in range(1,MR+1):
                            B[I] = Q[I][NLAST]
                        #--->196
                        break
                    else:
                        I2R = JR
                        NRODPR = NROD  
            #196
            JM2 = NLAST  - I2R
            if(JR == 1):#--->198
                for I in range(1,MR+1):
                    Q[I][NLAST] = np.double(0.000)
                #197
                #CONTINUE
                #--->202
            else:
                #198
                if (NROD == 0):#--->200
                    for I in range(1,MR+1):
                        II = IP+I
                        Q[I][NLAST] = P[II]
                    #199
                    IP = IP+MR
                    #---->202
                else:
                    #200
                    for I in range(1,MR+1):
                        Q[I][NLAST] = Q[I][NLAST]-Q[I][JM2]
                    #201
                    #CONTINUE
                    #--->202
            #202
            COSGEN (KR,1,np.double(0.500),FDEN,TCOS)
            COSGEN (LR,1,np.double(0.500),FDEN,TCOS(KR+1))
            if (LR == 0):#------>204
                for I in range(1,MR+1):
                    B[I] = FISTAG*B[I]
            #203
            #CONTINUE
            #204
            TRIX (KR,LR,MR,A,BB,C,B,TCOS,D,W)
            for I in range(1,MR+1):
                Q[I][NLAST] = Q[I][NLAST]+B
            #CONTIUE
            #205
            NLASTP = NLAST
            #206
            JSTEP = JR
            JR = I2R
            I2R = I2R/2
            while(JR > 0):
                #
                if (JR != 0):#--->222
                    if(MIXBND == 1):#---->207
                        #207
                        JSTART = 1+JR
                        #---->209
                    if(MIXBND == 2):#---->208
                        #208
                        JSTART = JR
                    #209
                    KR = KR-JR
                    if (NLAST+JR < N):#----->210
                        KR = KR-JR
                        NLAST = NLAST+JR
                        JSTOP = NLAST-JSTEP
                        #---->211
                    else:
                        #210
                        JSTOP = NLAST-JR
                        #--->211
                    #211
                    LR = KR-JR
                    COSGEN (JR,1,np.double(0.500),np.double(0.000),TCOS)
                    for J in range(JSTART,JSTOP,JSTEP):
                        JM2 = J-JR
                        JP2 = J+JR
                        if(J == JR):#--->213
                            for I in range(1,MR+1):
                                B[I] = Q[I][J]+Q[I][JP2]
                            #212
                            #CONTINUE
                            #215
                        else:
                            #213
                            for I in range(1,MR+1):
                                B[I] = Q[I][J]+Q[I][JM2]+Q[I][JP2]
                            #214
                            #CONTINUE
                        #215
                        if(JR != 1):#---->217
                            for I in range(1,MR+1):
                                Q[I][J] = np.double(0.00)
                            #216
                            #CONTINUE
                            #---->219
                        else:
                            #217
                            JM1 = J-I2R
                            JP1 = J+I2R
                            for I in range(1,MR+1):
                                Q[I][J] = np.double(0.500)*(Q[I][J]-Q[I][JM1]-Q[I][JP1])
                            #218
                            #CONTINUE
                        #219
                        TRIX (JR,0,MR,A,BB,C,B,TCOS,D,W)
                        for I in range(1,MR+1):
                            Q[I][J] = Q[I][J]+B[I]
                            #220
                            #CONTINUE
                    #221
                    #CONTINUE
                    NROD = 1
                    if(NLAST+I2R > N):
                        NROD = 0
                    if(NLASTP != NLAST):
                        #--->194
                        #194
                        J = NLAST - JR
                        for I in range(1,MR+1):
                            B[I] = Q[I][NLAST]+Q[I][J]
                        #195
                        #CONTINUE
                        #--->196
                        #196
                        JM2 = NLAST  - I2R
                        if(JR == 1):#--->198
                            for I in range(1,MR+1):
                                Q[I][NLAST] = np.double(0.000)
                            #197
                            #CONTINUE
                            #--->202
                        else:
                            #198
                            if (NROD == 0):#--->200
                                for I in range(1,MR+1):
                                    II = IP+I
                                    Q[I][NLAST] = P[II]
                                #199
                                IP = IP+MR
                                #---->202
                            else:
                                #200
                                for I in range(1,MR+1):
                                    Q[I][NLAST] = Q[I][NLAST]-Q[I][JM2]
                                #201
                                #CONTINUE
                                #--->202
                        #202
                        COSGEN (KR,1,np.double(0.500),FDEN,TCOS)
                        COSGEN (LR,1,np.double(0.500),FDEN,TCOS(KR+1))
                        if (LR == 0):#------>204
                            for I in range(1,MR+1):
                                B[I] = FISTAG*B[I]
                        #203
                        #CONTINUE
                        #204
                        TRIX (KR,LR,MR,A,BB,C,B,TCOS,D,W)
                        for I in range(1,MR+1):
                            Q[I][NLAST] = Q[I][NLAST]+B
                        #CONTIUE
                        #205
                        NLASTP = NLAST
                        #206
                        JSTEP = JR
                        JR = I2R
                        I2R = I2R/2
                    else:
                        #206
                        JSTEP = JR
                        JR = I2R
                        I2R = I2R/2
                else:
                    break
            #222
            W[1] = IPSTOR
            return














    



    if(ISTAG == 2):#---->103
        #103
        if(N > 3):#--->155
            while(NR > 1):
                #104
                JR = 2*I2R
                NROD = 0
                if ((NR/2)*2 == NR):
                    NROD=0
                if(MIXBND == 1):
                        #105
                        JSTART = 1
                if (MIXBND == 2):
                    #106
                    JSTART = JR
                    NROD = 1-NROD
                #107
                JSTOP = NLAST-JR
                if (NROD == 0):
                    JSTOP=JSTOP-I2R
                COSGEN(I2R,1,np.double(0.500),np.double(0.00),TCOS)
                I2RBY2= I2R/2
                if(JSTOP < JSTART):#--->108
                    J = JR
                    #--->116
                #108
                else:
                    for J in range(JSTART,JSTOP,JR):
                        JP1 = J+I2RBY2
                        JP2 = J+I2R
                        JP3 = JP2+I2RBY2
                        JM1 = J-I2RBY2
                        JM2 = J-I2R
                        JM3 = JM2-I2RBY2
                        if (J == 1):#-->109
                            JM1 = JP1
                            JM2 = JP2
                            JM3 = JP3
                        else:
                            #109
                            if(I2R == 1):#---->111
                                if(J == 1):
                                    JM2 = JP2
                                for I in range(1,MR+1):
                                    B[I] = np.double(2.00)*Q[I][J]
                                    Q[I][J] = Q[I][JM2]+Q[I][JP2]
                                #110
                                #--->113
                            else:
                                #111
                                for I in range(1,MR+1):
                                    FI = Q[I][J]
                                    Q[I][J] = Q[I][J]-Q[I][JM1]-Q[I][JP1]+Q[I][JM2]+Q[I][JP2]
                                    B[I] = FI+Q[I][J]-Q[I][JM3]-Q[I][JP3]
                                #112
                                #CONTINUE
                            #113
                            TRIX(I2R,0,MR,A,BB,C,B,TCOS,D,W)
                            for I in range(1,MR+1):
                                Q[I][J] = Q[I][J]+B[I]
                            #114
                            #COTINUE
                            #115
                            #CONTINUE
                            J = JSTOP+JR
                #116
                NLAST = J
                JM1 = J-I2RBY2
                JM2 = J-I2R
                JM3 = JM2-I2RBY2
                if(NROD != 0):#---->128
                    if(I2R == 1):#--->118
                        for I in range(1,MR+1):
                            B[I] = FISTAG*Q[I][J]
                            Q[I][J] = Q[I][JM2]
                        #117
                        #CONTINUE
                        #---->126
                    else:
                        #118
                        for I in range(1,MR+1):
                            B[I] = Q[I][J]+np.double(0.500)*(Q[I][JM2]-Q[I][JM1]-Q[I][JM3])
                        #119
                        #CONTINUE
                        if(NRODPR == 0):#--->121
                            for I in range(1,MR+1):
                                II = IP+I
                                Q[I][J] = Q[I][JM2]+P[II]
                                #120
                                #CONTINUE
                                IP = IP-MR
                                #--->123
                        else:
                            #121
                            for I in range(1,MR+1):
                                Q[I][J] = Q[I][J]-Q[I][JM1]+Q[I][JM2]
                            #122
                            #CONTINUE
                        #123
                        if(LR != 0):#---->124
                            COSGEN(LR,1,0,np.double(0.500),FDEN,TCOS(KR+1))
                            #---->126
                        else:
                            #124
                            for I in range(1,MR+1):
                                B[I] = FISTAG*B[I]
                            #125
                            #CONTINUE
                            #---->126
                    #126
                    COSGEN(KR,1,np.double(0.500),FDEN,TCOS)
                    TRIX(KR,LR,MR,A,BB,C,B,TCOS,D,W)
                    for I in range(1,MR+1):
                        Q[I][J] = Q[I][J]+B[I]
                    #127
                    KR = KR+I2R
                    #---->151
                else:
                    #128
                    JP1 = J+I2RBY2
                    JP2 = J+I2R
                    if(I2R == 1):#---->135
                        for I in range(1,MR+1):
                            B[I] = Q[I][J]
                        TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                        IP = 0
                        IPSTOR = MR
                        if(ISTAG == 1):#--->133
                            #133
                            for I in range(1,MR+1):
                                P[I] = B[I]
                                Q[I][J] = Q[I][JM2]+np.double(2.00)*Q[I][JP2]+np.double(3.00)*B[I]
                            #134
                            #COTINUE
                            #---->150
                        if(ISTAG == 2):#--->130
                            #130
                            for I in range(1,MR+1):
                                P[I] = B[I]
                                B[I] = B[I]+Q[I][N]
                            #131
                            #CONTINUE
                            TCOS[1] =  np.double(1.00)
                            TCOS[2] = np.double(0.00)
                            TRIX(1,1,MR,A,BB,C,B,TCOS,D,W)
                            for I in range(1,MR+1):
                                Q[I][J] = Q[I][JM2]+P[I]+B[I]
                            #132
                            #CONTINUE
                            #---->150
                    else:
                        #135
                        for I in range(1,MR+1):
                            B[I] = Q[I][J]+np.double(0.500)*(Q[I][JM2]-Q[I][JM1]-Q[I][JM3])
                        #136
                        if (NRODPR == 0):#--->138
                            for I in range(1,MR+1):
                                II = IP+I
                                B[I] = B[I]+P[II]
                            #137
                            #CONTINUE
                            #--->140
                        else:
                            #138
                            for I in range(1,MR+1):
                                B[I] = B[I]+Q[I][JP2]-Q[I][JP1]
                            #139
                            #CONTINUE
                        #140
                        TRIX(I2R,0,MR,A,BB,C,B,TCOS,D,W)
                        IP = IP+MR
                        IPSTOR = MAX0(IPSTOR,IP+MR)
                        for I in range(1,MR+1):
                            II = IP+I
                            P[II] = B[I]+np.double(0.500)*(Q[I][J]-Q[I][JM1]-Q[I][JP1])
                            B[I] =P[II]+Q[I][JP2]
                        #CONTINUE
                        #141
                        if (LR != 0):#--->142
                            COSGEN(LR,1,np.double(0.500),FDEN,TCOS(I2R+1))
                            MERGE(TCOS,0,I2R,I2R,LR,KR)
                            #---->144
                        else:
                            #142
                            for I in range(1,MR+1):
                                II = KR+I
                                TCOS[II] = TCOS[I]
                            #COTINUE
                            #143
                            #--->144
                        #144
                        COSGEN(KR,1,np.double(0.500),FDEN,TCOS)
                        if(LR == 0):#--->145
                            if(ISTAG == 1):#--->146
                                #146
                                for I in range(1,MR+1):
                                    B[I] = FISTAG*B[I]
                                #CONTINUE
                                #147
                                #---->148
                            if(ISTAG == 2): #----->145
                                #145
                                TRIX(KR,KR,MR,A,BB,C,B,TCOS,D,W)
                                #---->148
                        else:
                            #145
                            TRIX(KR,KR,MR,A,BB,C,B,TCOS,D,W)
                            #--->148

                        #148
                        for I in range(1,MR+1):
                            II = IP+I
                            Q[I][J] = Q[I][JM2]+P[II]+B[I]
                        #149
                        #COTINUE
                    #150
                    LR = KR
                    KR = KR+JR
                    #CONTINUE
                #151
                if(MIXBND == 1):#--->152
                    #152
                    NR = (NLAST-1)/JR+1
                    if(NR <= 3):#---155
                        #155
                        JM1 = J-I2R
                        JP1 = J+I2R
                        JM2 = NLAST-I2R
                        if(NR == 2):#---->184
                            #184
                            if(N == 2):#--->188
                                for I in range(1,MR+1):
                                    B[I] = Q[I][1]
                                #185
                                #CONTINUE
                                TCOS[1] = 0
                                TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                                for I in range(1,MR+1):
                                    Q[I][1] = B[I]
                                    B[I] = np.double(2.00)*(Q[I][2]+B[I])*FISTAG
                                #186
                                #CONTINUE
                                TCOS[1] = -FISTAG
                                TCOS[2] = np.double(2.00)
                                TRIX(2,0,MR,A,BB,C,B,TCOS,D,W)
                                for I in range(1,MR+1):
                                    Q[I][1] = Q[I][1]+B[I]
                                #187
                                #CONTINUE
                                JR = 1
                                I2R = 0
                                #--->194
                            else:
                                #188
                                for I in range(1,MR+1):
                                    II = IP+1
                                    B3[I] = np.double(0.00)
                                    B[I] = Q[I][1]+np.double(2.00)*P[II]
                                    Q[I][1] = np.double(0.500)*Q[I][1]-Q[I][JM1]
                                    B2[I] = np.double(2.00)*(Q[I][1]+Q[I][NLAST])
                                #189
                                K1 = KR+JR-1
                                TCOS[K1+1] = np.double(-2.00)
                                K4 = K1+3-ISTAG
                                COSGEN(KR+ISTAG-2,1,np.double(0.00),FNUM,TCOS(K4))
                                K4 = K1+KR+1
                                COSGEN(JR-1,1,np.double(0.00),np.double(1.00),TCOS(K4))
                                MERGE(TCOS,K1,KR,K1+KR,JR-1,0)
                                COSGEN(KR,1,np.double(0.500),FDEN,TCOS(K1+1))
                                K2 = KR
                                K4 = K1+K2+1
                                COSGEN(LR,1,np.double(0.500),FDEN,TCOS(K4))
                                K3 = LR
                                K4 = 0
                                TRI3(MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
                                for I in range(1,MR+1):
                                    B[I] = B[I]+B2[I]
                                #COTINUE
                                #190
                                TCOS[1] = np.double(2.00)
                                TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                                for I in range(1,MR+1):
                                    Q[I][1] = Q[I][NLAST]
                                #191
                                # --->194
                            #194
                            J = NLAST - JR
                            for I in range(1,MR+1):
                                B[I] = Q[I][NLAST]+Q[I][J]
                            #195
                            #CONTINUE
                            #--->196
                        if(LR != 0):#--->170
                            #170
                            for I in range(1,MR+1):
                                B[I] = np.double(0.500)*Q[I][1]-Q[I][JM1]+Q[I][J]
                            #171
                            if(NROD == 0):#---->173
                                for I in range(1,MR+1):
                                    II = IP+I
                                    B[I] = B[I]+P[II]
                                    #172
                                    #CONTINUE
                                    #--->175
                            else:
                                #173
                                for I in range(1,MR+1):
                                    B[I] = B[I]+Q[I][NLAST]-Q[I][JM2]
                                #174
                                #CONTINUE
                                #--->175
                            #175
                            for I in range(1,MR+1):
                                T = np.double(0.500)*(Q[I][J]-Q[I][JM1]-Q[I][JP1])
                                Q[I][J] = T
                                B2[I] = Q[I][NLAST]+T
                                B3[I] = Q[I][1]+np.double(2.00)*T
                            #176
                            #177
                            K1 = KR+2*JR-1
                            K2 = KR+JR
                            TCOS[K1+1] = np.double(-2.00)
                            K4 = K1+3-ISTAG
                            COSGEN(K2+ISTAG-2,1,np.double(0.00),FNUM,TCOS(K4))
                            K4 = K1+K2+1
                            COSGEN(JR-1,1,np.double(0.00),np.double(1.00),TCOS(K4))
                            MERGE(TCOS,K1,K2,K1+K2,JR-1,0)
                            K3 = K1+K2+LR
                            COSGEN(JR,1,np.double(0.500),np.double(0.00),TCOS(K3+1))
                            K4 = K3+JR+1
                            COSGEN(KR,1,np.double(0.500),FDEN,TCOS(K4))
                            MERGE(TCOS,K3,JR,K3+JR,KR,K1)
                            if(LR != 0):#---->178
                                COSGEN(LR,1,np.double(0.500),FDEN,TCOS(K4))
                                MERGE(TCOS,K3,JR,K3+JR,LR,K3-LR)
                                COSGEN(KR,1,np.double(0.500),FDEN,TCOS(K4))
                            #178
                            K3 = KR
                            K4 = KR
                            TRI3(MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
                            for I in range(1,MR+1):
                                B[I] = B[I]+B2[I]+B3[I]
                            #179
                            TCOS[1] = np.double(2.00)
                            TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                            for I in range(1,MR+1):
                                Q[I][J] = Q[I][J] +B[I]
                                B[I] = Q[I][1]+np.double(2.00)*Q[I][J]
                            #180
                            COSGEN(JR,1,np.double(0.500),np.double(0.00),TCOS)
                            TRIX(JR,0,MR,A,BB,C,B,TCOS,D,W)
                            if(JR == 1):#--->182
                                for I in range(1,MR+1):
                                    Q[I][1] = B[I]
                                #181
                                #CONTINUE
                                #---->194
                            else:
                                #182
                                for I in range(1,MR+1):
                                    Q[I][1] = np.double(0.500)*Q[I][JM1]+B[I]
                                #183
                                #CONTINUE
                                #--->194
                            #194
                            J = NLAST - JR
                            for I in range(1,MR+1):
                                B[I] = Q[I][NLAST]+Q[I][J]
                            #195
                            #CONTINUE
                            #--->196
                        if(N != 3):#--->161
                            #161
                            if(ISTAG == 1):#--->162
                                #162
                                for I in range(1,MR+1):
                                    B[I] = Q[I][J]+np.double(0.500)*Q[I][1]-Q[I][JM1]+Q[I][NLAST]-Q[I][JM2]
                                #163
                                COSGEN(JR,1,np.double(0.500),np.double(0.00),TCOS)
                                TRIX(JR,0,MR,A,BB,C,B,TCOS,D,W)
                                for I in range(1,MR+1):
                                    Q[I][J] = np.double(0.500)*(Q[I][J]-Q[I][JM1]-Q[I][JP1])+B[I]
                                    B[I] = Q[I][1]+np.double(2.00)*Q[I][NLAST]+np.double(4.00)*Q[I][J]
                                #164
                                JR2 = 2*JR
                                COSGEN(JR,1,np.double(0.00),np.double(0.00),TCOS)
                                for I in range(1,JR):
                                    I1 = JR+1
                                    I2 = JR+1-I
                                    TCOS[I1] = -TCOS[I2]
                                #165
                                #CONTINUE
                                TRIX(JR2,0,MR,A,BB,C,B,TCOS,D,W)
                                for I in range(1,MR+1):
                                    Q[I][J] = Q[I][J]+B[I]
                                    B[I] = Q[I][1]+np.double(2.00)*Q[I][J]
                                #166
                                COSGEN(KR,1,np.double(0.500),np.double(0.00),TCOS)
                                TRIX(KR,0,MR,A,BB,C,B,TCOS,D,W)
                                for I in range(1,MR+1):
                                    Q[I][1] =np.double(0.500)*Q[I][1]-Q[I][JM1]+B[I]
                                #COTINUE
                                #167
                                #---->194
                            if(ISTAG == 2):#--->170
                                #170
                                for I in range(1,MR+1):
                                    B[I] = np.double(0.500)*Q[I][1]-Q[I][JM1]+Q[I][J]
                                #171
                                if(NROD == 0):#---->173
                                    for I in range(1,MR+1):
                                        II = IP+I
                                        B[I] = B[I]+P[II]
                                        #172
                                        #CONTINUE
                                        #--->175
                                else:
                                    #173
                                    for I in range(1,MR+1):
                                        B[I] = B[I]+Q[I][NLAST]-Q[I][JM2]
                                    #174
                                    #CONTINUE
                                    #--->175
                                #175
                                for I in range(1,MR+1):
                                    T = np.double(0.500)*(Q[I][J]-Q[I][JM1]-Q[I][JP1])
                                    Q[I][J] = T
                                    B2[I] = Q[I][NLAST]+T
                                    B3[I] = Q[I][1]+np.double(2.00)*T
                                #176
                                #177
                                K1 = KR+2*JR-1
                                K2 = KR+JR
                                TCOS[K1+1] = np.double(-2.00)
                                K4 = K1+3-ISTAG
                                COSGEN(K2+ISTAG-2,1,np.double(0.00),FNUM,TCOS(K4))
                                K4 = K1+K2+1
                                COSGEN(JR-1,1,np.double(0.00),np.double(1.00),TCOS(K4))
                                MERGE(TCOS,K1,K2,K1+K2,JR-1,0)
                                K3 = K1+K2+LR
                                COSGEN(JR,1,np.double(0.500),np.double(0.00),TCOS(K3+1))
                                K4 = K3+JR+1
                                COSGEN(KR,1,np.double(0.500),FDEN,TCOS(K4))
                                MERGE(TCOS,K3,JR,K3+JR,KR,K1)
                                if(LR != 0):#---->178
                                    COSGEN(LR,1,np.double(0.500),FDEN,TCOS(K4))
                                    MERGE(TCOS,K3,JR,K3+JR,LR,K3-LR)
                                    COSGEN(KR,1,np.double(0.500),FDEN,TCOS(K4))
                                #178
                                K3 = KR
                                K4 = KR
                                TRI3(MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
                                for I in range(1,MR+1):
                                    B[I] = B[I]+B2[I]+B3[I]
                                #179
                                TCOS[1] = np.double(2.00)
                                TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                                for I in range(1,MR+1):
                                    Q[I][J] = Q[I][J] +B[I]
                                    B[I] = Q[I][1]+np.double(2.00)*Q[I][J]
                                #180
                                COSGEN(JR,1,np.double(0.500),np.double(0.00),TCOS)
                                TRIX(JR,0,MR,A,BB,C,B,TCOS,D,W)
                                if(JR == 1):#--->182
                                    for I in range(1,MR+1):
                                        Q[I][1] = B[I]
                                    #181
                                     #CONTINUE
                                    #---->194
                                else:
                                    #182
                                    for I in range(1,MR+1):
                                        Q[I][1] = np.double(0.500)*Q[I][JM1]+B[I]
                                    #183
                                    #CONTINUE
                                    #--->194
                            #194
                            J = NLAST - JR
                            for I in range(1,MR+1):
                                B[I] = Q[I][NLAST]+Q[I][J]
                            #195
                            #CONTINUE
                            #--->196
                        else:
                            #156
                            for I in range(1,MR+1):
                                B[I] = Q[I][2]
                            #CONTINUE
                            #157
                            TCOS[1] = np.double(0.00)
                            TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                            for I in range(1,MR+1):
                                Q[I][2] = B[I]
                                B[I] = np.double(4.00)*B[I][1]+np.double(2.00)*Q[I][3]
                            #158
                            #CONTINUE
                            TCOS[1] = np.double(-2.00)
                            TCOS[2] = np.double(2.00)
                            I1 = 0
                            I2 = 0
                            TRIX(I1,I2,MR,A,BB,C,B,TCOS,D,W)
                            for I in range(1,MR+1):
                                Q[I][2] = Q[I][2] + B[I]
                                B[I] = Q[I][1] + np.double(2.00)*Q[I][2]
                            #159
                            TCOS[1] = np.double(0.00)
                            TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                            for I in range(1,MR+1):
                                Q[I][1] = B[I]
                            #160
                            JR = 1
                            I2R = 0
                            #--->194
                            #194
                            J = NLAST - JR
                            for I in range(1,MR+1):
                                B[I] = Q[I][NLAST]+Q[I][J]
                            #195
                            #CONTINUE
                            #--->196
                        break
                    else:
                        I2R = JR
                        NRODPR = NROD
                if(MIXBND == 2):#---->153
                    #153
                    NR = NLAST/JR
                    if(NR <= 1):#--->192
                        #192
                        for I in range(1,MR+1):
                            B[I] = Q[I][NLAST]
                        #--->196
                        break
                    else:
                        I2R = JR
                        NRODPR = NROD
        else:
            #155
            JM1 = J-I2R
            JP1 = J+I2R
            JM2 = NLAST-I2R
            if(NR == 2):#---->184
                #184
                if(N == 2):#--->188
                    for I in range(1,MR+1):
                        B[I] = Q[I][1]
                    #185
                    #CONTINUE
                    TCOS[1] = 0
                    TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                    for I in range(1,MR+1):
                        Q[I][1] = B[I]
                        B[I] = np.double(2.00)*(Q[I][2]+B[I])*FISTAG
                    #186
                    #CONTINUE
                    TCOS[1] = -FISTAG
                    TCOS[2] = np.double(2.00)
                    TRIX(2,0,MR,A,BB,C,B,TCOS,D,W)
                    for I in range(1,MR+1):
                        Q[I][1] = Q[I][1]+B[I]
                    #187
                    #CONTINUE
                    JR = 1
                    I2R = 0
                    #--->194
                else:
                    #188
                    for I in range(1,MR+1):
                        II = IP+1
                        B3[I] = np.double(0.00)
                        B[I] = Q[I][1]+np.double(2.00)*P[II]
                        Q[I][1] = np.double(0.500)*Q[I][1]-Q[I][JM1]
                        B2[I] = np.double(2.00)*(Q[I][1]+Q[I][NLAST])
                    #189
                    K1 = KR+JR-1
                    TCOS[K1+1] = np.double(-2.00)
                    K4 = K1+3-ISTAG
                    COSGEN(KR+ISTAG-2,1,np.double(0.00),FNUM,TCOS(K4))
                    K4 = K1+KR+1
                    COSGEN(JR-1,1,np.double(0.00),np.double(1.00),TCOS(K4))
                    MERGE(TCOS,K1,KR,K1+KR,JR-1,0)
                    COSGEN(KR,1,np.double(0.500),FDEN,TCOS(K1+1))
                    K2 = KR
                    K4 = K1+K2+1
                    COSGEN(LR,1,np.double(0.500),FDEN,TCOS(K4))
                    K3 = LR
                    K4 = 0
                    TRI3(MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
                    for I in range(1,MR+1):
                        B[I] = B[I]+B2[I]
                    #COTINUE
                    #190
                    TCOS[1] = np.double(2.00)
                    TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                    for I in range(1,MR+1):
                        Q[I][1] = Q[I][NLAST]
                    #191
                    # --->194
                #194
                J = NLAST - JR
                for I in range(1,MR+1):
                    B[I] = Q[I][NLAST]+Q[I][J]
                #195
                #CONTINUE
                #--->196
            if(LR != 0):#--->170
                #170
                for I in range(1,MR+1):
                    B[I] = np.double(0.500)*Q[I][1]-Q[I][JM1]+Q[I][J]
                #171
                if(NROD == 0):#---->173
                    for I in range(1,MR+1):
                        II = IP+I
                        B[I] = B[I]+P[II]
                        #172
                        #CONTINUE
                        #--->175
                else:
                    #173
                    for I in range(1,MR+1):
                        B[I] = B[I]+Q[I][NLAST]-Q[I][JM2]
                    #174
                    #CONTINUE
                    #--->175
                #175
                for I in range(1,MR+1):
                    T = np.double(0.500)*(Q[I][J]-Q[I][JM1]-Q[I][JP1])
                    Q[I][J] = T
                    B2[I] = Q[I][NLAST]+T
                    B3[I] = Q[I][1]+np.double(2.00)*T
                #176
                #177
                K1 = KR+2*JR-1
                K2 = KR+JR
                TCOS[K1+1] = np.double(-2.00)
                K4 = K1+3-ISTAG
                COSGEN(K2+ISTAG-2,1,np.double(0.00),FNUM,TCOS(K4))
                K4 = K1+K2+1
                COSGEN(JR-1,1,np.double(0.00),np.double(1.00),TCOS(K4))
                MERGE(TCOS,K1,K2,K1+K2,JR-1,0)
                K3 = K1+K2+LR
                COSGEN(JR,1,np.double(0.500),np.double(0.00),TCOS(K3+1))
                K4 = K3+JR+1
                COSGEN(KR,1,np.double(0.500),FDEN,TCOS(K4))
                MERGE(TCOS,K3,JR,K3+JR,KR,K1)
                if(LR != 0):#---->178
                    COSGEN(LR,1,np.double(0.500),FDEN,TCOS(K4))
                    MERGE(TCOS,K3,JR,K3+JR,LR,K3-LR)
                    COSGEN(KR,1,np.double(0.500),FDEN,TCOS(K4))
                #178
                K3 = KR
                K4 = KR
                TRI3(MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
                for I in range(1,MR+1):
                    B[I] = B[I]+B2[I]+B3[I]
                #179
                TCOS[1] = np.double(2.00)
                TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                for I in range(1,MR+1):
                    Q[I][J] = Q[I][J] +B[I]
                    B[I] = Q[I][1]+np.double(2.00)*Q[I][J]
                #180
                COSGEN(JR,1,np.double(0.500),np.double(0.00),TCOS)
                TRIX(JR,0,MR,A,BB,C,B,TCOS,D,W)
                if(JR == 1):#--->182
                    for I in range(1,MR+1):
                        Q[I][1] = B[I]
                    #181
                    #CONTINUE
                    #---->194
                else:
                    #182
                    for I in range(1,MR+1):
                        Q[I][1] = np.double(0.500)*Q[I][JM1]+B[I]
                    #183
                    #CONTINUE
                    #--->194
                #194
                J = NLAST - JR
                for I in range(1,MR+1):
                    B[I] = Q[I][NLAST]+Q[I][J]
                #195
                #CONTINUE
                #--->196
            if(N != 3):#--->161
                #161
                if(ISTAG == 1):#--->162
                    #162
                    for I in range(1,MR+1):
                        B[I] = Q[I][J]+np.double(0.500)*Q[I][1]-Q[I][JM1]+Q[I][NLAST]-Q[I][JM2]
                    #163
                    COSGEN(JR,1,np.double(0.500),np.double(0.00),TCOS)
                    TRIX(JR,0,MR,A,BB,C,B,TCOS,D,W)
                    for I in range(1,MR+1):
                        Q[I][J] = np.double(0.500)*(Q[I][J]-Q[I][JM1]-Q[I][JP1])+B[I]
                        B[I] = Q[I][1]+np.double(2.00)*Q[I][NLAST]+np.double(4.00)*Q[I][J]
                    #164
                    JR2 = 2*JR
                    COSGEN(JR,1,np.double(0.00),np.double(0.00),TCOS)
                    for I in range(1,JR):
                        I1 = JR+1
                        I2 = JR+1-I
                        TCOS[I1] = -TCOS[I2]
                    #165
                    #CONTINUE
                    TRIX(JR2,0,MR,A,BB,C,B,TCOS,D,W)
                    for I in range(1,MR+1):
                        Q[I][J] = Q[I][J]+B[I]
                        B[I] = Q[I][1]+np.double(2.00)*Q[I][J]
                    #166
                    COSGEN(KR,1,np.double(0.500),np.double(0.00),TCOS)
                    TRIX(KR,0,MR,A,BB,C,B,TCOS,D,W)
                    for I in range(1,MR+1):
                        Q[I][1] =np.double(0.500)*Q[I][1]-Q[I][JM1]+B[I]
                    #COTINUE
                    #167
                    #---->194
                if(ISTAG == 2):#--->170
                    #170
                    for I in range(1,MR+1):
                        B[I] = np.double(0.500)*Q[I][1]-Q[I][JM1]+Q[I][J]
                    #171
                    if(NROD == 0):#---->173
                        for I in range(1,MR+1):
                            II = IP+I
                            B[I] = B[I]+P[II]
                            #172
                            #CONTINUE
                            #--->175
                    else:
                        #173
                        for I in range(1,MR+1):
                            B[I] = B[I]+Q[I][NLAST]-Q[I][JM2]
                        #174
                        #CONTINUE
                        #--->175
                    #175
                    for I in range(1,MR+1):
                        T = np.double(0.500)*(Q[I][J]-Q[I][JM1]-Q[I][JP1])
                        Q[I][J] = T
                        B2[I] = Q[I][NLAST]+T
                        B3[I] = Q[I][1]+np.double(2.00)*T
                    #176
                    #177
                    K1 = KR+2*JR-1
                    K2 = KR+JR
                    TCOS[K1+1] = np.double(-2.00)
                    K4 = K1+3-ISTAG
                    COSGEN(K2+ISTAG-2,1,np.double(0.00),FNUM,TCOS(K4))
                    K4 = K1+K2+1
                    COSGEN(JR-1,1,np.double(0.00),np.double(1.00),TCOS(K4))
                    MERGE(TCOS,K1,K2,K1+K2,JR-1,0)
                    K3 = K1+K2+LR
                    COSGEN(JR,1,np.double(0.500),np.double(0.00),TCOS(K3+1))
                    K4 = K3+JR+1
                    COSGEN(KR,1,np.double(0.500),FDEN,TCOS(K4))
                    MERGE(TCOS,K3,JR,K3+JR,KR,K1)
                    if(LR != 0):#---->178
                        COSGEN(LR,1,np.double(0.500),FDEN,TCOS(K4))
                        MERGE(TCOS,K3,JR,K3+JR,LR,K3-LR)
                        COSGEN(KR,1,np.double(0.500),FDEN,TCOS(K4))
                    #178
                    K3 = KR
                    K4 = KR
                    TRI3(MR,A,BB,C,K,B,B2,B3,TCOS,D,W,W2,W3)
                    for I in range(1,MR+1):
                        B[I] = B[I]+B2[I]+B3[I]
                    #179
                    TCOS[1] = np.double(2.00)
                    TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                    for I in range(1,MR+1):
                        Q[I][J] = Q[I][J] +B[I]
                        B[I] = Q[I][1]+np.double(2.00)*Q[I][J]
                    #180
                    COSGEN(JR,1,np.double(0.500),np.double(0.00),TCOS)
                    TRIX(JR,0,MR,A,BB,C,B,TCOS,D,W)
                    if(JR == 1):#--->182
                        for I in range(1,MR+1):
                            Q[I][1] = B[I]
                        #181
                        #CONTINUE
                        #---->194
                    else:
                        #182
                        for I in range(1,MR+1):
                            Q[I][1] = np.double(0.500)*Q[I][JM1]+B[I]
                        #183
                        #CONTINUE
                        #--->194
                #194
                J = NLAST - JR
                for I in range(1,MR+1):
                    B[I] = Q[I][NLAST]+Q[I][J]
                #195
                #CONTINUE
                #--->196
            else:
                #156
                for I in range(1,MR+1):
                    B[I] = Q[I][2]
                #CONTINUE
                #157
                TCOS[1] = np.double(0.00)
                TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                for I in range(1,MR+1):
                    Q[I][2] = B[I]
                    B[I] = np.double(4.00)*B[I][1]+np.double(2.00)*Q[I][3]
                #158
                #CONTINUE
                TCOS[1] = np.double(-2.00)
                TCOS[2] = np.double(2.00)
                I1 = 0
                I2 = 0
                TRIX(I1,I2,MR,A,BB,C,B,TCOS,D,W)
                for I in range(1,MR+1):
                    Q[I][2] = Q[I][2] + B[I]
                    B[I] = Q[I][1] + np.double(2.00)*Q[I][2]
                #159
                TCOS[1] = np.double(0.00)
                TRIX(1,0,MR,A,BB,C,B,TCOS,D,W)
                for I in range(1,MR+1):
                    Q[I][1] = B[I]
                #160
                JR = 1
                I2R = 0
                #--->194
                #194
                J = NLAST - JR
                for I in range(1,MR+1):
                    B[I] = Q[I][NLAST]+Q[I][J]
                #195
                #CONTINUE
                #--->196
        #196
        JM2 = NLAST  - I2R
        if(JR == 1):#--->198
            for I in range(1,MR+1):
                Q[I][NLAST] = np.double(0.000)
            #197
            #CONTINUE
            #--->202
        else:
            #198
            if (NROD == 0):#--->200
                for I in range(1,MR+1):
                    II = IP+I
                    Q[I][NLAST] = P[II]
                #199
                IP = IP+MR
                #---->202
            else:
                #200
                for I in range(1,MR+1):
                    Q[I][NLAST] = Q[I][NLAST]-Q[I][JM2]
                #201
                #CONTINUE
                #--->202
        #202
        COSGEN (KR,1,np.double(0.500),FDEN,TCOS)
        COSGEN (LR,1,np.double(0.500),FDEN,TCOS(KR+1))
        if (LR == 0):#------>204
            for I in range(1,MR+1):
                B[I] = FISTAG*B[I]
        #203
        #CONTINUE
        #204
        TRIX (KR,LR,MR,A,BB,C,B,TCOS,D,W)
        for I in range(1,MR+1):
            Q[I][NLAST] = Q[I][NLAST]+B
        #CONTIUE
        #205
        NLASTP = NLAST
        #206
        JSTEP = JR
        JR = I2R
        I2R = I2R/2
        while(JR > 0):
            #
            if (JR != 0):#--->222
                if(MIXBND == 1):#---->207
                    #207
                    JSTART = 1+JR
                    #---->209
                if(MIXBND == 2):#---->208
                    #208
                    JSTART = JR
                #209
                KR = KR-JR
                if (NLAST+JR < N):#----->210
                    KR = KR-JR
                    NLAST = NLAST+JR
                    JSTOP = NLAST-JSTEP
                    #---->211
                else:
                    #210
                    JSTOP = NLAST-JR
                    #--->211
                #211
                LR = KR-JR
                COSGEN (JR,1,np.double(0.500),np.double(0.000),TCOS)
                for J in range(JSTART,JSTOP,JSTEP):
                    JM2 = J-JR
                    JP2 = J+JR
                    if(J == JR):#--->213
                        for I in range(1,MR+1):
                            B[I] = Q[I][J]+Q[I][JP2]
                        #212
                        #CONTINUE
                        #215
                    else:
                        #213
                        for I in range(1,MR+1):
                            B[I] = Q[I][J]+Q[I][JM2]+Q[I][JP2]
                        #214
                        #CONTINUE
                    #215
                    if(JR != 1):#---->217
                        for I in range(1,MR+1):
                            Q[I][J] = np.double(0.00)
                        #216
                        #CONTINUE
                        #---->219
                    else:
                        #217
                        JM1 = J-I2R
                        JP1 = J+I2R
                        for I in range(1,MR+1):
                            Q[I][J] = np.double(0.500)*(Q[I][J]-Q[I][JM1]-Q[I][JP1])
                        #218
                        #CONTINUE
                    #219
                    TRIX (JR,0,MR,A,BB,C,B,TCOS,D,W)
                    for I in range(1,MR+1):
                        Q[I][J] = Q[I][J]+B[I]
                        #220
                        #CONTINUE
                #221
                #CONTINUE
                NROD = 1
                if(NLAST+I2R > N):
                    NROD = 0
                if(NLASTP != NLAST):
                    #--->194
                    #194
                    J = NLAST - JR
                    for I in range(1,MR+1):
                        B[I] = Q[I][NLAST]+Q[I][J]
                    #195
                    #CONTINUE
                    #--->196
                    #196
                    JM2 = NLAST  - I2R
                    if(JR == 1):#--->198
                        for I in range(1,MR+1):
                            Q[I][NLAST] = np.double(0.000)
                        #197
                        #CONTINUE
                        #--->202
                    else:
                        #198
                        if (NROD == 0):#--->200
                            for I in range(1,MR+1):
                                II = IP+I
                                Q[I][NLAST] = P[II]
                            #199
                            IP = IP+MR
                            #---->202
                        else:
                            #200
                            for I in range(1,MR+1):
                                Q[I][NLAST] = Q[I][NLAST]-Q[I][JM2]
                            #201
                            #CONTINUE
                            #--->202
                    #202
                    COSGEN (KR,1,np.double(0.500),FDEN,TCOS)
                    COSGEN (LR,1,np.double(0.500),FDEN,TCOS(KR+1))
                    if (LR == 0):#------>204
                        for I in range(1,MR+1):
                            B[I] = FISTAG*B[I]
                    #203
                    #CONTINUE
                    #204
                    TRIX (KR,LR,MR,A,BB,C,B,TCOS,D,W)
                    for I in range(1,MR+1):
                        Q[I][NLAST] = Q[I][NLAST]+B
                    #CONTIUE
                    #205
                    NLASTP = NLAST
                    #206
                    JSTEP = JR
                    JR = I2R
                    I2R = I2R/2
                else:
                    #206
                    JSTEP = JR
                    JR = I2R
                    I2R = I2R/2
            else:
                break
        #222
        W[1] = IPSTOR
        return





def GENBUN (NPEROD,N,MPEROD,M,A,B,C,IDIMY,Y,IERROR,W):
    Y = np.zeros((IDIMY,1))
    W = np.zeros((1))
    B = np.zeros((1))
    A = np.zeros((1))
    C = np.zeros((1))

    IERROR = 0
    if(M <= 2):
        IERROR = 1
    if(N <= 2):
        IERROR = 2
    if(IDIMY < M):
        IERROR = 3
    if(NPEROD < 0 or NPEROD > 4):
        IERROR = 4
    if(MPEROD < 0 or MPEROD > 4):
        IERROR = 5
    if(MPEROD != 1):#---->102
        for I in range(2,M):
            if(A[I] != C[0]):#--->103
                #103
                IERROR = 6
            if(C[I] != C[0]):#---->103
                #103
                IERROR = 6
            if(B[I] != B[0]):#---->103
                #103
                IERROR = 6
        #101
        #---->104
        
    else:
        #102
        if (A[0] != 0 or C[M] != 0):
            IERROR = 7
        #--->104
    #104
    MP1 = M+1
    IWBA = MP1
    IWBB = IWBA+M
    IWBC = IWBB+M
    IWB2 = IWBC+M
    IWB3 = IWB2+M
    IWW1 = IWB3+M
    IWW2 = IWW1+M
    IWW3 = IWW2+M
    IWD = IWW3+M
    IWTCOS = IWD+M
    IWP = IWTCOS+4*N
    for I in range(1,M):
        K = IWBA+I-1
        W[K] = -A[I]
        K = IWBB+I-1
        W[K] = -C[I]
        K = IWBB+I-1
        W[K] = np.double(2.000)-B[I]
        for J in range(1,N):
            Y[I][J] = -Y[I][J]
        #105
        #CONTINUE
    #106
    MP = MPEROD+1
    NP = MPEROD+1
    if(MP == 1):#--->114
        #114
        MH = (M+1)/2
        MHM1 = MH-1
        MODD = 1
        if(MH*2 == M):
            MODD = 2
        for J in range(1,N):
            for I in range(1,MHM1):
                MHPI = MH+I
                MHMI = MH-I
                W[I] = Y[MHMI][J]-Y[MHPI][J]
                W[MHPI] = Y[MHM1][J]+Y [MHPI][J]
            #115
            #CONTINUE
            W[MH] = np.double(2.00)*Y[MH][J]
            if(MODD == 1):#--->117
                #117
                for I in range(1,M):
                    Y[I][J] = W[I]
                #118
                #CONTINUE
            if(MODD == 2):#--->116
                #116
                W[M] = np.doube(2.00)*Y[M][J]
                #117
                for I in range(1,M):
                    Y[I][J] = W[I]
                #118
                #CONTINUE
            #119
            K = IWBC+MHM1-1
            I = IWBA+MHM1
            W[K] = np.double(0.00)
            W[I] = np.double(0.00)
            W[K+1] = np.double(2.00)*W[K+1]
            if(MODD == 1):#--->120
                #120
                K = IWBB+MHM1-1
                W[K] = W[K] - W[I-1]
                W[IWBC-1] = W[IWBC-1]+W[IWBB-1]
                #---->122
                #122
                #--->107???
            if(MODD == 2):#---->121
                #121
                W[IWBB-1] = W[K+1]
                #122
                #--->107????
            #107
            if(NP == 1):#--->108
                #108
                POISP2 (M,N,W[IWBA],W[IWBB],W[IWBC],Y,IDIMY,W,W[IWB2],W[IWB3],W[IWW1],W[IWW2],W[IWW3],W[IWD],W[IWTCOS],W[IWP])
                #--->112
                while (IREV == 1 ):
                    #112
                    IPSTOR = W[IWW1]
                    IREV = 2
                    if(NPEROD != 4):#--->124
                        #113
                        if(MP == 1):#--->127
                            #127
                            for J in range(1,N):
                                for I in range(1,MHM1):
                                    MHMI = MH-I
                                    MHMI = MH+I
                                    W[MHMI] = np.double(0.500)*(Y[MHPI][J]+Y[I][J])
                                    W[MHPI] = np.double(0.500)*(Y[MHPI][J]-Y[I][J])
                                #128
                                #CONTINUE
                                W[MH] = np.double(0.500)*Y[MH][J]
                                if(MODD == 1):#--->130
                                    #130
                                    #CONTINUE
                                    for I in range(1,M):
                                        Y[i][J] = W[I]
                                    #131
                                    #CONTINUE
                                    break
                                if(MODD == 2):#---->129
                                    #129
                                    W[M] = np.float(0.500)*Y[M][J]
                                    #130
                                    for I in range(1,M):
                                        Y[i][J] = W[I]
                                    #131
                                    #CONTINUE
                                    break
                        if(MP == 2):#--->133
                            break
                            #end
                    else:
                        #124
                        for J in range(1,NBY2):
                            MSKIP = N+1-J
                            for I in range(1,M):
                                A1 = Y[I][J]
                                Y[I][J] = Y[I][MSKIP]
                                Y[I][MSKIP] = A1
                            #125
                        #126
                        if(IREV == 1):#--->110
                            #110
                            if(NP == 3):#--->110
                                POISN2 (M,N,1,2,W[IWBA],W[IWBB],W[IWBC],Y,IDIMY,W,W[IWB2],W[IWB3],W[IWW1],W[IWW2],W[IWW3],W[IWD],W[IWTCOS],W[IWP])
                                #---->112
                        if(IREV == 2):#--->113
                            #113
                            if(MP == 1):#--->127
                                #127
                                for J in range(1,N):
                                    for I in range(1,MHM1):
                                        MHMI = MH-I
                                        MHMI = MH+I
                                        W[MHMI] = np.double(0.500)*(Y[MHPI][J]+Y[I][J])
                                        W[MHPI] = np.double(0.500)*(Y[MHPI][J]-Y[I][J])
                                    #128
                                    #CONTINUE
                                    W[MH] = np.double(0.500)*Y[MH][J]
                                    if(MODD == 1):#--->130
                                        #130
                                        #CONTINUE
                                        for I in range(1,M):
                                            Y[i][J] = W[I]
                                        #131
                                        #CONTINUE
                                        break
                                    if(MODD == 2):#---->129
                                        #129
                                        W[M] = np.float(0.500)*Y[M][J]
                                        #130
                                        for I in range(1,M):
                                            Y[i][J] = W[I]
                                        #131
                                        #CONTINUE
                                        break
                            if(MP == 2):#--->133
                                break
                                #end

            if(NP == 2):#--->109
                #109
                POISD2 (M,N,1,W[IWBA],W[IWBB],W[IWBC],Y,IDIMY,W,W[IWW1],W[IWD],W[IWTCOS],W[IWP])
                #--->112
                while (IREV == 1 ):
                    #112
                    IPSTOR = W[IWW1]
                    IREV = 2
                    if(NPEROD != 4):#--->124
                        #113
                        if(MP == 1):#--->127
                            #127
                            for J in range(1,N):
                                for I in range(1,MHM1):
                                    MHMI = MH-I
                                    MHMI = MH+I
                                    W[MHMI] = np.double(0.500)*(Y[MHPI][J]+Y[I][J])
                                    W[MHPI] = np.double(0.500)*(Y[MHPI][J]-Y[I][J])
                                #128
                                #CONTINUE
                                W[MH] = np.double(0.500)*Y[MH][J]
                                if(MODD == 1):#--->130
                                    #130
                                    #CONTINUE
                                    for I in range(1,M):
                                        Y[i][J] = W[I]
                                    #131
                                    #CONTINUE
                                    break
                                if(MODD == 2):#---->129
                                    #129
                                    W[M] = np.float(0.500)*Y[M][J]
                                    #130
                                    for I in range(1,M):
                                        Y[i][J] = W[I]
                                    #131
                                    #CONTINUE
                                    break
                        if(MP == 2):#--->133
                            break
                            #end
                    else:
                        #124
                        for J in range(1,NBY2):
                            MSKIP = N+1-J
                            for I in range(1,M):
                                A1 = Y[I][J]
                                Y[I][J] = Y[I][MSKIP]
                                Y[I][MSKIP] = A1
                            #125
                        #126
                        if(IREV == 1):#--->110
                            #110
                            if(NP == 3):#--->110
                                POISN2 (M,N,1,2,W[IWBA],W[IWBB],W[IWBC],Y,IDIMY,W,W[IWB2],W[IWB3],W[IWW1],W[IWW2],W[IWW3],W[IWD],W[IWTCOS],W[IWP])
                                #---->112
                        if(IREV == 2):#--->113
                            #113
                            if(MP == 1):#--->127
                                #127
                                for J in range(1,N):
                                    for I in range(1,MHM1):
                                        MHMI = MH-I
                                        MHMI = MH+I
                                        W[MHMI] = np.double(0.500)*(Y[MHPI][J]+Y[I][J])
                                        W[MHPI] = np.double(0.500)*(Y[MHPI][J]-Y[I][J])
                                    #128
                                    #CONTINUE
                                    W[MH] = np.double(0.500)*Y[MH][J]
                                    if(MODD == 1):#--->130
                                        #130
                                        #CONTINUE
                                        for I in range(1,M):
                                            Y[i][J] = W[I]
                                        #131
                                        #CONTINUE
                                        break
                                    if(MODD == 2):#---->129
                                        #129
                                        W[M] = np.float(0.500)*Y[M][J]
                                        #130
                                        for I in range(1,M):
                                            Y[i][J] = W[I]
                                        #131
                                        #CONTINUE
                                        break
                            if(MP == 2):#--->133
                                break
                                #end
            if(NP == 3):#--->110
                POISN2 (M,N,1,2,W[IWBA],W[IWBB],W[IWBC],Y,IDIMY,W,W[IWB2],W[IWB3],W[IWW1],W[IWW2],W[IWW3],W[IWD],W[IWTCOS],W[IWP])
                #---->112
                while (IREV == 1 ):
                    #112
                    IPSTOR = W[IWW1]
                    IREV = 2
                    if(NPEROD != 4):#--->124
                        #113
                        if(MP == 1):#--->127
                            #127
                            for J in range(1,N):
                                for I in range(1,MHM1):
                                    MHMI = MH-I
                                    MHMI = MH+I
                                    W[MHMI] = np.double(0.500)*(Y[MHPI][J]+Y[I][J])
                                    W[MHPI] = np.double(0.500)*(Y[MHPI][J]-Y[I][J])
                                #128
                                #CONTINUE
                                W[MH] = np.double(0.500)*Y[MH][J]
                                if(MODD == 1):#--->130
                                    #130
                                    #CONTINUE
                                    for I in range(1,M):
                                        Y[i][J] = W[I]
                                    #131
                                    #CONTINUE
                                    break
                                if(MODD == 2):#---->129
                                    #129
                                    W[M] = np.float(0.500)*Y[M][J]
                                    #130
                                    for I in range(1,M):
                                        Y[i][J] = W[I]
                                    #131
                                    #CONTINUE
                                    break
                        if(MP == 2):#--->133
                            break
                            #end
                    else:
                        #124
                        for J in range(1,NBY2):
                            MSKIP = N+1-J
                            for I in range(1,M):
                                A1 = Y[I][J]
                                Y[I][J] = Y[I][MSKIP]
                                Y[I][MSKIP] = A1
                            #125
                        #126
                        if(IREV == 1):#--->110
                            #110
                            if(NP == 3):#--->110
                                POISN2 (M,N,1,2,W[IWBA],W[IWBB],W[IWBC],Y,IDIMY,W,W[IWB2],W[IWB3],W[IWW1],W[IWW2],W[IWW3],W[IWD],W[IWTCOS],W[IWP])
                                #---->112
                        if(IREV == 2):#--->113
                            #113
                            if(MP == 1):#--->127
                                #127
                                for J in range(1,N):
                                    for I in range(1,MHM1):
                                        MHMI = MH-I
                                        MHMI = MH+I
                                        W[MHMI] = np.double(0.500)*(Y[MHPI][J]+Y[I][J])
                                        W[MHPI] = np.double(0.500)*(Y[MHPI][J]-Y[I][J])
                                    #128
                                    #CONTINUE
                                    W[MH] = np.double(0.500)*Y[MH][J]
                                    if(MODD == 1):#--->130
                                        #130
                                        #CONTINUE
                                        for I in range(1,M):
                                            Y[i][J] = W[I]
                                        #131
                                        #CONTINUE
                                        break
                                    if(MODD == 2):#---->129
                                        #129
                                        W[M] = np.float(0.500)*Y[M][J]
                                        #130
                                        for I in range(1,M):
                                            Y[i][J] = W[I]
                                        #131
                                        #CONTINUE
                                        break
                            if(MP == 2):#--->133
                                break
                                #end
            if(NP == 4):#--->111
                POISN2 (M,N,1,1,W[IWBA],W[IWBB],W[IWBC],Y,IDIMY,W,W[IWB2],W[IWB3],W[IWW1],W[IWW2],W[IWW3],W[IWD],W[IWTCOS],W[IWP])
                #--->112
                while (IREV == 1 ):
                    #112
                    IPSTOR = W[IWW1]
                    IREV = 2
                    if(NPEROD != 4):#--->124
                        #113
                        if(MP == 1):#--->127
                            #127
                            for J in range(1,N):
                                for I in range(1,MHM1):
                                    MHMI = MH-I
                                    MHMI = MH+I
                                    W[MHMI] = np.double(0.500)*(Y[MHPI][J]+Y[I][J])
                                    W[MHPI] = np.double(0.500)*(Y[MHPI][J]-Y[I][J])
                                #128
                                #CONTINUE
                                W[MH] = np.double(0.500)*Y[MH][J]
                                if(MODD == 1):#--->130
                                    #130
                                    #CONTINUE
                                    for I in range(1,M):
                                        Y[i][J] = W[I]
                                    #131
                                    #CONTINUE
                                    break
                                if(MODD == 2):#---->129
                                    #129
                                    W[M] = np.float(0.500)*Y[M][J]
                                    #130
                                    for I in range(1,M):
                                        Y[i][J] = W[I]
                                    #131
                                    #CONTINUE
                                    break
                        if(MP == 2):#--->133
                            break
                            #end
                    else:
                        #124
                        for J in range(1,NBY2):
                            MSKIP = N+1-J
                            for I in range(1,M):
                                A1 = Y[I][J]
                                Y[I][J] = Y[I][MSKIP]
                                Y[I][MSKIP] = A1
                            #125
                        #126
                        if(IREV == 1):#--->110
                            #110
                            if(NP == 3):#--->110
                                POISN2 (M,N,1,2,W[IWBA],W[IWBB],W[IWBC],Y,IDIMY,W,W[IWB2],W[IWB3],W[IWW1],W[IWW2],W[IWW3],W[IWD],W[IWTCOS],W[IWP])
                                #---->112
                        if(IREV == 2):#--->113
                            #113
                            if(MP == 1):#--->127
                                #127
                                for J in range(1,N):
                                    for I in range(1,MHM1):
                                        MHMI = MH-I
                                        MHMI = MH+I
                                        W[MHMI] = np.double(0.500)*(Y[MHPI][J]+Y[I][J])
                                        W[MHPI] = np.double(0.500)*(Y[MHPI][J]-Y[I][J])
                                    #128
                                    #CONTINUE
                                    W[MH] = np.double(0.500)*Y[MH][J]
                                    if(MODD == 1):#--->130
                                        #130
                                        #CONTINUE
                                        for I in range(1,M):
                                            Y[i][J] = W[I]
                                        #131
                                        #CONTINUE
                                        break
                                    if(MODD == 2):#---->129
                                        #129
                                        W[M] = np.float(0.500)*Y[M][J]
                                        #130
                                        for I in range(1,M):
                                            Y[i][J] = W[I]
                                        #131
                                        #CONTINUE
                                        break
                            if(MP == 2):#--->133
                                break
                                #end
            if(NP == 5):#---->123 
                #123
                IREV = 1
                NBY2 = N/2
                #124
                for J in range(1,NBY2):
                    MSKIP = N+1-J
                    for I in range(1,M):
                        A1 = Y[I][J]
                        Y[I][J] = Y[I][MSKIP]
                        Y[I][MSKIP] = A1
                    #125
                    #CONTINUE
                #126
                #CONTINUE
                while(IREV == 1):
                    if(IREV == 1):#--->110
                        #110
                        if(NP == 3):#--->110
                            POISN2 (M,N,1,2,W[IWBA],W[IWBB],W[IWBC],Y,IDIMY,W,W[IWB2],W[IWB3],W[IWW1],W[IWW2],W[IWW3],W[IWD],W[IWTCOS],W[IWP])
                            #---->112
                        #112
                        IPSTOR = W[IWW1]
                        IREV = 2
                        if(NPEROD == 4):#--->124
                            #124
                            for J in range(1,NBY2):
                                MSKIP = N+1-J
                            for I in range(1,M):
                                A1 = Y[I][J]
                                Y[I][J] = Y[I][MSKIP]
                                Y[I][MSKIP] = A1
                            #125
                            #CONTINUE
                        else:
                            #113
                            if(MP == 1):#--->127
                                #127
                                for J in range(1,N):
                                    for I in range(1,MHM1):
                                        MHMI = MH-I
                                        MHMI = MH+I
                                        W[MHMI] = np.double(0.500)*(Y[MHPI][J]+Y[I][J])
                                        W[MHPI] = np.double(0.500)*(Y[MHPI][J]-Y[I][J])
                                    #128
                                #CONTINUE
                                        W[MH] = np.double(0.500)*Y[MH][J]
                                    if(MODD == 1):#--->130
                                        #130
                                        #CONTINUE
                                        for I in range(1,M):
                                            Y[i][J] = W[I]
                                        #131
                                        #CONTINUE
                                        break
                                    if(MODD == 2):#---->129
                                        #129
                                        W[M] = np.float(0.500)*Y[M][J]
                                        #130
                                        for I in range(1,M):
                                            Y[i][J] = W[I]
                                        #131
                                        #CONTINUE
                                        break
                            if(MP == 2):#--->133
                                break
                                #end

                    if(IREV == 2):#--->113
                        #113
                        if(MP == 1):#--->127
                            #127
                            for J in range(1,N):
                                for I in range(1,MHM1):
                                    MHMI = MH-I
                                    MHMI = MH+I
                                    W[MHMI] = np.double(0.500)*(Y[MHPI][J]+Y[I][J])
                                    W[MHPI] = np.double(0.500)*(Y[MHPI][J]-Y[I][J])
                                #128
                            #CONTINUE
                                    W[MH] = np.double(0.500)*Y[MH][J]
                                if(MODD == 1):#--->130
                                    #130
                                    #CONTINUE
                                    for I in range(1,M):
                                        Y[i][J] = W[I]
                                    #131
                                    #CONTINUE
                                    break
                                if(MODD == 2):#---->129
                                    #129
                                    W[M] = np.float(0.500)*Y[M][J]
                                    #130
                                    for I in range(1,M):
                                        Y[i][J] = W[I]
                                    #131
                                    #CONTINUE
                                    break
                        if(MP == 2):#--->133
                            break
                            #end











    if(MP == 2):#--->107
        #107
        if(NP == 1):#--->108
            #108
            POISP2 (M,N,W[IWBA],W[IWBB],W[IWBC],Y,IDIMY,W,W[IWB2],W[IWB3],W[IWW1],W[IWW2],W[IWW3],W[IWD],W[IWTCOS],W[IWP])
            #--->112
            while (IREV == 1 ):
                #112
                IPSTOR = W[IWW1]
                IREV = 2
                if(NPEROD != 4):#--->124
                    #113
                    if(MP == 1):#--->127
                        #127
                        for J in range(1,N):
                            for I in range(1,MHM1):
                                MHMI = MH-I
                                MHMI = MH+I
                                W[MHMI] = np.double(0.500)*(Y[MHPI][J]+Y[I][J])
                                W[MHPI] = np.double(0.500)*(Y[MHPI][J]-Y[I][J])
                            #128
                            #CONTINUE
                            W[MH] = np.double(0.500)*Y[MH][J]
                            if(MODD == 1):#--->130
                                #130
                                #CONTINUE
                                for I in range(1,M):
                                    Y[i][J] = W[I]
                                #131
                                #CONTINUE
                                break
                            if(MODD == 2):#---->129
                                #129
                                W[M] = np.float(0.500)*Y[M][J]
                                #130
                                for I in range(1,M):
                                    Y[i][J] = W[I]
                                #131
                                #CONTINUE
                                break
                    if(MP == 2):#--->133
                        break
                        #end
                else:
                    #124
                    for J in range(1,NBY2):
                        MSKIP = N+1-J
                        for I in range(1,M):
                            A1 = Y[I][J]
                            Y[I][J] = Y[I][MSKIP]
                            Y[I][MSKIP] = A1
                        #125
                    #126
                    if(IREV == 1):#--->110
                        #110
                        if(NP == 3):#--->110
                            POISN2 (M,N,1,2,W[IWBA],W[IWBB],W[IWBC],Y,IDIMY,W,W[IWB2],W[IWB3],W[IWW1],W[IWW2],W[IWW3],W[IWD],W[IWTCOS],W[IWP])
                            #---->112
                    if(IREV == 2):#--->113
                        #113
                        if(MP == 1):#--->127
                            #127
                            for J in range(1,N):
                                for I in range(1,MHM1):
                                    MHMI = MH-I
                                    MHMI = MH+I
                                    W[MHMI] = np.double(0.500)*(Y[MHPI][J]+Y[I][J])
                                    W[MHPI] = np.double(0.500)*(Y[MHPI][J]-Y[I][J])
                                #128
                                #CONTINUE
                                W[MH] = np.double(0.500)*Y[MH][J]
                                if(MODD == 1):#--->130
                                    #130
                                    #CONTINUE
                                    for I in range(1,M):
                                        Y[i][J] = W[I]
                                    #131
                                    #CONTINUE
                                    break
                                if(MODD == 2):#---->129
                                    #129
                                    W[M] = np.float(0.500)*Y[M][J]
                                    #130
                                    for I in range(1,M):
                                        Y[i][J] = W[I]
                                    #131
                                    #CONTINUE
                                    break
                        if(MP == 2):#--->133
                            break
                            #end

        if(NP == 2):#--->109
            #109
            POISD2 (M,N,1,W[IWBA],W[IWBB],W[IWBC],Y,IDIMY,W,W[IWW1],W[IWD],W[IWTCOS],W[IWP])
            #--->112
            while (IREV == 1 ):
                #112
                IPSTOR = W[IWW1]
                IREV = 2
                if(NPEROD != 4):#--->124
                    #113
                    if(MP == 1):#--->127
                        #127
                        for J in range(1,N):
                            for I in range(1,MHM1):
                                MHMI = MH-I
                                MHMI = MH+I
                                W[MHMI] = np.double(0.500)*(Y[MHPI][J]+Y[I][J])
                                W[MHPI] = np.double(0.500)*(Y[MHPI][J]-Y[I][J])
                            #128
                            #CONTINUE
                            W[MH] = np.double(0.500)*Y[MH][J]
                            if(MODD == 1):#--->130
                                #130
                                #CONTINUE
                                for I in range(1,M):
                                    Y[i][J] = W[I]
                                #131
                                #CONTINUE
                                break
                            if(MODD == 2):#---->129
                                #129
                                W[M] = np.float(0.500)*Y[M][J]
                                #130
                                for I in range(1,M):
                                    Y[i][J] = W[I]
                                #131
                                #CONTINUE
                                break
                    if(MP == 2):#--->133
                        break
                        #end
                else:
                    #124
                    for J in range(1,NBY2):
                        MSKIP = N+1-J
                        for I in range(1,M):
                            A1 = Y[I][J]
                            Y[I][J] = Y[I][MSKIP]
                            Y[I][MSKIP] = A1
                        #125
                    #126
                    if(IREV == 1):#--->110
                        #110
                        if(NP == 3):#--->110
                            POISN2 (M,N,1,2,W[IWBA],W[IWBB],W[IWBC],Y,IDIMY,W,W[IWB2],W[IWB3],W[IWW1],W[IWW2],W[IWW3],W[IWD],W[IWTCOS],W[IWP])
                            #---->112
                    if(IREV == 2):#--->113
                        #113
                        if(MP == 1):#--->127
                            #127
                            for J in range(1,N):
                                for I in range(1,MHM1):
                                    MHMI = MH-I
                                    MHMI = MH+I
                                    W[MHMI] = np.double(0.500)*(Y[MHPI][J]+Y[I][J])
                                    W[MHPI] = np.double(0.500)*(Y[MHPI][J]-Y[I][J])
                                #128
                                #CONTINUE
                                W[MH] = np.double(0.500)*Y[MH][J]
                                if(MODD == 1):#--->130
                                    #130
                                    #CONTINUE
                                    for I in range(1,M):
                                        Y[i][J] = W[I]
                                    #131
                                    #CONTINUE
                                    break
                                if(MODD == 2):#---->129
                                    #129
                                    W[M] = np.float(0.500)*Y[M][J]
                                    #130
                                    for I in range(1,M):
                                        Y[i][J] = W[I]
                                    #131
                                    #CONTINUE
                                    break
                        if(MP == 2):#--->133
                            break
                            #end
        if(NP == 3):#--->110
            POISN2 (M,N,1,2,W[IWBA],W[IWBB],W[IWBC],Y,IDIMY,W,W[IWB2],W[IWB3],W[IWW1],W[IWW2],W[IWW3],W[IWD],W[IWTCOS],W[IWP])
            #---->112
            while (IREV == 1 ):
                #112
                IPSTOR = W[IWW1]
                IREV = 2
                if(NPEROD != 4):#--->124
                    #113
                    if(MP == 1):#--->127
                        #127
                        for J in range(1,N):
                            for I in range(1,MHM1):
                                MHMI = MH-I
                                MHMI = MH+I
                                W[MHMI] = np.double(0.500)*(Y[MHPI][J]+Y[I][J])
                                W[MHPI] = np.double(0.500)*(Y[MHPI][J]-Y[I][J])
                            #128
                            #CONTINUE
                            W[MH] = np.double(0.500)*Y[MH][J]
                            if(MODD == 1):#--->130
                                #130
                                #CONTINUE
                                for I in range(1,M):
                                    Y[i][J] = W[I]
                                #131
                                #CONTINUE
                                break
                            if(MODD == 2):#---->129
                                #129
                                W[M] = np.float(0.500)*Y[M][J]
                                #130
                                for I in range(1,M):
                                    Y[i][J] = W[I]
                                #131
                                #CONTINUE
                                break
                    if(MP == 2):#--->133
                        break
                        #end
                else:
                    #124
                    for J in range(1,NBY2):
                        MSKIP = N+1-J
                        for I in range(1,M):
                            A1 = Y[I][J]
                            Y[I][J] = Y[I][MSKIP]
                            Y[I][MSKIP] = A1
                        #125
                    #126
                    if(IREV == 1):#--->110
                        #110
                        if(NP == 3):#--->110
                            POISN2 (M,N,1,2,W[IWBA],W[IWBB],W[IWBC],Y,IDIMY,W,W[IWB2],W[IWB3],W[IWW1],W[IWW2],W[IWW3],W[IWD],W[IWTCOS],W[IWP])
                            #---->112
                    if(IREV == 2):#--->113
                        #113
                        if(MP == 1):#--->127
                            #127
                            for J in range(1,N):
                                for I in range(1,MHM1):
                                    MHMI = MH-I
                                    MHMI = MH+I
                                    W[MHMI] = np.double(0.500)*(Y[MHPI][J]+Y[I][J])
                                    W[MHPI] = np.double(0.500)*(Y[MHPI][J]-Y[I][J])
                                #128
                                #CONTINUE
                                W[MH] = np.double(0.500)*Y[MH][J]
                                if(MODD == 1):#--->130
                                    #130
                                    #CONTINUE
                                    for I in range(1,M):
                                        Y[i][J] = W[I]
                                    #131
                                    #CONTINUE
                                    break
                                if(MODD == 2):#---->129
                                    #129
                                    W[M] = np.float(0.500)*Y[M][J]
                                    #130
                                    for I in range(1,M):
                                        Y[i][J] = W[I]
                                    #131
                                    #CONTINUE
                                    break
                        if(MP == 2):#--->133
                            break
                            #end
        if(NP == 4):#--->111
            POISN2 (M,N,1,1,W[IWBA],W[IWBB],W[IWBC],Y,IDIMY,W,W[IWB2],W[IWB3],W[IWW1],W[IWW2],W[IWW3],W[IWD],W[IWTCOS],W[IWP])
            #--->112
            while (IREV == 1 ):
                #112
                IPSTOR = W[IWW1]
                IREV = 2
                if(NPEROD != 4):#--->124
                    #113
                    if(MP == 1):#--->127
                        #127
                        for J in range(1,N):
                            for I in range(1,MHM1):
                                MHMI = MH-I
                                MHMI = MH+I
                                W[MHMI] = np.double(0.500)*(Y[MHPI][J]+Y[I][J])
                                W[MHPI] = np.double(0.500)*(Y[MHPI][J]-Y[I][J])
                            #128
                            #CONTINUE
                            W[MH] = np.double(0.500)*Y[MH][J]
                            if(MODD == 1):#--->130
                                #130
                                #CONTINUE
                                for I in range(1,M):
                                    Y[i][J] = W[I]
                                #131
                                #CONTINUE
                                break
                            if(MODD == 2):#---->129
                                #129
                                W[M] = np.float(0.500)*Y[M][J]
                                #130
                                for I in range(1,M):
                                    Y[i][J] = W[I]
                                #131
                                #CONTINUE
                                break
                    if(MP == 2):#--->133
                        break
                        #end
                else:
                    #124
                    for J in range(1,NBY2):
                        MSKIP = N+1-J
                        for I in range(1,M):
                            A1 = Y[I][J]
                            Y[I][J] = Y[I][MSKIP]
                            Y[I][MSKIP] = A1
                        #125
                    #126
                    if(IREV == 1):#--->110
                        #110
                        if(NP == 3):#--->110
                            POISN2 (M,N,1,2,W[IWBA],W[IWBB],W[IWBC],Y,IDIMY,W,W[IWB2],W[IWB3],W[IWW1],W[IWW2],W[IWW3],W[IWD],W[IWTCOS],W[IWP])
                            #---->112
                    if(IREV == 2):#--->113
                        #113
                        if(MP == 1):#--->127
                            #127
                            for J in range(1,N):
                                for I in range(1,MHM1):
                                    MHMI = MH-I
                                    MHMI = MH+I
                                    W[MHMI] = np.double(0.500)*(Y[MHPI][J]+Y[I][J])
                                    W[MHPI] = np.double(0.500)*(Y[MHPI][J]-Y[I][J])
                                #128
                                #CONTINUE
                                W[MH] = np.double(0.500)*Y[MH][J]
                                if(MODD == 1):#--->130
                                    #130
                                    #CONTINUE
                                    for I in range(1,M):
                                        Y[i][J] = W[I]
                                    #131
                                    #CONTINUE
                                    break
                                if(MODD == 2):#---->129
                                    #129
                                    W[M] = np.float(0.500)*Y[M][J]
                                    #130
                                    for I in range(1,M):
                                        Y[i][J] = W[I]
                                    #131
                                    #CONTINUE
                                    break
                        if(MP == 2):#--->133
                            break
                            #end
        if(NP == 5):#---->123 
            #123
            IREV = 1
            NBY2 = N/2
            #124
            for J in range(1,NBY2):
                MSKIP = N+1-J
                for I in range(1,M):
                    A1 = Y[I][J]
                    Y[I][J] = Y[I][MSKIP]
                    Y[I][MSKIP] = A1
                #125
                #CONTINUE
            #126
            #CONTINUE
            while(IREV == 1):
                if(IREV == 1):#--->110
                    #110
                    if(NP == 3):#--->110
                        POISN2 (M,N,1,2,W[IWBA],W[IWBB],W[IWBC],Y,IDIMY,W,W[IWB2],W[IWB3],W[IWW1],W[IWW2],W[IWW3],W[IWD],W[IWTCOS],W[IWP])
                        #---->112
                    #112
                    IPSTOR = W[IWW1]
                    IREV = 2
                    if(NPEROD == 4):#--->124
                        #124
                        for J in range(1,NBY2):
                            MSKIP = N+1-J
                        for I in range(1,M):
                            A1 = Y[I][J]
                            Y[I][J] = Y[I][MSKIP]
                            Y[I][MSKIP] = A1
                        #125
                        #CONTINUE
                    else:
                        #113
                        if(MP == 1):#--->127
                            #127
                            for J in range(1,N):
                                for I in range(1,MHM1):
                                    MHMI = MH-I
                                    MHMI = MH+I
                                    W[MHMI] = np.double(0.500)*(Y[MHPI][J]+Y[I][J])
                                    W[MHPI] = np.double(0.500)*(Y[MHPI][J]-Y[I][J])
                                #128
                            #CONTINUE
                                    W[MH] = np.double(0.500)*Y[MH][J]
                                if(MODD == 1):#--->130
                                    #130
                                    #CONTINUE
                                    for I in range(1,M):
                                        Y[i][J] = W[I]
                                    #131
                                    #CONTINUE
                                    break
                                if(MODD == 2):#---->129
                                    #129
                                    W[M] = np.float(0.500)*Y[M][J]
                                    #130
                                    for I in range(1,M):
                                        Y[i][J] = W[I]
                                    #131
                                    #CONTINUE
                                    break
                        if(MP == 2):#--->133
                            break
                            #end

                if(IREV == 2):#--->113
                    #113
                    if(MP == 1):#--->127
                        #127
                        for J in range(1,N):
                            for I in range(1,MHM1):
                                MHMI = MH-I
                                MHMI = MH+I
                                W[MHMI] = np.double(0.500)*(Y[MHPI][J]+Y[I][J])
                                W[MHPI] = np.double(0.500)*(Y[MHPI][J]-Y[I][J])
                            #128
                        #CONTINUE
                                W[MH] = np.double(0.500)*Y[MH][J]
                            if(MODD == 1):#--->130
                                #130
                                #CONTINUE
                                for I in range(1,M):
                                    Y[i][J] = W[I]
                                #131
                                #CONTINUE
                                break
                            if(MODD == 2):#---->129
                                #129
                                W[M] = np.float(0.500)*Y[M][J]
                                #130
                                for I in range(1,M):
                                    Y[i][J] = W[I]
                                #131
                                #CONTINUE
                                break
                    if(MP == 2):#--->133
                        break
                        #end
    
    #132
    #133
    #CONTINUE
    W[0] = IPSTOR+IWP-1




import math

def COSGEN(BN,JUMP,FNUM,FEDN,A):
    A = np.zeros((1))
    

    PI = PIMACH(DUM)
    if(N != 0):
        if(IJUMP == 1):
            NP1 = N+1
            Y = PI /(np.double(N) + FDEN)
            for I in range(1,N):
                X = np.double(NP1 -I)-FNUM
                A[I] = np.double(2.00)*math.cos(X*Y)
        else:
            K3 = N/IJUMP+1
            K4 = K3-1
            PIBYN = PI/np.double(N+IJUMP)
            for K in range(1,IJUMP):
                K1 = (K-1)*K3
                K5 = (K-1)*K4
                for I in range(1,K4):
                    X = K1+I
                    K2 = K5+I
                    A[K2] = np.double(2.00)*math.cos(X*PIBYN)

    return



#ni idea de como se declara:
#SUBROUTINE MINSO4 (USOL, IDMN, ZN, ZM, PERTB)

def MINSO4 (USOL, IDMN, ZN, ZM, PERTB):
    

    ISTR = 1
    IFNL = K
    JSTR = 1
    JFNL = L


    UTE = 0.0
    ETE = 0.0
    for I in range (IS,MS): 
        II = I-IS+1
        for J in range (JS,NS):

            JJ= J-JS+1
            ETE =  ETE+ZM[II]*ZN[JJ]
            UTE = UTE+USOL[I][J]+ZM[II]*ZN[JJ]

    PERTRB = UTE/ETE

    for I in range (ISTR,IFNL):
        for J in range (JSTR,JFNL):
            USOL[I][J] = USOL[I][J]-PERTRB

    return


def POISP2(M,N,A,BB,C,Q,IDIMQ,B,B2,B3,W,W2,W3,D,TCOS,P):
    A = np.zeros((1))
    BB = np.zeros((1))
    C = np.zeros((1))
    Q = np.zeros((IDIMQ,1))
    B = np.zeros((1))
    B2 = np.zeros((1))
    B3 = np.zeros((1))
    W = np.zeros((1))
    W2 = np.zeros((1))
    W3 = np.zeros((1))
    D = np.zeros((1))
    TCOS = np.zeros((1))
    P = np.zeros((1))


    MR = M
    NR = (N+1)/2
    NRM1 = NR-1
    if(2*NR == N):#--->107
        for J in range(1,NRM1):
            NRMJ = NR-J
            NRPJ = NR+J
            for I in range(1,MR)+1:
                S = Q[I][NRMJ]-Q[I][NRPJ]
                T = Q[I][NRMJ]+Q[I][NRPJ]
                Q[I][NRMJ] = S
                Q[I][NRPJ] = T
                #101
                #CONTINUE
        #102
        for I in range(1,MR+1):
            Q[I][NR] = np.double(2.00)*Q[I][NR]
            Q[I][N] = np.double(2.00)*Q[I][N]
        #103
        #CONTINUE
        POISD2(MR,NRM1,1,A,BB,C,Q,IDIMQ,B,W,D,TCOS,P)
        IPSTOR = W[1]
        POISN2(MR,NR+1,1,1,A,BB,C,Q(1,NR),IDIMQ,B,B2,B3,W,W2,W3,D,1,TCOS,P)
        IPSTOR = MAX0(IPSTOR,int(W[1]))
        for J in range(1,NRM1):
            NRMJ = NR-J
            NRPJ = NR+J
            for I in range(1,MR+1):
                S = np.double(0.500)*(Q[I][NRPJ]+Q[I][NRMJ])
                T = np.double(0.500)*(Q[I][NRPJ]-Q[I][NRMJ])
            #CONTINUE
            #104
        #105
        for I in range(1,MR+1):
            Q[I][NR] = np.double(0.500)*Q[I][NR]
            Q[I][N] = np.double(0.500)*Q[I][N]
        #106
        #--->118
    else:
        #107
        for J in range(1,NRM1):
            NRPJ = N+1-J
            for I in range(1,MR):
                S = Q[I][J]-Q[I][NRPJ]
                T = Q[I][J]+Q[I][NRPJ]
                Q[I][J] = S
                Q[I][NRPJ] = T
                #CONTINUE
                #108
        #109
        for I in range(1,MR+1):
            Q[I][NR] = np.double(2.00)*Q[I][NR]
        #110
        LH = NRM1/2
        for J in range(1,LH):
            NRMJ = NR-J
            for I in range(1,MR+1):
                S = Q[I][J]
                Q[I][J] = Q[I][NRMJ]
                Q[I][NRMJ] = S
            #COTINUE
            #111
        #112
        POISD2 (MR,NRM1,2,A,BB,C,Q,IDIMQ,B,W,D,TCOS,P)
        IPSTOR = W[1]
        POISN2 (MR,NR,2,1,A,BB,C,Q(1,NR),IDIMQ,B,B2,B3,W,W2,W3,D,1,TCOS,P)
        IPSTOR = MAX0(IPSTOR,int(W[1]))
        for J in range(1,NRM1):
            NRMJ = NR+J
            for I in range(1,MR+1):
                S = np.double(0.500)*(Q[I][NRPJ]+Q[I][J])
                T = np.double(0.500)*(Q[I][NRPJ]-Q[I][J])
                Q[I][NRPJ] = T
                Q[I][J] = S
            #CONTINUE
            #113
        #114
        #CONTINUE
        for I in range(1,MR+1):
            Q[I][NR] = np.double(0.500)*Q[I][NR]
        #115
        #CONTINUE
        for I in range(1,LH):
            NRMJ = NR -J
            for I in range(1,MR+1):
                S = Q[I][J]
                Q[I][J] = Q[I][NRMJ]
                Q[I][NRMJ] = S
            #116
            #CONTINUE
        #117
        #COTINUE
    #118
    #CONTINUE
    W[1] = IPSTOR
    return









                        
            


                




                                

                             
                    

                            
                    
                        



                            

                            







                    

                        
                           










    







            







