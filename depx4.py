import numpy as np
def DEPX4(IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,D,N,NBDCND,BDC,BDD,COFX,GRHS,USOL,IDMN,W,PERTRB,IERROR):
    W=[]
    CHKPR4(IORDER,A,B,M,MBDCND,C,D,N,NBDCND,COFX,IDMN,IERROR)
    if IERROR != 0:
        return
    L=N+1
    if NBDCND == 0:
        L=N
    K=M+1
    L=N+1
#   ESTIMATE LOG BASE 2 OF N
    LOG2N=int(np.log(np.float32(N+1))/np.log(np.float32(2.0))+np.float32(0.50))
    LENGTH=4*(N+1)+(10+LOG2N)*(M+1)
    IERROR = 11
    LINPUT = int(W[0]+0.5)
    LOUTPT = LENGTH+6*(K+L)+1
    W.append(np.float32(LOUTPT))
#   IF (LOUTPT .GT. LINPUT) RETURN
    IERROR = 0
#   SET WORK SPACE INDICES
    I1 = LENGTH+2
    I2 = I1+L
    I3 = I2+L
    I4 = I3+L
    I5 = I4+L
    I6 = I5+L
    I7 = I6+L
    I8 = I7+K
    I9 = I8+K
    I10 = I9+K
    I11 = I10+K
    I12 = I11+K
    I13 = 2
    SPELIA4(IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,D,N,NBDCND,BDC,BDD,COFX,W(I1),W(I2),W(I3),W(I4),W(I5),W(I6),W(I7),W(I8),W(I9),W(I10),W(I11),W(I12),GRHS,USOL,IDMN,W(I13),PERTRB,IERROR)
    return 
def SPELI4(IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,D,N,NBDCND,BDC,BDD,COFX,AN,BN,CN,DN,UN,ZN,AM,BM,CM,DM,UM,ZM,GRHS,USOL,IDMN,W,PERTRB,IERROR):
    global KSWX,KSWY,K,L
    global AIT,BIT,CIT,DIT
    global MIT,NIT,IS,MS
    global JS,NS,DLX,DLY
    global TDLX3,TDLY3,DLX4,DLY4
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
    if (KSWX==2 or KSWX==3):
        f=0
def CHKPR4(IORDER,A,B,M,MBDCND,C,D,N,NBDCND,COFX,IDMN,IERROR):
    IERROR=1
    if (A>=B and C>=D):
        return
    IERROR=2
    if (MBDCND<0 or MBDCND>4):
        return
    IERROR=3
    if (NBDCND<0 or NBDCND>4):
        return
    IERROR=5
    if (IDMN<7):
        return
    IERROR=6
    if (M>(IDMN-1)or M<6):
        return
    IERROR=7
    if (N<5):
        return
    IERROR=8
    if (IORDER!=2 and IORDER!=4):
        return
    DLX=(B-A)/float(M)
    for I in range(2,M):
        XI=A+float(I-1)*DLX
        COFX(XI,AI,BI,CI)
    if (AI<float(0.00)):
        IERROR=10
    IERROR=0
def CHKSN4(MBDCND,NBDCND,ALPHA,BETA,COFX,SINGLR):
    SINGLR=False
    if ((MBDCND!=0 and MBDCND!=3)or(NBDCND!=0 and NBDCND!=3)):
        return
    if (MBDCND==3):
        if (ALPHA!=float(0.00)or BETA!=float(0.00)):
            return
    for I in range(IS,MS):
        XI=AIT+(I-1)*DLX
        nptd.COFX(XI,AI,BI,CI)
        if (CI != 0.0):
            return
    SINGLR= True
    return
def DEFE4(COFX,IDMN,USOL,GRHS):
    for I in range(IS,MS):
        XI=AIT+(I-1)*DLX
        COFX(XI,AI,BI,CI)
        for J in range(JS,NS):
            DX4(USOL,IDMN,I,J,UXXX,UXXXX)
            DY4(USOL,IDMN,I,J,UYYY,UYYYY)
            TX=AI*UXXXX/12.0+BI*UXXX/6.0
            TY=UYYYY/12.0
            if (KSWX !=1 or (I<1 and I>K)):
                TX=AI/3.0*(UXXXX/4.0+UXXX/DLX)
            if (KSWY!=1 or (J<1 and J>L)):
                TY = (UYYYY/4.0+UYYY/DLY)/3.0
            GRHS(I,J)=GRHS(I,J)+DLY**2*(DLX**2*TX+DLY**2*TY)
    for I in range(IS,MS):
        for J in range(JS,NS):
            USOL(I,J)=GRHS(I,J)
    return


def MINSO4(USOL,IDMN,ZN,ZM,PERTB):
    ISTR=1
    IFNL=k
    JSTR=1
    JFNL=L
    UTE=float(0.000)
    ETE=float(0.000)
    for I in range(IS,MS):
        II=I-IS+1
        for J in range(JS,NS):
            JJ=J-JS+1
            ETE=ETE+ZM(II)*ZN(JJ)
            UTE=UTE+USOL(I,J)*ZM(II)*ZN(JJ)
    PERTB=UTE/ETE
    for I in range(ISTR, IFNL):
        for J in range(JSTR, JFNL):
            USOL(I,J)=USOL(I,J)-PERTB
    return
def GENBUN(NPEROD,N,MPEROD,M,A,B,C,IDIMY,Y,IERROR,W):
    IERROR=0
    if(M<=2):
        IERROR=1
    if(N<=2):
        IERROR=2
    if(IDIMY<M):
        IERROR=3
    if(NPEROD<0 or NPEROD>4):
        IERROR=4
    if(MPEROD<0 or MPEROD>1):
        IERROR=5
    if(MPEROD==1):
        if(A[1]!=0 or C[M]!=0):
            IERROR=7
        else:
            if(IERROR!=0):
                return
    for I in range(2,M):
        if(A[I]!=C[1]):
            IERROR=6
        if(C[I]!=C[1]):
            IERROR=6
        if(B[I]!=B[1]):
            IERROR=6
    MP1=M+1
    IWBA=MP1
    IWBB=IWBA+M
    IWBC=IWBB+M
    IWB2=IWBC+M
    IWB3=IWB2+M
    IWW1=IWB3+M
    IWW2=IWW1+M
    IWW3=IWW2+M
    IWD=IWW3+M
    IWTCOS=IWD+M
    IWP=IWTCOS+4*N
    for I in range(1,M):
        K=IWBA+I-1
        W[K]=-A[I]
        K=IWBC+I-1
        W[K]=-C[I]
        K=IWBB+I-1
        W[K]=float(2.00)-B[I]
        for J in range(1,N):
            Y[I][J]=-Y[I][J]
    MP=MPEROD+1
    NP=NPEROD+1
    POISP2(M,N,W(IWBA),W(IWBB),W(IWBC),Y,IDIMY,W,W(IWB2),W(IWB3),W(IWW1),W(IWW2),W(IWW3),W(IWD),W(IWTCOS), W(IWP))
    POISD2 (M,N,1,W(IWBA),W(IWBB),W(IWBC),Y,IDIMY,W,W(IWW1),W(IWD),W(IWTCOS),W(IWP))
    POISN2 (M,N,1,2,W(IWBA),W(IWBB),W(IWBC),Y,IDIMY,W,W(IWB2),W(IWB3),W(IWW1),W(IWW2),W(IWW3),W(IWD),W(IWTCOS),W(IWP))
    POISN2 (M,N,1,1,W(IWBA),W(IWBB),W(IWBC),Y,IDIMY,W,W(IWB2),W(IWB3),W(IWW1),W(IWW2),W(IWW3),W(IWD),W(IWTCOS),W(IWP))
    IPSTOR=W[IWW1]
    IREV=2
    if(NPEROD==4):
        for J in range(1,NBY2):
            MSKIP=N+1-J
            for I in range(1,M):
                A1=Y[I][J]
                Y[I][J]=Y[I][MSKIP]
                Y[I][MSKIP]=A1
                continue
            continue
    MH = (M+1)/2
    MHM1 = MH-1
    MODD = 1
    if(MH==M):
        MODD=2
    for J in range(1,N):
        for I in range(1,MHM1):
            MHPI=MH+I
            MHMI=MH-I
            W[I]=Y[MHMI][J]-Y[MHPI][J]
            W[MHPI]=Y[MHMI][J]+Y[MHPI][J]
            continue
        W[MH]=float(2.00)*Y[MH][J]
        W[M]=float(2.00)*Y[M][J]
        continue
    for I in range(1,M):
        Y[I][J]=W[I]
        continue
    K=IWBC+MHM1-1
    I = IWBA+MHM1
    W[K] = float(0.00)
    W[I] = float(0.00)
    W[K+1] = float(2.00)*W[K+1]
    K=IWBB+MHM1-1
    W[K]=W[K]-W[I-1]
    W[IWBC-1]=W[IWBC-1]+W[IWBB-1]
    W[IWBB-1]=W[K+1]
    IREV = 1
    NBY2 = N/2
    for J in range(1,NBY2):
        MSKIP=N+1-J
        for I in range(1,M):
            A1=Y[I][J]
            Y[I][J]=Y[I][MSKIP]
            Y[I][MSKIP]=A1
            continue
        continue
    for J in range(1,N):
        for I in range(1,MHM1):
            MHMI=MH-I
            MHPI = MH+I
            W[MHMI] =float(0.50)*(Y[MHPI][J]+Y[I][J])
            W[MHPI] =float(0.50)*(Y[MHPI][J]-Y[I][J])
            continue
        W[MH]=float(0.50)*Y[MH][J]
        W[M]=float(0.50)*Y[M][J]
        continue
    for I in range(1,M):
        Y[I][J]=W[I]
        continue
    W[1]=IPSTOR+IWP-1
    return
def MERGE (TCOS,I1,M1,I2,M2,I3):
    J1 = 1
    J2 = 1
    J = I3