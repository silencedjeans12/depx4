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
