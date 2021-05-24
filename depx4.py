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
    r=SPELIA4(IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,D,N,NBDCND,BDC,BDD,COFX,W(I1),W(I2),W(I3),W(I4),W(I5),W(I6),W(I7),W(I8),W(I9),W(I10),W(I11),W(I12),GRHS,USOL,IDMN,W(I13),PERTRB,IERROR)
    return r
# Aca van a ir todos tus metodos Andres