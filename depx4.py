import numpy as np
def DEPX4(IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,D,N,NBDCND,BDC,BDD,COFX,GRHS,USOL,IDMN,W,PERTRB,IERROR):
    GRHS=np.zeros((IDMN,1))
    USOL=np.zeros((IDMN,1))
    BDA=np.zeros((1))
    BDB=np.zeros((1))
    BDC=np.zeros((1))
    BDD=np.zeros((1))
    W=np.zeros((1))
    CHKPR4(IORDER,A,B,M,MBDCND,C,D,N,NBDCND,COFX,IDMN,IERROR)
    if IERROR != 0:
        return
    L=N+1
    if NBDCND == 0:
        L=N
    K=M+1
    L=N+1
#   ESTIMATE LOG BASE 2 OF N
    
    LOG2N=int(np.log(np.double(N+1))/np.log(np.double(2.0))+np.double(0.50))
    LENGTH=4*(N+1)+(10+LOG2N)*(M+1)
    IERROR = 11
    LINPUT = int(W[0]+0.5)
    LOUTPT = LENGTH+6*(K+L)+1
    W[0]=(np.double(LOUTPT))
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
    SPELIA4(IORDER,A,B,M,MBDCND,BDA,ALPHA,BDB,BETA,C,D,N,NBDCND,BDC,BDD,COFX,W[I1],W[I2],W[I3],W[I4],W[I5],W[I6],W[I7],W[I8],W[I9],W[I10],W[I11],W[I12],GRHS,USOL,IDMN,W[I13],PERTRB,IERROR) 
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
    DLY=(DIT-CIT)/np.double(N)
    for I in range(2,M):
        for J in range(2,N):
            USOL[I][J]=DLY**2**GRHS[I][J]
            
        
    if (KSWX==2 or KSWX==3):
        f=0
def CHKPR4(IORDER,A,B,M,MBDCND,C,D,N,NBDCND,COFX,IDMN,IERROR):
#
#     THIS PROGRAM CHECKS THE INPUT PARAMETERS FOR ERRORS
#
#
#
#     CHECK DEFINITION OF SOLUTION REGION
#
    IERROR=1
    if (A>=B or C>=D):
        return
#
#     CHECK BOUNDARY SWITCHES
#
    IERROR=2
    if (MBDCND<0 or MBDCND>4):
        return
    IERROR=3
    if (NBDCND<0 or NBDCND>4):
        return
#
#     CHECK FIRST DIMENSION IN CALLING ROUTINE
#
    IERROR=5
    if (IDMN<7):
        return
#
#   CHECK M
#
    IERROR=6
    if (M>(IDMN-1)or M<6):
        return
#
#   CHECK N
#
    IERROR=7
    if (N<5):
        return
#
#   CHECK IORDER
#
    IERROR=8
    if (IORDER!=2 and IORDER!=4):
        return
#
#     CHECK INTL
#
#
#     CHECK THAT EQUATION IS ELLIPTIC
#
    DLX=(B-A)/np.double(M)
    for I in range(2,M):
        XI=A+np.double(I-1)*DLX
        COFX(XI,AI,BI,CI)
    if (AI<np.double(0.00)):
        IERROR=10
#
#     NO ERROR FOUND
#
    IERROR=0
def CHKSN4(MBDCND,NBDCND,ALPHA,BETA,COFX,SINGLR):
#
#     THIS SUBROUTINE CHECKS IF THE PDE   SEPX4
#     MUST SOLVE IS A SINGULAR OPERATOR
#
    SINGLR=False
#
#     CHECK IF THE BOUNDARY CONDITIONS ARE
#     ENTIRELY PERIODIC AND/OR MIXED
#
    if ((MBDCND!=0 and MBDCND!=3)or(NBDCND!=0 and NBDCND!=3)):
        return
#
#     CHECK THAT MIXED CONDITIONS ARE PURE NEUMAN
#
    if (MBDCND==3):
        if (ALPHA!=np.double(0.00)or BETA!=np.double(0.00)):
            return
#
#     CHECK THAT NON-DERIVATIVE COEFFICIENT FUNCTIONS
#     ARE ZERO
#
    for I in range(IS,MS):
        XI=AIT+(I-1)*DLX
        COFX(XI,AI,BI,CI)
        if (CI != 0.0):
            return
    SINGLR= True
def DEFE4(COFX,IDMN,USOL,GRHS):
#
#     THIS SUBROUTINE FIRST APPROXIMATES THE TRUN1ATION ERROR GIVEN BY
#     TRUN1(X,Y)=DLX**2*TX+DLY**2*TY WHERE
#     TX=AFUN(X)*UXXXX/12.0+BFUN(X)*UXXX/6.0 ON THE INTERIOR AND
#     AT THE BOUNDARIES IF PERIODIC(HERE UXXX,UXXXX ARE THE THIRD
#     AND FOURTH PARTIAL DERIVATIVES OF U WITH RESPECT TO X).
#     TX IS OF THE FORM AFUN(X)/3.0*(UXXXX/4.0+UXXX/DLX)
#     AT X=A OR X=B IF THE BOUNDARY CONDITION THERE IS MIXED.
#     TX=0.0 ALONG SPECIFIED BOUNDARIES.  TY HAS SYMMETRIC FORM
#     IN Y WITH X,AFUN(X),BFUN(X) REPLACED BY Y,DFUN(Y),EFUN(Y).
#     THE SECOND ORDER SOLUTION IN USOL IS USED TO APPROXIMATE
#     (VIA SECOND ORDER FINITE DIFFERENCING) THE TRUN1ATION ERROR
#     AND THE RESULT IS ADDED TO THE RIGHT HAND SIDE IN GRHS
#     AND THEN TRANSFERRED TO USOL TO BE USED AS A NEW RIGHT
#     HAND SIDE WHEN CALLING BLKTRI FOR A FOURTH ORDER SOLUTION.
#
    GRHS=np.zeros((IDMN,1))
    USOL=np.zeros((IDMN,1))
#
#
#     COMPUTE TRUN1ATION ERROR APPROXIMATION OVER THE ENTIRE MESH
#
    for I in range(IS,MS):
        XI=AIT+(I-1)*DLX
        COFX(XI,AI,BI,CI)
        for J in range(JS,NS):
#
#     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT (XI,YJ)
#
            DX4(USOL,IDMN,I,J,UXXX,UXXXX)
            DY4(USOL,IDMN,I,J,UYYY,UYYYY)
            TX=AI*UXXXX/12.0+BI*UXXX/6.0
            TY=UYYYY/12.0
#
#     RESET FORM OF TRUN1ATION IF AT BOUNDARY WHICH IS NON-PERIODIC
#
            if (KSWX !=1 or (I<1 and I>K)):
                TX=AI/3.0*(UXXXX/4.0+UXXX/DLX)
            if (KSWY!=1 or (J<1 and J>L)):
                TY = (UYYYY/4.0+UYYY/DLY)/3.0
            GRHS[I][J]=GRHS[I][J]+DLY**2*(DLX**2*TX+DLY**2*TY)
#
#     RESET THE RIGHT HAND SIDE IN USOL
#       
    for I in range(IS,MS):
        for J in range(JS,NS):
            USOL[I][J]=GRHS[I][J]

def MINSO4(USOL,IDMN,ZN,ZM,PERTB):
    USOL=np.zeros((IDMN,1))
    ZN=np.zeros((1))
    ZM=np.zeros((1))
    ISTR=1
    IFNL=k
    JSTR=1
    JFNL=L
    UTE=np.double(0.000)
    ETE=np.double(0.000)
    for I in range(IS,MS):
        II=I-IS+1
        for J in range(JS,NS):
            JJ=J-JS+1
            ETE=ETE+ZM[II]*ZN[JJ]
            UTE=UTE+USOL[I][J]*ZM[II]*ZN[JJ]
            
        
    PERTB=UTE/ETE
    for I in range(ISTR, IFNL):
        for J in range(JSTR, JFNL):
            USOL[I][J]=USOL[I][J]-PERTB
            
        
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
        W[K]=np.double(2.00)-B[I]
        for J in range(1,N):
            Y[I][J]=-Y[I][J]
            
        
    MP=MPEROD+1
    NP=NPEROD+1

    #Inico respaldo
    if MP ==1:
        MH=(M+1)/2
        MHM1 = MH-1
        MODD = 1
        if (MH*2==M):
            MODD=2
        for J in range(1,N):
            for I in range(1,MHM1):
                MHPI=MH+I
                MHMI=MH-I
                W[I]=Y[MHMI][J]-Y[MHPI][J]
                W[MHPI]=Y[MHMI][J]+Y[MHPI][J]
                
            W[MH]=np.double(2.00)*Y[MH][J]
            if MODD==1:
                for I in range(1,M):
                    Y[I][J]=W[I]
                    
                
            if MODD==2:
                W[M]=np.double(2.00)*Y[M][J]
                for I in range(1,M):
                    Y[I][J]=W[I]
                    
                
            K= IWBC+MHM1-1
            I = IWBA+MHM1
            W[K] = np.double(0.000)
            W[I] = np.double(0.000)
            W[K+1] = np.double(2.00)*W[K+1]
            if MODD==1:
                K=IWBB+MHM1-1
                W[K]=W[K]-W[I-1]
                W[IWBC-1]=W[IWBC-1]+W[IWBB-1]
            if MODD==2:
                W[IWBB-1]=W[K+1]
            
    if MP==2:
        if NP==1:
            POISP2(M,N,W[IWBA],W[IWBB],W[IWBC],Y,IDIMY,W,W[IWB2],W[IWB3],W[IWW1],W[IWW2],W[IWW3],W[IWD],W[IWTCOS],W[IWP])
        if NP==2:
            POISD2(M,N,1,W[IWBA],W[IWBB],W[IWBC],Y,IDIMY,W,W[IWW1],W[IWD],W[IWTCOS],W[IWP])
        if NP==3:
            POISN2(M,N,1,2,W[IWBA],W[IWBB],W[IWBC],Y,IDIMY,W,W[IWB2],W[IWB3],W[IWW1],W[IWW2],W[IWW3],W[IWD],W[IWTCOS],W[IWP])
        if NP==4:
            POISN2(M,N,1,1,W[IWBA],W[IWBB],W[IWBC],Y,IDIMY,W,W[IWB2],W[IWB3],W[IWW1],W[IWW2],W[IWW3],W[IWD],W[IWTCOS],W[IWP])
        IPSTOR=W[IWW1]
        IREV=2
        if NPEROD==4:
            for J in range(1,NBY2):
                MSKIP=N+1-J
                for I in range(1,M):
                    A1=Y[I][J]
                    Y[I][J]=Y[I][MSKIP]
                    Y[I][MSKIP]=A1
                    
                
        if NP==5:
            IREV = 1
            NBY2 = N/2
        for J in range(1,NBY2):
            MSKIP=N+1-J
            for I in range(1,M):
                A1=Y[I][J]
                Y[I][J]=Y[I][MSKIP]
                Y[I][MSKIP]=A1
                
            
#Fnal
    
    
    
    
    IPSTOR=W[IWW1]
    IREV=2
    if(NPEROD==4):
        for J in range(1,NBY2):
            MSKIP=N+1-J
            for I in range(1,M):
                A1=Y[I][J]
                Y[I][J]=Y[I][MSKIP]
                Y[I][MSKIP]=A1
                
            
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
            
        W[MH]=np.double(2.00)*Y[MH][J]
        W[M]=np.double(2.00)*Y[M][J]
        
    for I in range(1,M):
        Y[I][J]=W[I]
        
    K=IWBC+MHM1-1
    I = IWBA+MHM1
    W[K] = np.double(0.00)
    W[I] = np.double(0.00)
    W[K+1] = np.double(2.00)*W[K+1]
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
            
        
    for J in range(1,N):
        for I in range(1,MHM1):
            MHMI=MH-I
            MHPI = MH+I
            W[MHMI] =np.double(0.50)*(Y[MHPI][J]+Y[I][J])
            W[MHPI] =np.double(0.50)*(Y[MHPI][J]-Y[I][J])
            
        W[MH]=np.double(0.50)*Y[MH][J]
        W[M]=np.double(0.50)*Y[M][J]
        
    for I in range(1,M):
        Y[I][J]=W[I]
        
    W[1]=IPSTOR+IWP-1
    return
def POISD2 (MR,NR,ISTAG,BA,BB,BC,Q,IDIMQ,B,W,D,TCOS,P):
    Q=np.zeros((IDIMQ,1))
    BA=np.zeros((1))
    BB=np.zeros((1))
    BC=np.zeros((1))
    TCOS=np.zeros((1))
    B=np.zeros((1))
    D=np.zeros((1))
    W=np.zeros((1))
    P=np.zeros((1))
    M=MR
    N=NR
    JSH=0
    FI=1.0/np.double(ISTAG)
    IP=-M
    IPSTOR=0
    if ISTAG==1:#--->101
        KR=0
        IRREG=1
        if N>1:#---->106
            return
    elif ISTAG==2:#--->102
        KR=1
        JSTSAV=1
        IRREG=2
    if ISTAG==1:#---->101
        #101
        KR=0
        IRREG=1
        if N>1:#----->106
            #106
            LR=0
            for I in range(1,M):
                P[I]=np.double(0.00)
                
            NUN=N
            JST=1
            JSP=N
            L=2*JST
            NODD=2-2*((NUN+1)/2)+NUN
            while IRREG==1:

                if NODD==1:#------>110
                    #110
                    JSP=JSP-JST
                    if IRREG!=1:
                        JSP=JSP-L
                if NODD==2:#---->109
                    #109
                    JSP=JSP-L
                COSGEN(JST,1,np.double(0.50),np.double(0.00),TCOS)
                if L<JSP:#---->118
                    for J in range(L,JSP,L):
                        JM1 = J-JSH
                        JP1 = J+JSH
                        JM2 = J-JST
                        JP2 = J+JST
                        JM3 = JM2-JSH
                        JP3 = JP2+JSH
                        if JST!=1:
                            for I in range(1,M):
                                T = Q[I][J]-Q[I][JM1]-Q[I][JP1]+Q[I][JM2]+Q[I][JP2]
                                B[I] = T+Q[I][J]-Q[I][JM3]-Q[I][JP3]
                                Q[I][J] = T
                        else:
                            for I in range(1,M):
                                B[I]=np.double(2.00)*Q[I][J]
                                Q[I][J]=Q[I][JM2]+Q[I][JP2]
                        TRIX(JST,0,M,BA,BB,BC,B,TCOS,D,W)
                        for I in range(1,M):
                            Q[I][J]=Q[I][J]+B[I]
                            
                        
                if NODD==1:
                    if IRREG==1:
                        NUN=NUN/2
                        NODDPR=NODD
                        JSH=JST
                        JST=2*JST
                        if NUN>=2:#---->108
                            L = 2*JST
                            NODD = 2-2*((NUN+1)/2)+NUN
                        else:
                            J=JSP
                            for I in range(1,M):
                                B[I]=Q[I][J]
                        if IRREG==1:#--->154
                            #154
                            COSGEN(JST,1,np.double(0.50),np.double(0.00),TCOS)
                            IDEG=JST
                        if IRREG==2:#-->155
                            #155
                            KR=LR+JST
                            COSGEN(KR,JSTSAV,np.double(0.00),FI,TCOS)
                            COSGEN(LR,JSTSAV,np.double(0.00),FI,TCOS(KR+1))
                            IDEG=KR
                        TRIX(IDEG,LR,M,BA,BB,BC,B,TCOS,D,W)
                        JM1=J-JSH
                        JP1 = J+JSH            
                        if IRREG==1:#--->157
                            #157
                            for I in range(1,M):
                                Q[I][J]=np.double(.50)*(Q[I][J]-Q[I][JM1]-Q[I][JP1])+B[I]
                        if IRREG==2:#--->159
                            #159
                            if NODDPR==1:#--->160
                                #160
                                for I in range(1,M):
                                    IP1=IP+I
                                    Q[I][J]=P[IP1]+B[I]
                                IP=IP-M        
                            if NODDPR==2:#--->162
                                #162
                                for I in range(1,M):
                                    Q[I][J]=Q[I][J]-Q[I][JM1]+B[I]
                            JST = JST/2
                            JSH = JST/2
                            NUN = 2*NUN
                            if NUN<N:
                                for J in range(JST,N,L):
                                    JM1 = J-JSH
                                    JP1 = J+JSH
                                    JM2 = J-JST
                                    JP2 = J+JST 
                                    if J<JST:#ELSE DE 166
                                        for I in range(1,M):
                                            B[I]=Q[I][J]+Q[I][JP2]
                                        COSGEN(JST,1,np.double(0.50),np.double(0.00),TCOS)
                                        IDEG = JST
                                        JDEG = 0
                                    else:
                                        if JP2<=N:#168
                                            for I in range(1,M):
                                                B[I]=Q[I][J]+Q[I][JM2]+Q[I][JP2]
                                            COSGEN(JST,1,np.double(0.50),np.double(0.00),TCOS)
                                            IDEG = JST
                                            JDEG = 0
                                        else:
                                            for I in range(1,M):
                                                B[I]=Q[I][J]+Q[I][JM2]
                                            if JST<JSTSAV:
                                                IRREG=1
                                            if IRREG==1:
                                                COSGEN(JST,1,np.double(0.50),np.double(0.00),TCOS)
                                                IDEG = JST
                                                JDEG = 0
                                            if IRREG==2:
                                                if(J+L>N):
                                                    LR=LR-JST
                                                KR=JST+LR
                                                COSGEN(KR,JSTSAV,np.double(0.00),FI,TCOS)
                                                COSGEN(LR,JSTSAV,np.double(0.00),FI,TCOS[KR+1])
                                                IDEG=KR
                                                JDEG=LR
                                    TRIX(IDEG,JDEG,M,BA,BB,BC,B,TCOS,D,W)
                                    if JST>1:
                                        if JP2>N:
                                            if IRREG==1:#--->175
                                                for I in range(1,M):
                                                    Q[I][J]=np.double(0.50)*(Q[I][J]-Q[I][JM1]-Q[I][JP1])+B[I]
                                            if IRREG==2:#---->178
                                                #178
                                                if J+JSH>N:#--->180
                                                    for I in range(1,M):
                                                        Q[I][J]=B[I]+Q[I][J]-Q[I][JM1]
                                                        
                                                    
                                                else:
                                                    for I in range(1,M):
                                                        IP1=IP+I
                                                        Q[I][J]=B[I]+P[IP1]
                                                        
                                                    IP=IP-M
                                                L=L/2
                                                
                                                      


                                       

        else:
            TCOS[0]=np.double(0.00)
            for I in range(1,M):
                B[I]=Q[I][1]
                
            TRIX(1,0,M,BA,BB,BC,B,TCOS,D,W)
            for I in range(1,M):
                Q[I][1]=B[I]
                
            W[0]=IPSTOR
#continuacion.   
    
        
def MERGE (TCOS,I1,M1,I2,M2,I3):
    TCOS=np.zeros((1))
    J1 = 1
    J2 = 1
    J = I3
    if(M1==0):
        K=J+J2+1
        for J in range(J2,M2):
            M=K+J
            L=J+I2
            TCOS[M]=TCOS[L]
    else:
        if(M2==0):
            K=J+J1+1
            for J in range(J1,M1):
                M=K+J
                L=J+I1
                TCOS[M]=TCOS[L]
                
        else:
            J = J+1
            L = J1+I1
            X = TCOS[L]
            L = J2+I2
            Y = TCOS[L]
            TCOS[J]=X
            J1=J1+1
    return
def TRIX (IDEGBR,IDEGCR,M,A,B,C,Y,TCOS,D,W):
    A=np.zeros((1))
    B=np.zeros((1))
    C=np.zeros((1))
    Y=np.zeros((1))
    TCOS=np.zeros((1))
    D=np.zeros((1))
    W=np.zeros((1))
    MM1=M-1
    FB=IDEGBR+1
    FC=IDEGCR+1
    L=FB/FC
    LINT=1
    for K in range(1,IDEGBR):
        X=TCOS[K]
        if K!=L:
            Z=np.double(1.00)/(B[1]-X)
            D[1]=C[1]*Z
            Y[1]=Y[1]*Z
        else:
            I=IDEGBR+LINT
            XX=X-TCOS[I]
            for I in range(1,M):
                W[I]=Y[I]
                Y[I]=XX*Y[I]
            Z=np.double(1.00)/(B[1]-X)
            D[1]=C[1]*Z
            Y[1]=Y[1]*Z   
        for I in range(2,MM1):
            Z=np.double(1.00)/(B[I]-X-A[I]*D[I-1])
            D[I]=C[I]*Z
            Y[I]=(Y[I]-A[I]*Y[I-1])*Z
            Z=B[M]-X-A[M]*D[MM1]
            if Z!=0:
                Y[M]=(Y[M]-A[M]*Y[MM1])/Z
                
            else:
                Y[M]=np.double(0.00)
            for IP in range(1,MM1):
                I=M-IP
                Y[I]=Y[I]-D[I]*Y[I+1]
                
    if K==L:
        for I in range(1,M):
            Y[I]=Y[I]+W[I]
            
        LINT=LINT+1
        L=(np.double(LINT)*FB)/FC
    return
