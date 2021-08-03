# Metodos de programacion de Yaritza
def DY4(U,IDMN,I,J,UYYY,UYYYY):
    if(J>2 and J<(L-1)):
        #
        # COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS ON THE INTERIOR
        #
        UYYY = (-U[I][J-2]+np.double(2.00)*U[I][J-1]-np.double(2.00)*U[I][J+1]+U[I][J+2])/TDLY3
        UYYYY = (U[I][J-2]-np.double(4.00)*U[I][J-1]+np.double(6.00)*U[I][J]-np.double(4.00)*U[I][J+1]+U[I][J+2])/DLY4
        return
    if(J==1):
        # COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=C
        if(KSWY==1):
            #
            #PERIODIC AT X=A
            #
            UYYY = (-U[I][L-2]+np.double(2.00)*U[I][L-1]-np.double(2.00)*U[I][2]+U[I][3])/TDLY3
            UYYYY = (U[I][L-2]-np.double(4.00)*U[I][L-1]+np.double(6.00)*U[I][1]-np.double(4.00)*U[I][2]+U[I][3])/DLY4 
            return
        else:
            #
            UYYY = (-np.double(5.00)*U[I][1]+np.double(18.00)*U[I][2]-np.double(24.00)*U[I][3]+np.double(14.00)*U[I][4]-np.double(3.00)*U[I][5])/TDLY3
            UYYYY = (np.double(3.00)*U[I][1]-np.double(14.00)*U[I][2]+np.double(26.00)*U[I][3]-np.double(24.00)*U[I][4]+np.double(11.00)*U[I][5]-np.double(2.00)*U[I][6])/DLY4
            return
    if(J==2):
        if(KSWY==1):
            UYYY = (-U[I][L-1]+np.double(2.00)*U[I][1]-np.double(2.00)*U[I][3]+U[I][4])/TDLY3
            UYYYY=(U[I][L-1]-np.double(4.00)*U[I][1]+np.double(6.00)*U[I][2]-np.double(4.00)*U[I][3]+U[I][4])/DLY4
            return
        else:
            UYYY = (-np.double(3.00)*U[I][1]+np.double(10.00)*U[I][2]-np.double(12.00)*U[I][3]+np.double(6.00)*U[I][4]-U[I][5])/TDLY3
            UYYYY = (np.double(2.00)*U[I][1]-np.double(9.00)*U[I][2]+np.double(16.00)*U[I][3]-np.double(14.00)*U[I][4]+np.double(6.00)*U[I][5]-U[I][6])/DLY4
            return
    if(J==L-1):
        if(KSWY==1):
           UYYY = (-U[I][L-3]+np.double(2.00)*U[I][L-2]-np.double(2.00)*U[I][1]+U[I][2])/TDLY3 
           UYYYY = (U[I][L-3]-np.double(4.00)*U[I][L-2]+np.double(6.00)*U[I][L-1]-np.double(4.00)*U[I][1]+U[I][2])/DLY4
           return
        else:
            UYYY = (U[I][L-4]-np.double(6.00)*U[I][L-3]+np.double(12.00)*U[I][L-2]-np.double(10.00)*U[I][L-1]+np.double(3.00)*U[I][L])/TDLY3
            UYYYY = (-U[I][L-5]+np.double(6.00)*U[I][L-4]-np.double(14.00)*U[I][L-3]+np.double(16.00)*U[I][L-2]-np.double(9.00)*U[I][L-1]+np.double(2.00)*U[I][L])/DLY4
            return
    if(J==L):
        UYYY = -(np.double(3.00)*U[I][L-4]-np.double(14.00)*U[I][L-3]+np.double(24.00)*U[I][L-2]-np.double(18.00)*U[I][L-1]+np.double(5.00)*U[I][L])/TDLY3
        UYYYY = (-2.0*U[I][L-5]+11.0*U[I][L-4]-24.0*U[I][L-3]+26.0*U[I][L-2]-14.0*U[I][L-1]+3.0*U[I][L])/DLY4
        return
def DX4(U,IDMN,I,J,UXXX,UXXXX):
    if(I>2 and I<(K-1)):
        UXXX = (-U[I-2][J]+np.double(2.00)*U[I-1][J]-np.double(2.00)*U[I+1][J]+U[I+2][J])/TDLX3
        UXXXX = (U[I-2][J]-np.double(4.00)*U[I-1][J]+np.double(6.00)*U(I,J)-np.double(4.00)*U[I+1][J]+U[I+2][J])/DLX4
        return
    if(I==1):
        if(KSWX==1):
            UXXX = (-U[K-2][J]+np.double(2.00)*U[K-1][J]-np.double(2.00)*U[2][J]+U[3][J])/(TDLX3)
            UXXXX = (U[K-2][J]-np.double(4.00)*U[K-1][J]+np.double(6.00)*U[1][J]-np.double(4.00)*U[2][J]+U[3][J])/DLX4
            return
        else:
            UXXX = (-np.double(5.00)*U[1][J]+np.double(18.00)*U[2][J]-np.double(24.00)*U[3][J]+np.double(14.00)*U[4][J]-np.double(3.00)*U[5][J])/(TDLX3)
            UXXXX = (np.double(3.00)*U[1][J]-np.double(14.00)*U[2][J]+np.double(26.00)*U[3][J]-np.double(24.00)*U[4][J]+np.double(11.00)*U[5][J]-np.double(2.00)*U(6,J))/DLX4
            return
    if(I==2):
        if(KSWX==1):
            UXXX = (-U[K-1][J]+np.double(2.00)*U[1][J]-np.double(2.00)*U[3][J]+U[4][J])/(TDLX3)
            UXXXX = (U[K-1][J]-np.double(4.00)*U[1][J]+np.double(6.00)*U[2][J]-np.double(4.00)*U[3][J]+U[4][J])/DLX4
            return
        else:
            UXXX = (-np.double(3.00)*U[1][J]+np.double(10.00)*U[2][J]-np.double(12.00)*U[3][J]+np.double(6.00)*U[4][J]-U[5][J])/TDLX3
            UXXXX = (np.double(2.00)*U[1][J]-np.double(9.00)*U[2][J]+np.double(16.00)*U[3][J]-np.double(14.00)*U[4][J]+np.double(6.00)*U[5][J]-U[6][J])/DLX4
            return
    if(I==K-1):
        if(KSWX==1):
            UXXX = (-U[K-3][J]+np.double(2.00)*U[K-2][J]-np.double(2.00)*U[1][J]+U[2][J])/TDLX3
            UXXXX = (U[K-3][J]-np.double(4.00)*U[K-2][J]+np.double(6.00)*U[k-1][J]-np.double(4.00)*U[1][J]+U[2][J])/DLX4
            return
        else:
            UXXX = (U[K-4][J]-np.double(6.00)*U[K-3][J]+np.double(12.00)*U[K-2][J]-np.double(10.00)*U[k-1][J]+np.double(3.00)*U[K][J])/TDLX3
            UXXXX = (-U[K-5][J]+np.double(6.00)*U[K-4][J]-np.double(14.00)*U[K-3][J]+np.double(16.00)*U[K-2][J]-np.double(9.00)*U[k-1][J]+np.double(2.00)*U[K][J])/DLX4
            return
    if(I==K):
        UXXX = -(np.double(3.00)*U[K-4][J]-np.double(14.00)*U[K-3][J]+np.double(24.00)*U[K-2][J]-np.double(18.00)*U[k-1][J]+np.double(5.00)*U[K][J])/TDLX3
        UXXXX = (-np.double(2.00)*U[K-5][J]+np.double(11.00)*U[K-4][J]-np.double(24.00)*U[K-3][J]+np.double(26.00)*U[K-2][J]-np.double(14.00)*U[k-1][J]+np.double(3.00)*U[K][J])/DLX4
        return
def ORTHO4(USOL,IDMN,ZN,ZM,PERTRB):
    ISTR=IS
    IFNL=MS
    JSTR=JS
    JFNL=NS

    UTE=np.double(0.00)
    ETE=np.double(0.00)
    for I in range(IS,MS):
        II=I-IS+1
        for J in range(JS,NS):
            JJ=J-JS+1
            ETE=ETE+ZM[II]*ZN[JJ]
            UTE=UTE+USOL[I][J]*ZM[II]*ZN[JJ]
    
    PERTRB=UTE/ETE
    for I in range(ISTR,IFNL):
        for J in range(JSTR,JFNL):
            USOL[I][J]=USOL[I][J]-PERTRB
    return
