# Metodos de programacion de Yaritza
def DY4(U,IDMN,I,J,UYYY,UYYYY):
    if(J>2 and J<(L-1)):
        #
        # COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS ON THE INTERIOR
        #
        UYYY = (-U[I][J-2]+float(2.00)*U[I][J-1]-float(2.00)*U[I][J+1]+U[I][J+2])/TDLY3
        UYYYY = (U[I][J-2]-float(4.00)*U[I][J-1]+float(6.00)*U[I][J]-float(4.00)*U[I][J+1]+U[I][J+2])/DLY4
        return
    if(J==1):
        # COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=C
        if(KSWY==1):
            #
            #PERIODIC AT X=A
            #
            UYYY = (-U[I][L-2]+float(2.00)*U[I][L-1]-float(2.00)*U[I][2]+U[I][3])/TDLY3
            UYYYY = (U[I][L-2]-float(4.00)*U[I][L-1]+float(6.00)*U[I][1]-float(4.00)*U[I][2]+U[I][3])/DLY4 
            return
        else:
            #
            UYYY = (-float(5.00)*U[I][1]+float(18.00)*U[I][2]-float(24.00)*U[I][3]+float(14.00)*U[I][4]-float(3.00)*U[I][5])/TDLY3
            UYYYY = (float(3.00)*U[I][1]-float(14.00)*U[I][2]+float(26.00)*U[I][3]-float(24.00)*U[I][4]+float(11.00)*U[I][5]-float(2.00)*U[I][6])/DLY4
            return
    if(J==2):
        if(KSWY==1):
            UYYY = (-U[I][L-1]+float(2.00)*U[I][1]-float(2.00)*U[I][3]+U[I][4])/TDLY3
            UYYYY=(U[I][L-1]-float(4.00)*U[I][1]+float(6.00)*U[I][2]-float(4.00)*U[I][3]+U[I][4])/DLY4
            return
        else:
            UYYY = (-float(3.00)*U[I][1]+float(10.00)*U[I][2]-float(12.00)*U[I][3]+float(6.00)*U[I][4]-U[I][5])/TDLY3
            UYYYY = (float(2.00)*U[I][1]-float(9.00)*U[I][2]+float(16.00)*U[I][3]-float(14.00)*U[I][4]+float(6.00)*U[I][5]-U[I][6])/DLY4
            return
    if(J==L-1):
        if(KSWY==1):
           UYYY = (-U[I][L-3]+float(2.00)*U[I][L-2]-float(2.00)*U[I][1]+U[I][2])/TDLY3 
           UYYYY = (U[I][L-3]-float(4.00)*U[I][L-2]+float(6.00)*U[I][L-1]-float(4.00)*U[I][1]+U[I][2]))/DLY4
           return
        else:
            UYYY = (U[I][L-4]-float(6.00)*U[I][L-3]+float(12.00)*U[I][L-2]-float(10.00)*U[I][L-1]+float(3.00)*U[I][L])/TDLY3
            UYYYY = (-U(I,L-5)+float(6.00)*U[I][L-4]-float(14.00)*U[I][L-3]+float(16.00)*U[I][L-2]-float(9.00)*U[I][L-1]+float(2.00)*U[I][L])/DLY4
            return
    if(J==L):
        UYYY = -(float(3.00)*U[I][L-4]-float(14.00)*U[I][L-3]+float(24.00)*U[I][L-2]-float(18.00)*U[I][L-1]+float(5.00)*U[I][L])/TDLY3
        UYYYY = (-2.0*U(I,L-5)+11.0*U[I][L-4]-24.0*U[I][L-3]+26.0*U[I][L-2]-14.0*U[I][L-1]+3.0*U[I][L])/DLY4
        return