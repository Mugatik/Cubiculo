!  Created by Aitor Urresti Gonzalez on 15/02/11.


module props
implicit none
save

real, allocatable::TempH(:,:)
real::ks,kl,Ts,Tl
real::cs,cm,cl

contains

real function TH(input)
implicit none

!Funci�n que calcula la temperatura que corresponde a una entalp�a


real::input
integer::pos1,N

N=size(TempH,1)

if (input<=TempH(1,2)) then
    TH=TempH(1,1)-(TempH(1,2)-Input)*(TempH(2,1)-TempH(1,1))/(TempH(2,2)-TempH(1,2))
elseif (input>=TempH(N,2)) then
    TH=TempH(N,1)+(input-TempH(N,2))*(TempH(N,1)-TempH(N-1,1))/(TempH(N,2)-TempH(N-1,2))
else
    pos1=minloc(abs(TempH(:,2)-input),1)
    if (TempH(pos1,2)<=input) then
        TH=TempH(pos1,1)+(TempH(pos1+1,1)-TempH(pos1,1))/(TempH(pos1+1,2)-TempH(pos1,2))*(Input-TempH(pos1,2))
        else
        TH=TempH(pos1-1,1)+(TempH(pos1,1)-TempH(pos1-1,1))/(TempH(pos1,2)-TempH(pos1-1,2))*(Input-TempH(pos1-1,2))
    endif
endif

end function

real function HT(input)
implicit none

!Funci�n que calcula la entalp�a que corresponde a una temperatura

real::input
integer::pos1,N

N=size(TempH,1)

if (input<=TempH(1,1)) then
    HT=TempH(1,2)-(TempH(1,1)-Input)*(TempH(2,2)-TempH(1,2))/(TempH(2,1)-TempH(1,1))
elseif (input>=TempH(N,1)) then
    HT=TempH(N,2)+(input-TempH(N,1))*(TempH(N,2)-TempH(N-1,2))/(TempH(N,1)-TempH(N-1,1))
else
    pos1=1
    do while (input>=TempH(pos1+1,1))
        pos1=pos1+1
    end do
    HT=TempH(pos1,2)+(Input-TempH(pos1,1))*(TempH(pos1+1,2)-TempH(pos1,2))/(TempH(pos1+1,1)-TempH(pos1,1))
endif

end function

real function k(T1,T2)
implicit none

!Funci�n que calcula la conductividad intermedia entre los valores de dos temperaturas

real::T1,T2
real::k1,k2

if ((T1<=Ts).and.(T2<=Ts)) then !Estamos en estado s�lido
    k=ks
elseif((T1>=Tl).and.(T2>=Tl)) then  !Estamos en estado l�quido
    k=kl
else
    if(T1<=Ts)then
        k1=ks
    elseif(T1>=Tl)then
        k1=kl
    else
        k1=ks+(T1-Ts)*(kl-ks)/(Tl-Ts)
    endif
    if(T2<=Ts)then
        k2=ks
    elseif(T2>=Tl)then
        k2=kl
    else
        k2=ks+(T2-Ts)*(kl-ks)/(Tl-Ts)
    endif
    k=2*k1*k2/(k1+k2)
endif

end function

real function HTLin(input)
implicit none

!Funci�n que calcula la entalp�a de una temperatura usando la relaci�n linealizada

real::input

if(input<=Ts)then
    HTLin=input*cs
elseif(input<=Tl)then
    HTLin=Ts*cs+(input-Ts)*cm
else
    HTLin=Ts*cs+(Tl-Ts)*cm+(input-Tl)*cl
endif

end function

real function THLin(input)
implicit none

!Funci�n que calcula la temperatura que corresponde a una temperatura usando la relaci�n linealizada

real::input
real::hs,hl

hs=cs*ts
hl=hs+cm*(tl-ts)

if(input<=hs)then
    THLin=input/cs
elseif(input<=Hl)then
    THLin=Ts+(input-hs)/cm
else
    THLin=Tl+(input-hl)/cl
endif

end function

real function PhiThetaLin(input)
implicit none

!Funci�n que obtiene la entalp�a adimensional correspondiente a una temperatura adimensional usanlo la relaci�n lienal

real::input
real::Thetal

Thetal=cs/cm

if(input<=0)then
    PhiThetaLin=input
elseif(input<=Thetal)then
    PhiThetaLin=input/Thetal
else
    PhiThetaLin=1+cl*(input-Thetal)/cs
endif

end function

real function ThetaPhiLin(input)
implicit none

!Funci�n que obtiene la tempratura adimensional correspondiente a una entalp�a adimensional usanlo la relaci�n lienal

real::input
real::Thetal

Thetal=cs/cm

if(input<=0)then
    ThetaPhiLin=input
elseif(input<=1)then
    ThetaPhiLin=input*Thetal
else
    ThetaPhiLin=(Input-1)*cs/cl+Thetal
endif

end function

end module
