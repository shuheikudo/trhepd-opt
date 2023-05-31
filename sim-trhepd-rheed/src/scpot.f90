!*******************************************************************
!   RHEED multi slice method :  subroutines for bulk & surf
!    v.3    90/4/12   v.4:2014/4    T.Hanada
!   contains asfcrr asfcef scpot
!*******************************************************************
! # of gaussians for real atomic scattering factor : nabr
! # of gaussians for imaginary atomic scattering factor : nabi
!**********************************************************
!       energy correction of atomic scattering factor
!**********************************************************
        subroutine asfcrr(be,wn,nelm,mab,nabr,nabi,a,b,c,dchg &
                         ,dth,dtk,dtz,elm,tem,smesh)
        implicit none
        integer :: nelm,mab,nabr,nabi
        real(8) :: a(mab,nelm),b(mab,nelm),c(mab,nelm),dchg(nelm)
        real(8) :: dtz(nelm),dth(nelm),dtk(nelm),elm(nelm),be,wn,tem,smesh

        real(8) :: s,rel,rat,vib,dwz
        integer :: i,j
        real(8), parameter :: pi=atan(1d0)*4d0, c2m=511.001d0, ek=.262466d0

        wn=sqrt(1d3*be*ek*(1d0+0.5d0*be/c2m))
        s=4d0*pi/smesh
        rel=(1d0+be/c2m)*s
        rat=sqrt(100d0/be*(1d0+50d0/c2m)/(1d0+be*0.5d0/c2m))*s
        vib=11407d0*tem
        do i=1,nelm
          ! real part
          if (tem > 1d-5) then       ! Debye T
            dwz=vib/(dtz(i)*dtz(i)*elm(i))                 ! B
            dth(i)=sqrt((vib+vib)/elm(i))/(4d0*pi*dth(i))  ! sqrt(<u^2>)
            dtk(i)=sqrt((vib+vib)/elm(i))/(4d0*pi*dtk(i))  ! sqrt(<u^2>)
          else if (tem < -1d-5) then ! sqrt(<u^2>)
            dwz=(8d0*pi*pi)*dtz(i)*dtz(i)                  ! B
          else                       ! B = 8 pi^2 <u^2>
            dwz=dtz(i)                         ! B
            dth(i)=sqrt(dth(i)*0.5d0)/(pi+pi)  ! sqrt(<u^2>)
            dtk(i)=sqrt(dtk(i)*0.5d0)/(pi+pi)  ! sqrt(<u^2>)
          endif
          ! dth(i) & dtk(i) will be used only in the subroutine asfcef.

          do j=1,nabr                   ! a <= 0 for positron
            c(j,i)=b(j,i)/(-16d0*pi*pi)
            a(j,i)=a(j,i)*sqrt(4d0*pi/(b(j,i)+dwz))*rel
            b(j,i)=-4d0*pi*pi/(b(j,i)+dwz)
          end do
          if (abs(dchg(i)) > 1d-4) dchg(i)=dchg(i)*1.88972619d0*rel
          ! imaginary part (single electron & phonon excitation)
          if (nabi >= 2) then
            do j=nabr+1,nabr+nabi
              c(j,i)=b(j,i)/(-16d0*pi*pi)
              a(j,i)=a(j,i)*sqrt(4d0*pi/b(j,i))*rat
              b(j,i)=-4d0*pi*pi/b(j,i)
            end do
          endif
        end do
        return
        end
!**********************************************************
!       a.s.f. coefficient
!**********************************************************
        subroutine asfcef(nv,nelm,igh,igk,ghx,ghy,gky &
                         ,mab,nabr,nabi,a,c,asf,dth,dtk,Ghk)
        implicit none
        integer,intent(in) :: igh(nv),igk(nv),nv,nelm,mab,nabr,nabi
        real(8), intent(in) :: dth(nelm),dtk(nelm),ghx,ghy,gky
        real(8), intent(in) :: a(mab,nelm),c(mab,nelm)
        real(8), intent(out) :: Ghk(nv), asf(nv, mab, nelm)

        integer :: i, j
        real(8) :: s(nv), sj(nv)

        s = (ghx*igh)**2 + (ghy*igh+gky*igk)**2
        do j=1,nelm
          sj=-0.5d0*((dth(j)*ghx*igh)**2+ &
                  (dth(j)*ghy*igh+dtk(j)*gky*igk)**2)
          ! real part
          do i=1,nabr                  ! a <= 0 for positron
            asf(:,i,j)=a(i,j)*exp(c(i,j)*s(:)+sj(:)) ! Debye-Waller
          end do
          ! imaginary part
          if (nabi >= 2) then
            do i=nabr+1,nabr+nabi
              asf(:,i,j)=a(i,j)*exp(c(i,j)*s(:)+sj(:)*0.5d0) ! TDS
            end do
          endif
        end do

        Ghk = sqrt(s) ! for ion

        if (nabi == 1) then
          asf(1,nabr+1,:)=a(nabr+1,:) ! V(imag)=a*V(real), b is dummy
        endif
        return
        end



        ! compute structure factor for scpot
        subroutine scpot_stf(natm,nv,nsg,ma,mb,na,nb,gh,gk,x,y,dx,dy,st)
        implicit none
        integer, intent(in) :: natm, nv,nsg,ma,mb,na,nb
        real(8), intent(in) :: gh(nv),gk(nv)
        real(8), intent(in) :: x(natm),y(natm),dx,dy
        complex(8), intent(out) :: st(nv, natm)
        integer :: j
        do j=1,natm
          call strfac(nsg,nv,ma,mb,na,nb,gh,gk,x(j),y(j),dx,dy,st(1, j))
        end do
        end subroutine

!**********************************************************
!       scattering potential
!**********************************************************
        subroutine scpot(mv,nv,v,vi,natm,nelm,ielm,z,zi,mab,nabr,nabi,asf,b,dchg,rmesh,Ghk,ocr,st)
        implicit none
        integer, intent(in) :: mv,nv,natm,nelm,mab,nabr,nabi
        real(8), intent(in) :: asf(mv,mab,nelm),b(mab,nelm),dchg(nelm),Ghk(nv)
        real(8), intent(in) :: ocr(natm),z(natm),zi,rmesh
        integer, intent(in) :: ielm(natm)
        complex(8), intent(in) :: st(nv, natm)
        complex(8), intent(out) :: v(nv),vi(nv)


        ! local variables
        integer :: i,j,k,ie
        real(8) :: ocrmesh,z1,z2,ex,aex
        real(8) :: tv(nv), tvi(nv)
        complex(8) :: ct
        ! let exp(uf) -> 0
        real(8), parameter :: uf=-70d0

        do j=1,natm
          ocrmesh=rmesh*ocr(j)
          ie=ielm(j)
          z1=zi-z(j)
          z2=z1*z1
          ! ion
          if (abs(dchg(ie)) > 1d-4) then
            z1=-abs(z1)
            tv(1)=ocrmesh*z1*dchg(ie)
            do k=2,nv
              ex=Ghk(k)*z1
              if (ex > uf) then
                ex=exp(ex)*ocrmesh/Ghk(k)*dchg(ie)
                tv(k)=ex
              else
                tv(k)=0d0
              endif
            end do
          else
            tv = 0d0
          endif

          tvi = 0d0
          ! imaginary part : vi
          if (nabi >= 2) then
            do i=nabr+1,nabr+nabi
              ex=b(i,ie)*z2
              if (ex > uf) then
                ex=exp(ex)*ocrmesh
                tvi=asf(1:nv,i,ie)*ex
              endif
            end do
          endif

          ! real part : v
          do i=1,nabr
            ex=b(i,ie)*z2
            if (ex > uf) then
              ex=exp(ex)*ocrmesh
              if (nabi == 1) then             ! asf <= 0 for positron
                do k=1,nv
                  aex = asf(k,i,ie)*ex
                  tv(k)=tv(k)+aex
                  ! imaginary part : vi ( |asf(1,nabr+1,ie)| times v) 
                  tvi(k)=tvi(k)+abs(asf(1,nabr+1,ie)*aex)
                end do
              else
                tv=tv+asf(1:nv,i,ie)*ex
              end if
            endif
          end do
          !v = v + tv*st(:,j)
          !vi = vi + tvi*st(:,j)*(0d0,1d0)
          do k=1,nv
            ct = st(k,j)
            v(k) = v(k) + tv(k) * ct
            vi(k) = vi(k) + dcmplx(0d0, tvi(k))*ct
          end do
        end do ! j=1,natm
        return
        end

