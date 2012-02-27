MODULE grids
USE header, ONLY : NI,NJ,NK,yc,xc,dx,dy
!*** grids
  REAL*8 :: xi_c(0:NI+1),et_c(0:NJ+1),sg_c(0:NK+1)
  REAL*8 :: xi_f(0:NI),et_f(0:NJ),sg_f(0:NK)
  real*8 :: dyM(0:NJ+1)

CONTAINS
  SUBROUTINE ini_grids()
    IMPLICIT NONE
    INTEGER :: i,j,k

    DO i = 0, NI+1
       xi_c(i) = REAL(i)
    ENDDO
    xi_f(:) = xi_c(0:NI)+0.5d0
    DO i = 0, NJ+1
       et_c(i) = REAL(i)
    ENDDO
    et_f(:) = et_c(0:NJ)+0.5d0
    DO i = 0, Nk+1
       sg_c(i) = REAL(i)
    ENDDO
    sg_f(:) = sg_c(0:NK)+0.5d0



      yc(0) = -0.5*dy*1.d-3 
      do j=1,NJ+1 
         yc(j)= yc(j-1) + dy*1.d-3 
      end do 
      xc(0) = -0.5*dx*1.d-3 
      do i=1,NI+1 
         xc(i)= xc(i-1) + dx*1.d-3 
      end do 





  END SUBROUTINE ini_grids


END MODULE grids
