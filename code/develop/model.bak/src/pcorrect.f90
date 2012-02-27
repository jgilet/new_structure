subroutine pcorrect(pcorr) 
!----------------------------------------------------                   
  USE header
!     Correct the NH presure (add the correction pcorr to p)            
!     modified for periodicew bc                                        
!     --------------------------                                        
!                                                                       
      integer i,j,k 
      real*8 ::                                                 &
     &     pcorr(0:NI+1,0:NJ+1,0:NK+1)                                  
                                                                        
      do k=0,NK+1 
         do j=0,NJ+1 
            do i=0,NI+1 
               p(i,j,k) = p(i,j,k) + pcorr(i,j,k) 
            end do 
         end do 
      end do 
                                                                        
      return 

!commented by J.Wang, 06/29/2011
!It is not necessary to compute the following every time step.
!Error message should be added in the initionlization subroutine.
!     Boundary conditions                                               
!      if ((NJ.lt.3).or.(NK.lt.3)) then
!         write(6,*) 'NJ,NK too small, stop in pcorrect'
!         stop
!      end if
                                                                        
!     Extrapolation at solid boundaries                                 
      do k=1,NK 
         do i=1,NI 
            p(i,0,k)= 3.d0*(p(i,1,k) -p(i,2,k))+p(i,3,k) 
            p(i,NJ+1,k)= 3.d0*(p(i,NJ,k) -p(i,NJ-1,k))+p(i,NJ-2,k) 
         end do 
!     Periodic-ew boundaries
!J.Wang 06/29/2011
#ifdef periodic_ew
         do j=0,NJ+1 
            p(0,j,k)= p(NI,j,k) 
            p(NI+1,j,k)= p(1,j,k) 
         end do 
#endif
      end do 

      do j=0,NJ+1 
         do i=0,NI+1 
            p(i,j,NK+1)= 3.d0*(p(i,j,NK) -p(i,j,NK-1))+p(i,j,NK-2) 
            p(i,j,0)= 3.d0*(p(i,j,1) - p(i,j,2))+p(i,j,3) 
         end do 
      end do 
!                                                                       
      return 
      END                                           
