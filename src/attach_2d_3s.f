!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2018 Guido Dhondt
!              Copyright (C) 2020 Shunji Uno
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!     
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of 
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
      subroutine attach_2d_3s(pneigh,pnode,ratio,dist,xil,etl)
!
!     ataches node with coordinates in "pnode" to the face containing 
!     "nterms" nodes with coordinates in field "pneigh" (nterms < 9).
!     cave: the coordinates are stored in pneigh(1..3,*)
!
      implicit none
!
      integer i,j,k,imin,jmin,im,jm
!
      real*8 ratio(9),pneigh(3,9),pnode(3),a,xi(-1:1,-1:1),
     &  et(-1:1,-1:1),p(3),distmin,d1,dist,xil,etl,
     &  ximin,ximax,etmin,etmax
!
      intent(in) pneigh
!
      intent(inout) xil,etl,pnode,dist,ratio
!
      d1=1.d0
!
      xi(0,0)=0.d0
      et(0,0)=0.d0
!      call distattach_2d(xi(0,0),et(0,0),pneigh,pnode,a,p,
!     &     ratio,3)
      ratio(1)=0.d0
      ratio(2)=1.d0/2.d0
      ratio(3)=1.d0/2.d0
!
      p(1)=ratio(1)*pneigh(1,1)+ratio(2)*pneigh(1,2)+
     &   ratio(3)*pneigh(1,3)
      p(2)=ratio(1)*pneigh(2,1)+ratio(2)*pneigh(2,2)+
     &   ratio(3)*pneigh(2,3)
      p(3)=ratio(1)*pneigh(3,1)+ratio(2)*pneigh(3,2)+
     &   ratio(3)*pneigh(3,3)
!
      a=(pnode(1)-p(1))**2+(pnode(2)-p(2))**2+
     &   (pnode(3)-p(3))**2
!
      distmin=a
      imin=0
      jmin=0
!
      do k=1,8
!
!     initialisation
!
         d1=d1/10.d0
!     
         do i=-1,1
            do j=-1,1
               if((i.eq.0).and.(j.eq.0)) cycle
!
               xi(i,j)=xi(0,0)+i*d1
               et(i,j)=et(0,0)+j*d1
!
!              check whether inside the (-1,1)x(-1,1) domain
!
               if((xi(i,j).le.1.d0).and.
     &              (xi(i,j).ge.-1.d0).and.
     &              (et(i,j).le.1.d0).and.
     &              (et(i,j).ge.-1.d0)) then
!                  call distattach_2d(xi(i,j),et(i,j),pneigh,pnode,a,p,
!     &                 ratio,3)
                  if(xi(i,j)+et(i,j).le.0.d0) then
                     ratio(2)=(xi(i,j)+1.d0)/2.d0
                     ratio(3)=(et(i,j)+1.d0)/2.d0
                  else
                     ratio(2)=(1.d0-et(i,j))/2.d0
                     ratio(3)=(1.d0-xi(i,j))/2.d0
                  endif
                  ratio(1)=1.d0-ratio(2)-ratio(3)
!
                  p(1)=ratio(1)*pneigh(1,1)+ratio(2)*pneigh(1,2)+
     &               ratio(3)*pneigh(1,3)
                  p(2)=ratio(1)*pneigh(2,1)+ratio(2)*pneigh(2,2)+
     &               ratio(3)*pneigh(2,3)
                  p(3)=ratio(1)*pneigh(3,1)+ratio(2)*pneigh(3,2)+
     &               ratio(3)*pneigh(3,3)
!
                  a=(pnode(1)-p(1))**2+(pnode(2)-p(2))**2+
     &               (pnode(3)-p(3))**2
!     
!                 checking for smallest initial distance
!     
                  if(a.lt.distmin) then
                     distmin=a
                     imin=i
                     jmin=j
                  endif
               endif
!
            enddo
         enddo
!     
!     minimizing the distance from the face to the node
!     
         if(d1.lt.0.05d0) then
            ximin=xi(imin,jmin)-d1*20
            ximax=xi(imin,jmin)+d1*20
            etmin=et(imin,jmin)-d1*20
            etmax=et(imin,jmin)+d1*20
         else
            ximin=-1.d0
            ximax=1.d0
            etmin=-1.d0
            etmax=1.d0
         endif
         do
!     
!     exit if minimum found
!     
            if((imin.eq.0).and.(jmin.eq.0)) exit
!
!           new center of 3x3 matrix
!
            xi(0,0)=xi(imin,jmin)
            et(0,0)=et(imin,jmin)
!
            im=imin
            jm=jmin
!
            imin=0
            jmin=0
!     
            do i=-1,1
               do j=-1,1
                  if((i+im.lt.-1).or.(i+im.gt.1).or.
     &                 (j+jm.lt.-1).or.(j+jm.gt.1)) then
!
                     xi(i,j)=xi(0,0)+i*d1
                     et(i,j)=et(0,0)+j*d1
!
!              check whether inside the (-1,1)x(-1,1) domain
!
!                     if((xi(i,j).le.1.d0).and.
!     &                  (xi(i,j).ge.-1.d0).and.
!     &                  (et(i,j).le.1.d0).and.
!     &                  (et(i,j).ge.-1.d0)) then
                     if((xi(i,j).le.ximax).and.
     &                  (xi(i,j).ge.ximin).and.
     &                  (et(i,j).le.etmax).and.
     &                  (et(i,j).ge.etmin)) then
!                        call distattach_2d(xi(i,j),et(i,j),pneigh,
!     &                       pnode,a,p,ratio,3)
!
                        if(xi(i,j)+et(i,j).le.0.d0) then
                           ratio(2)=(xi(i,j)+1.d0)/2.d0
                           ratio(3)=(et(i,j)+1.d0)/2.d0
                        else
                           ratio(2)=(1.d0-et(i,j))/2.d0
                           ratio(3)=(1.d0-xi(i,j))/2.d0
                        endif
                        ratio(1)=1.d0-ratio(2)-ratio(3)
!
                        p(1)=ratio(1)*pneigh(1,1)+ratio(2)*pneigh(1,2)+
     &                     ratio(3)*pneigh(1,3)
                        p(2)=ratio(1)*pneigh(2,1)+ratio(2)*pneigh(2,2)+
     &                     ratio(3)*pneigh(2,3)
                        p(3)=ratio(1)*pneigh(3,1)+ratio(2)*pneigh(3,2)+
     &                     ratio(3)*pneigh(3,3)
!
                        a=(pnode(1)-p(1))**2+(pnode(2)-p(2))**2+
     &                     (pnode(3)-p(3))**2
!
!                       check for new minimum
!
                        if(a.lt.distmin) then
                           distmin=a
                           imin=i
                           jmin=j
                        endif
                     endif
!
                  endif
               enddo
            enddo
         enddo
      enddo
!
!      call distattach_2d(xi(0,0),et(0,0),pneigh,pnode,a,p,
!     &     ratio,3)
      if(xi(0,0)+et(0,0).le.0.d0) then
         ratio(2)=(xi(0,0)+1.d0)/2.d0
         ratio(3)=(et(0,0)+1.d0)/2.d0
      else
         ratio(2)=(1.d0-et(0,0))/2.d0
         ratio(3)=(1.d0-xi(0,0))/2.d0
      endif
      ratio(1)=1.d0-ratio(2)-ratio(3)
!
      p(1)=ratio(1)*pneigh(1,1)+ratio(2)*pneigh(1,2)+
     &   ratio(3)*pneigh(1,3)
      p(2)=ratio(1)*pneigh(2,1)+ratio(2)*pneigh(2,2)+
     &   ratio(3)*pneigh(2,3)
      p(3)=ratio(1)*pneigh(3,1)+ratio(2)*pneigh(3,2)+
     &   ratio(3)*pneigh(3,3)
!
      a=(pnode(1)-p(1))**2+(pnode(2)-p(2))**2+
     &   (pnode(3)-p(3))**2
!
      do i=1,3
        pnode(i)=p(i)
      enddo
!
      dist=dsqrt(a)
!
      if(xi(0,0)+et(0,0).le.0.d0) then
         xil=(xi(0,0)+1.d0)/2.d0
         etl=(et(0,0)+1.d0)/2.d0
      else
         xil=(1.d0-et(0,0))/2.d0
         etl=(1.d0-xi(0,0))/2.d0
      endif
!     
c      write(*,*) 'end attach_2d.f'
      return
      end
