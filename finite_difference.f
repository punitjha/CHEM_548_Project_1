      program finite_difference
c 	  *********************************************************************
c 	  this program solves the 1D Schrodinger equation for an arbitrary potential using the finite difference method.
c 	  LAPACK and BLAS routines are used to calculate the eigenvectors and eigevalues of the tridiagoanal matrix that is obtained 
c 	  the boundary conditions are coded in the finite-difference derivative, i.e the first and last terms of the mesh are dropped.
c  	  DSTEVX computes selected eigenvalues and, optionally, eigenvectors of a real symmetric tridiagonal matrix A.  Eigenvalues and
c  	  eigenvectors can be selected by specifying either a range of values c  or a range of indices for the desired eigenvalues.
c 	  ************************************************************************
      implicit none
      integer    points 
	  integer    states
      parameter (points=500, states=10) !setting the values
      integer    i,n,p   ! variable for looping
      real*8     xxx,fc  !stores the distance between two mesh points 
	  real*8     potent1(points), potent2(points) !the potential function array
	  real*8     pmin ! setting the minimum potential for plotting purposes 
      real*8     diag1(points),diag2(points) ! the single dimensional array  for the diagonal elements 
      real*8     subd1(points),subd2(points) ! the single dimensional array  for the diagonal elements 
      real*8     energy1(points),energy2(points)  ! array to store the energies /eigenvalues 
      real*8     vector1(points,states),vector2(points,states) ! ! array to store the eigenvectors for each state
	  real*8     de ! this is used in plot settings to scale the graph
      integer    m !required by the LAPACK routine to pass out the no. of eigenvalues calculated
      integer    iwork(5*points) !required by the LAPACK routine 
      integer    ifail(points) !required by the LAPACK routine 
      integer    info !required by the LAPACK routine to check if the diagonalization was successfull or not 
      real*8     abstol !required by the LAPACK routine to check for the convergence of the values 
      real*8     work(5*points) !required by the LAPACK routine 
      real*8     dlamch !required by the LAPACK routine 
c	  ********************************************************
c 	  the units are changed to atomic units and the schrodinger equation takes the form as taught in the class
c 	  we create the mesh points now.
c	  *****************************************************************************************************************************
      xxx=20d0/dble(points-1) ! this assumes that the particle is in a box of lenght 20 units
c     ***************************************************************************************************************************************************

      
c       ALL THE POTENTIALS USED IN THE PROJECT ARE DEFINED BELOW AND EXCEPT FOR THE MORSE POTENTIAL ALL OTHERS ARE COMMMENTED OUT


c	  **************************************************************
c	  defining the potential morse potential for the ground state  *
c     **************************************************************
      write(*,'("Morse Potential and Frank-Condon Principle")')
      do i=1,points
      potent1(i)=5d0*(1-exp(-0.4*(dble(i-dble(points/3))*xxx)))**2
      enddo




c	  **************************************************************
c	  defining the potential morse potential for the excited state *
c     **************************************************************
      do i=1,points
      potent2(i)=4d0*(1-exp(-0.4*(dble(i-dble(1.5*points/3))*xxx)))**2+6
      enddo



c	  ********************************************
c	  defining the harmonic osillator potential  *
c     ********************************************
c      write(*,'("harmonic osillator potential")')
c      do i=1,points
c      potent(i)=0.5d0*(exp((dble(i-1)-dble((points-1)/2))*xxx))**2
c      enddo



c	  *********************************************
c	  defining the particle finite well potential *
c     *********************************************
c      write(*,'("particle finite well potential")')
c      do i=1,points
c      if((i.le.50).or.(i.ge.(points-50))) then
c      potent(i)=7
c      else 
c	   potent(i)=0
c      end if 
c      enddo


c	  ******************************************
c	  defining the particle in a box potential *
c     ******************************************

c      write(*,'("particle in a box")')
c      do i=1,points
c      potent(i)=0d0
c      enddo



c	  ***************************************************************************
c	  defining the particle in a box with a rectangular potential in the middle *
c     ***************************************************************************
c      write(*,'("particle in a box with rectangular potential in the middle")')
c      do i=1,points
c      if((i.ge.90).and.(i.le.110)) then
c      potent(i)=60d0*xxx
c      else 
c      potent(i)=0d0
c      end if 
c      enddo



c	  **********************************************
c	  this is the arbitrary potential defined here *
c     **********************************************
c      write(*,'("arbitrary potential")')
c      do i=1,points
c	  if((i.le.(points/2))) then
c      potent(i)=dble(i)*xxx
c      else 
c      potent(i)=dble((points-i))*xxx
c      end if 
c      enddo


c     ********************************************************************************************************************************************************
c     the finite difference mehtod is now used to discretize the K.E operator 
c     the boundary conditions are included via dropping the psi terms outside the mess
      do i=1,points
        diag1(i)=potent1(i)+(2d0/xxx**2)/2
        subd1(i)=(-1d0/xxx**2)/2
      enddo
c	  ********************************************************
c	  ********************************************************
c     the finite difference mehtod is now used to discretize the K.E operator 
c     the boundary conditions are included via dropping the psi terms outside the mess
      do i=1,points
        diag2(i)=potent2(i)+(2d0/xxx**2)/2
        subd2(i)=(-1d0/xxx**2)/2
      enddo
c	  ********************************************************
c	  we now call the LAPACK pakage to diagonalise our tridiagonal matrix  
      abstol=2d0*dlamch('s')    ! this is used by the LAPACK subroutine dstevx and is the The absolute error tolerance for the eigenvalues.
      call dstevx('v','i',points,diag1,subd1, 0d0,1d0,1,states,abstol,
     &            m,energy1,vector1,points,work,iwork,ifail,info)   !data separated into two lines are .f files cannot take more than 72 columns input --fixed source
      if(info.ne.0) write(*,*)"diagonalization failed !!"
c	  ********************************************************
c	  ********************************************************
c	  we now call the LAPACK pakage to diagonalise our tridiagonal matrix  
      abstol=2d0*dlamch('s')    ! this is used by the LAPACK subroutine dstevx and is the The absolute error tolerance for the eigenvalues.
      call dstevx('v','i',points,diag2,subd2, 0d0,1d0,1,states,abstol,
     &            m,energy2,vector2,points,work,iwork,ifail,info)   !data separated into two lines are .f files cannot take more than 72 columns input --fixed source
      if(info.ne.0) write(*,*)"diagonalization failed !!"
c	  ********************************************************
c       prinitng the eigenenergies to the terminal
      write(*,'("eigenenergies_1:")')
      do n=1,states
        write(*,'(I4,F30.20)') n,energy1(n)
      enddo
c	  ********************************************************
c	  ********************************************************
c       prinitng the eigenenergies to the terminal
      write(*,'("eigenenergies_2:")')
      do n=1,states
        write(*,'(I4,F30.20)') n,energy2(n)
      enddo
c	  ********************************************************
c	  ********************************************************
c       prinitng the difference in enegergies.
      write(*,'("eigenenergies_differences:")')
      do n=1,states
      do   p=1,states
      write(*,'("v''''=",I4," to v''= ",I4," Energy diffference is
     .        : ",F30.25)')n-1,p-1,(energy2(p)-energy1(n))
c        write(*,'(I4,F30.20)') n,energy2(n)
      enddo
      enddo
c	  ********************************************************
c	  ********************************************************
      write(*,*)"Calculating the Frank-Condon Factors"
      do i=1,states
       do p=1,states
         fc=0
        do n=1,points
         fc=fc+(vector1(n,i)*vector2(n,p))
          enddo
      write(*,'("v''''=",I4," to v''= ",I4," Frank-Condon Coefficient is
     .        : ",F30.25)')i-1,p-1,fc**2
      open(12,file='bar.data')
      write(12,'(I2,F30.20,F30.25)')i-1,
     .      (energy2(p)-energy1(i)),fc**2
        enddo
      enddo
c     *************************************************************
c     generating the output into a file to be plotted by gnu plot
c     setting the minimum of potential for plotting range
c	  ********************************************************
      pmin=1d10
      do i=1,points
        if(potent1(i).lt.pmin) pmin=potent1(i)
      enddo
c	  ********************************************************
c     setting the average setting of energy levels for plot scalling
      de=(energy2(states)-energy2(1))/dble(states-1)
      open(10,file='plot_1.data')
      write(10,'("set yrange [",F13.8,":",F13.8,"]")')
     .      pmin,energy2(m)+de  !setting the y range while plotting 
c	  ********************************************************
c     writing the plot commands for potential and eigenfunctions
      de=0.15d0*sqrt(dble(points))*de
      write(10,'("plot ""plot_1.data"" index 0 u 1:2 w lp, ", A1)')
     .      char(92) !char(92) is recognized by fortran 77 as the \ to be used in gnuplot
      write(10,'("""plot_1.data"" index 1 u 1:2 w lp, ", A1)')
     .      char(92) !char(92) is recognized by fortran 77 as the \ to be used in gnuplot

      do n=1,states
      write(10,'("""plot_1.data"" index",I4," u 1:(",F12.8,"+",
     .      F12.8,"*$2) w l ,",A1)') n+1,energy1(n),de,char(92)
      enddo
c     **********************************************************************
      do n=1,states-1
      write(10,'("""plot_1.data"" index",I4," u 1:(",F12.8,"+",
     .      F12.8,"*$2) w l ,",A1)') n+states+1,energy2(n),de,char(92)
      enddo
c     ****************************************************************** !written separately since the last one does not have a / character
      write(10,'("""plot_1.data"" index",I4," u 1:(",F12.8,"+",
     .      F12.8,"*$2) w l")')  2*states+1,energy2(states),de   
c     writing the data to the file
c     including the potential in the file the comments in gnuplot start with #
      write(10,'("# potential 1")')   
      do i=1,points
        write(10,'(I6,E20.10)') i-1,potent1(i)
      enddo
      write(10,*) ! setting two space in the output file necessary for the data to be read by gnuplot
      write(10,*)
      write(10,'("# potential 2")')   
      do i=1,points
        write(10,'(I6,E20.10)') i-1,potent2(i)
      enddo
      write(10,*) ! setting two space in the output file necessary for the data to be read by gnuplot
      write(10,*)
c     printing the eigenstates to the file 
      do n=1,states
        do i=1,points
          write(10,'(I6,E20.10)') i-1,vector1(i,n)
        enddo
        write(10,'(I6,E20.10)') 0,vector1(1,n)  !this is for the straight line
        write(10,*)   ! setting two space in the output file necessary for the data to be read by gnuplot
        write(10,*)
      enddo
c     printing the eigenstates to the file 
      do n=1,states
        do i=1,points
          write(10,'(I6,E20.10)') i-1,vector2(i,n)
        enddo
        write(10,'(I6,E20.10)') 0,vector2(1,n)
        write(10,*)   ! setting two space in the output file necessary for the data to be read by gnuplot
        write(10,*)
      enddo
c	  ********************************************************
      end
c	  ********************************************************
c     linking the two LAPACK packages to the program
      include 'dstevx.f' 
      include 'dstevx-blas.f'

