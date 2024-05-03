      program AssociationKinetics
      implicit none
c>>  v2 is without oriental costrain because distance and orientation are correlated
      integer ID_num,ID_num2
      parameter (ID_num=100)

      integer nres_A_max
      parameter (nres_A_max=5000)
      integer nres_B_max
      parameter (nres_B_max=5000)
      integer natom_A_max
      parameter (natom_A_max=50000)
      integer natom_B_max
      parameter (natom_B_max=50000)
      integer nres_A,nres_B
      integer natom_A,natom_B
c      real*8 cell_range_x
c      parameter (cell_range_x=100.0)
c      real*8 cell_range_y
c      parameter (cell_range_y=100.0)
c      real*8 cell_range_z
c      parameter (cell_range_z=100.0)
      real*8 initial_clash
      parameter (initial_clash=10.0)
      real*8 dr,pai
      parameter (dr=3.8)
      parameter (pai=3.1415926)  
      integer simu_step
      parameter (simu_step=1000)
      integer num_trajec
      parameter (num_trajec=1000)
      real*8 temp_ini,temp_fin
      parameter (temp_ini=1.0,temp_fin=1.0) 
      real*8 time_step
      parameter (time_step=0.01) ! unit=1ns
      real*8 diffusion_constant_A,rot_D_A
c      parameter (diffusion_constant_A=11.5) ! unit=Angstron square per ns
c      parameter (rot_D_A=4.3) ! unit=degree per ns
      real*8 diffusion_constant_B,rot_D_B
c      parameter (diffusion_constant_B=12.2)
c      parameter (rot_D_B=5.2)
      real*8 dist_cutoff
      parameter (dist_cutoff=5.0)
      real*8 well_width
      parameter (well_width=2.0)
      real*8 k_constrain
      parameter (k_constrain=0.005)
      real*8 GO_depth
      parameter (GO_depth=5.0)
      real*8 rc
      parameter (rc=13.0)
      real*8 ElectricConstant
      parameter (ElectricConstant=562.0)
      real*8 DebyeLength
      real*8 W_S,W_E
      parameter (W_S=0.04,W_E=1.0)
cc
      integer output_flag_rec,output_flag_str,output_flag_trjrec
      integer output_flag_trjlog
      parameter (output_flag_rec=1,output_flag_str=0) 
      parameter (output_flag_trjrec=1,output_flag_trjlog=0)
      integer trj_output_freq
      parameter (trj_output_freq=100)
cc
      integer RB_res_num
      parameter (RB_res_num=4)
      real*8 dist_disturb
c      parameter (dist_disturb=30.0)
      real*8 angle_disturb
c      parameter (angle_disturb=30.0)
      real*8 dist_constrain
c      parameter (dist_constrain=15.0)
c      real*8 bond_thetapd,bond_thetapd_cutoff
c      parameter (bond_thetapd=180.0,bond_thetapd_cutoff=90.0)
c      real*8 bond_thetaot,bond_thetaot_cutoff
c      parameter (bond_thetaot=90.0,bond_thetaot_cutoff=90.0)
      integer numchainA_max,numchainB_max
      parameter (numchainA_max=10,numchainB_max=10)
ccc
      integer numchainA(ID_num)
      integer numchainB(ID_num)
      character*4 pdbid(ID_num)
      character*1 chainid_A(numchainA_max,ID_num)
      character*1 chainid_B(numchainB_max,ID_num)
      character*1 chainid_A2(numchainA_max)
      character*1 chainid_B2(numchainB_max)
      real*8 diffusion_constant_A_matr(ID_num),rot_D_A_matr(ID_num)
      real*8 diffusion_constant_B_matr(ID_num),rot_D_B_matr(ID_num)
      real*8 DebyeLength_matr(ID_num)
      real*8 IC_matr(ID_num)
      real*8 dist_constrain_matr(ID_num)
      real*8 angle_disturb_matr(ID_num)
      character*6 recna_A(natom_A_max),recna_B(natom_B_max)
      integer atno_A(natom_A_max),atno_B(natom_B_max)
      character*4 atna_A(natom_A_max),atna_B(natom_B_max)
      character*1 altloc_A(natom_A_max),altloc_B(natom_B_max)
      character*3 resna_A(natom_A_max),resna_B(natom_B_max)
      integer resseq_A(natom_A_max),resseq_B(natom_B_max)
      real*8 x_A(natom_A_max),y_A(natom_A_max),z_A(natom_A_max)
      real*8 x_B(natom_B_max),y_B(natom_B_max),z_B(natom_B_max)
      character*1 atom_chainid_A(natom_A_max)
      character*1 atom_chainid_B(natom_B_max)
      character*1 res_chainid_A(nres_A_max)
      character*1 res_chainid_B(nres_B_max)
      character*3 res_na_A(nres_A_max)
      character*3 res_na_B(nres_B_max)
      real*8 res_x_A(nres_A_max),res_y_A(nres_A_max),res_z_A(nres_A_max)
      real*8 res_x_B(nres_B_max),res_y_B(nres_B_max),res_z_B(nres_B_max)
      integer res_seq_A(nres_A_max),res_seq_B(nres_B_max)
      real*8 res_xsc_A(nres_A_max),res_ysc_A(nres_A_max)
      real*8 res_zsc_A(nres_A_max)
      real*8 res_xsc_B(nres_B_max),res_ysc_B(nres_B_max)
      real*8 res_zsc_B(nres_B_max)
      real*8 res_x_A1(nres_A_max),res_y_A1(nres_A_max)
      real*8 res_z_A1(nres_A_max)
      real*8 res_x_B1(nres_B_max),res_y_B1(nres_B_max)
      real*8 res_z_B1(nres_B_max)
      real*8 res_x_A0(nres_A_max),res_y_A0(nres_A_max)
      real*8 res_z_A0(nres_A_max)
      real*8 res_x_B0(nres_B_max),res_y_B0(nres_B_max)
      real*8 res_z_B0(nres_B_max)
      real*8 res_x_A_old(nres_A_max),res_y_A_old(nres_A_max)
      real*8 res_z_A_old(nres_A_max)
      real*8 res_x_B_old(nres_B_max),res_y_B_old(nres_B_max)
      real*8 res_z_B_old(nres_B_max)
      real*8 res_x_A_new(nres_A_max),res_y_A_new(nres_A_max)
      real*8 res_z_A_new(nres_A_max)
      real*8 res_x_B_new(nres_B_max),res_y_B_new(nres_B_max)
      real*8 res_z_B_new(nres_B_max)
      real*8 res_x_A_lcts(nres_A_max)
      real*8 res_y_A_lcts(nres_A_max)
      real*8 res_z_A_lcts(nres_A_max)
      real*8 res_x_B_lcts(nres_B_max)
      real*8 res_y_B_lcts(nres_B_max)
      real*8 res_z_B_lcts(nres_B_max)
      real*8 res_xsc_A1(nres_A_max),res_ysc_A1(nres_A_max)
      real*8 res_zsc_A1(nres_A_max)
      real*8 res_xsc_B1(nres_B_max),res_ysc_B1(nres_B_max)
      real*8 res_zsc_B1(nres_B_max)
      real*8 res_xsc_A0(nres_A_max),res_ysc_A0(nres_A_max)
      real*8 res_zsc_A0(nres_A_max)
      real*8 res_xsc_B0(nres_B_max),res_ysc_B0(nres_B_max)
      real*8 res_zsc_B0(nres_B_max)
      real*8 res_xsc_A_old(nres_A_max)
      real*8 res_ysc_A_old(nres_A_max)
      real*8 res_zsc_A_old(nres_A_max)
      real*8 res_xsc_B_old(nres_B_max)
      real*8 res_ysc_B_old(nres_B_max)
      real*8 res_zsc_B_old(nres_B_max)
      real*8 res_xsc_A_new(nres_A_max)
      real*8 res_ysc_A_new(nres_A_max)
      real*8 res_zsc_A_new(nres_A_max)
      real*8 res_xsc_B_new(nres_B_max)
      real*8 res_ysc_B_new(nres_B_max)
      real*8 res_zsc_B_new(nres_B_max)
      real*8 res_xsc_A_lcts(nres_A_max)
      real*8 res_ysc_A_lcts(nres_A_max)
      real*8 res_zsc_A_lcts(nres_A_max)
      real*8 res_xsc_B_lcts(nres_B_max)
      real*8 res_ysc_B_lcts(nres_B_max)
      real*8 res_zsc_B_lcts(nres_B_max)
      real*8 res_x_A_fin(nres_A_max,num_trajec)
      real*8 res_y_A_fin(nres_A_max,num_trajec)
      real*8 res_z_A_fin(nres_A_max,num_trajec)
      real*8 res_x_B_fin(nres_B_max,num_trajec)
      real*8 res_y_B_fin(nres_B_max,num_trajec)
      real*8 res_z_B_fin(nres_B_max,num_trajec)
      real*8 PB_x,PB_y,PB_z
      integer pair_list(natom_B_max,2)
      integer pair_num
      real*8 pair_dist(natom_B_max)
      integer respair_list(nres_B_max,2)
      integer respair_num
      real*8 respair_dist(nres_B_max)
      integer atompr(2,1000)
      real*8 xa(1000),ya(1000),za(1000),xb(1000),yb(1000),zb(1000)
      integer npair,npair_a,npair_b
      real*8 RMSD
      integer native_contact
      real*8 native_contact_ratio
      real*8 SolvantEffect
      real*8 ElectricPotential
      real*8 charge_A,charge_B
ccc
      real*8 simu_temp     
      integer mc_time_step
      integer n_t,id
      integer i,j,k,i2,j2
      real*8 temp_i,temp_j,temp_k
      real*8 theta,phi,psai,phai
      real*8 t(3,3)
      real*8 cm1_a_x
      real*8 cm1_a_y
      real*8 cm1_a_z
      real*8 cm1_b_x
      real*8 cm1_b_y
      real*8 cm1_b_z
      real*8 dist
      real*8 ene_vdw,ene_bond
      real*8 rij2,inverser,inverser6
      real*8 sum_diff,sum_diff2
      real*8 distance_step
      integer index
      integer pair_flag,pair_index
      real*8 ene_NativeContact
      integer formDimer_flag
      real*8 formDimer_time
      integer record_formDimer_flag(num_trajec)
      real*8 record_formDimer_time(num_trajec)
      real*8 ene_GO
      real*8 rij
      real*8 buribility_x_A(nres_A_max)
      real*8 buribility_y_A(nres_A_max)
      real*8 buribility_z_A(nres_A_max)
      real*8 buribility_r_A(nres_A_max)
      real*8 buribility_A(nres_A_max)
      real*8 buribility_x_B(nres_B_max)
      real*8 buribility_y_B(nres_B_max)
      real*8 buribility_z_B(nres_B_max)
      real*8 buribility_r_B(nres_B_max)
      real*8 buribility_B(nres_B_max)
      integer sc_num
      integer rescontact_flag
      character*17 serialnum
      real*8 dist_AB_COM
cc
      real*8 RB_A_x(RB_res_num)
      real*8 RB_A_y(RB_res_num)
      real*8 RB_A_z(RB_res_num)
      real*8 RB_B_x(RB_res_num)
      real*8 RB_B_y(RB_res_num)
      real*8 RB_B_z(RB_res_num)
      real*8 RB_A_x_00(RB_res_num)
      real*8 RB_A_y_00(RB_res_num)
      real*8 RB_A_z_00(RB_res_num)
      real*8 RB_A_x_01(RB_res_num)
      real*8 RB_A_y_01(RB_res_num)
      real*8 RB_A_z_01(RB_res_num)
      real*8 RB_B_x_00(RB_res_num)
      real*8 RB_B_y_00(RB_res_num)
      real*8 RB_B_z_00(RB_res_num)
      real*8 RB_B_x_01(RB_res_num)
      real*8 RB_B_y_01(RB_res_num)
      real*8 RB_B_z_01(RB_res_num)
      real*8 dist_cm1_IJ
      real*8 dist_temp,dist_vect
      real*8 vect_x,vect_y,vect_z
      real*8 rand_x,rand_y,rand_z
      integer clash_flag
      real*8 point_x(3),point_y(3),point_z(3)
      real*8 theta_pd,theta_ot
      integer flag
      character*100 fnam
      character*11 filename
      character*3 flag_end

      real rand3
      double precision r3            
      real rand4
      double precision r4
      real rand5
      double precision r5

      r3=5.0
      r4=5.0
      r5=5.0    

cccccccccccccccccccccccccccccccccccccccccccccccc

      do i=1,ID_num
         do j=1,numchainA_max
            chainid_A(j,i)='0'
         enddo
         do j=1,numchainB_max
            chainid_B(j,i)='0'
         enddo
      enddo

      do i=1,ID_num
         diffusion_constant_A_matr(i)=0
         rot_D_A_matr(i)=0
         diffusion_constant_B_matr(i)=0
         rot_D_B_matr(i)=0
         DebyeLength_matr(i)=0
         angle_disturb_matr(i)=0
      enddo

      call getarg(1,fnam)
      i=len_trim(fnam)
      filename=fnam(1:i)

      open(unit=10,file=filename//'.txt',status='old')
      index=0
      do i=1,ID_num
         read(10,*) flag_end
         if(flag_end.ne.'END')then
            index=index+1
         elseif(flag_end.eq.'END')then
            goto 200
         endif
      enddo
 200  continue
      close (10)
      
      ID_num2=index

      open(unit=10,file=filename//'.txt',status='old')
      do i=1,ID_num2
         read(10,*) pdbid(i),
     &        numchainA(i),(chainid_A(j,i),j=1,numchainA(i)),
     &        numchainB(i),(chainid_B(j,i),j=1,numchainB(i)),
     &        IC_matr(i),
     &        dist_constrain_matr(i)
      enddo
      close(10)


cccccccccccccccccccccccccccc

      do id=1,ID_num2



ccccccccccccccccccccccccccccc


         DebyeLength=(0.304/((IC_matr(id))**(0.5)))*10.0

         dist_constrain=dist_constrain_matr(id)
         dist_disturb=dist_constrain*2
         angle_disturb=0.5*(dist_constrain+dist_disturb)

         do j=1,numchainA_max
            chainid_A2(j)=chainid_A(j,id)
         enddo
         do j=1,numchainB_max
            chainid_B2(j)=chainid_B(j,id)
         enddo

c>>>  read chain A
         nres_A=0
         do j=1,5000
            res_na_A(j)='   '
            res_x_A(j)=0
            res_y_A(j)=0
            res_z_A(j)=0
            res_seq_A(j)=0
            res_chainid_A(j)=' '
         enddo
         natom_A=0
         do j=1,50000
            resna_A(j)='   '
            atna_A(j)='    '
            x_A(j)=0
            y_A(j)=0
            z_A(j)=0
            resseq_A(j)=0
            atno_A(j)=0
            atom_chainid_A(j)=' '
         enddo
         flag=0
    
         call readpdb(pdbid(id),numchainA(id),
     &        chainid_A2,
     &        nres_A,res_na_A,res_x_A,res_y_A,
     &        res_z_A,res_seq_A,natom_A,atno_A,
     &        resna_A,atna_A,x_A,y_A,z_A,resseq_A,
     &        res_chainid_A,atom_chainid_A,
     &        flag)

c         print*,id,nres_A

         diffusion_constant_A=26.647/((real(nres_A)*0.12)**(1.0/3.0))
         rot_D_A=17.33/(real(nres_A)*0.12)

c>>>  read chain B
         nres_B=0
         do j=1,5000
            res_na_B(j)='   '
            res_x_B(j)=0
            res_y_B(j)=0
            res_z_B(j)=0
            res_seq_B(j)=0
            res_chainid_B(j)=' '
         enddo
         natom_B=0
         do j=1,50000
            resna_B(j)='   '
            atna_B(j)='    '
            x_B(j)=0
            y_B(j)=0
            z_B(j)=0
            resseq_B(j)=0
            atno_B(j)=0
            atom_chainid_B(j)=' '
         enddo
         flag=0
    
         call readpdb(pdbid(id),numchainB(id),
     &        chainid_B2,
     &        nres_B,res_na_B,res_x_B,res_y_B,
     &        res_z_B,res_seq_B,natom_B,atno_B,
     &        resna_B,atna_B,x_B,y_B,z_B,resseq_B,
     &        res_chainid_B,atom_chainid_B,
     &        flag)

c         print*,id,nres_B
         
         diffusion_constant_B=26.647/((real(nres_B)*0.12)**(1.0/3.0))
         rot_D_B=17.33/(real(nres_B)*0.12)


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c>>          check interacting residues
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         do i=1,nres_A
            res_xsc_A(i)=0
            res_ysc_A(i)=0
            res_zsc_A(i)=0
         enddo
         do i=1,nres_A
            if(res_na_A(i).eq.'GLY')then
               res_xsc_A(i)=res_x_A(i)
               res_ysc_A(i)=res_y_A(i)
               res_zsc_A(i)=res_z_A(i)
            else
               sc_num=0
               if(res_na_A(i).eq.'ARG')then
                  do j=1,natom_A
                     if((resseq_A(j).eq.res_seq_A(i)).and.
     &                    (atom_chainid_A(j).eq.res_chainid_A(i)))then
                        if((atna_A(j).eq.' NH1').OR.
     &                       (atna_A(j).eq.' NH2'))then
                           sc_num=sc_num+1
                           res_xsc_A(i)=res_xsc_A(i)+x_A(j)
                           res_ysc_A(i)=res_ysc_A(i)+y_A(j)
                           res_zsc_A(i)=res_zsc_A(i)+z_A(j)
                        endif
                     endif
                  enddo              
               elseif(res_na_A(i).eq.'LYS')then
                  do j=1,natom_A
                     if((resseq_A(j).eq.res_seq_A(i)).and.
     &                    (atom_chainid_A(j).eq.res_chainid_A(i)))then
                        if(atna_A(j).eq.' NZ ')then
                           sc_num=sc_num+1
                           res_xsc_A(i)=res_xsc_A(i)+x_A(j)
                           res_ysc_A(i)=res_ysc_A(i)+y_A(j)
                           res_zsc_A(i)=res_zsc_A(i)+z_A(j)
                        endif
                     endif
                  enddo          
               elseif(res_na_A(i).eq.'ASP')then
                  do j=1,natom_A
                     if((resseq_A(j).eq.res_seq_A(i)).and.
     &                    (atom_chainid_A(j).eq.res_chainid_A(i)))then
                        if((atna_A(j).eq.' OD1').OR.
     &                       (atna_A(j).eq.' OD2'))then
                           sc_num=sc_num+1
                           res_xsc_A(i)=res_xsc_A(i)+x_A(j)
                           res_ysc_A(i)=res_ysc_A(i)+y_A(j)
                           res_zsc_A(i)=res_zsc_A(i)+z_A(j)
                        endif
                     endif
                  enddo              
               elseif(res_na_A(i).eq.'GLU')then
                  do j=1,natom_A
                     if((resseq_A(j).eq.res_seq_A(i)).and.
     &                    (atom_chainid_A(j).eq.res_chainid_A(i)))then
                        if((atna_A(j).eq.' OE1').OR.
     &                       (atna_A(j).eq.' OE2'))then
                           sc_num=sc_num+1
                           res_xsc_A(i)=res_xsc_A(i)+x_A(j)
                           res_ysc_A(i)=res_ysc_A(i)+y_A(j)
                           res_zsc_A(i)=res_zsc_A(i)+z_A(j)
                        endif
                     endif
                  enddo             
               elseif(res_na_A(i).eq.'HIS')then
                  do j=1,natom_A
                     if((resseq_A(j).eq.res_seq_A(i)).and.
     &                    (atom_chainid_A(j).eq.res_chainid_A(i)))then
                        if(atna_A(j).eq.' NE2')then
                           sc_num=sc_num+1
                           res_xsc_A(i)=res_xsc_A(i)+x_A(j)
                           res_ysc_A(i)=res_ysc_A(i)+y_A(j)
                           res_zsc_A(i)=res_zsc_A(i)+z_A(j)
                        endif
                     endif
                  enddo                  
               else
                  do j=1,natom_A
                     if((resseq_A(j).eq.res_seq_A(i)).and.
     &                    (atom_chainid_A(j).eq.res_chainid_A(i)))then
                        if((atna_A(j).ne.' CA ')
     &                       .AND.(atna_A(j).ne.' N  ')
     &                    .AND.(atna_A(j).ne.' O  ')
     &                       .AND.(atna_A(j).ne.' C  ')
     &                       .AND.(atna_A(j).ne.' H  '))then
                           sc_num=sc_num+1
                           res_xsc_A(i)=res_xsc_A(i)+x_A(j)
                           res_ysc_A(i)=res_ysc_A(i)+y_A(j)
                           res_zsc_A(i)=res_zsc_A(i)+z_A(j)
                        endif
                     endif
                  enddo
               endif
               if(sc_num.ne.0)then
                  res_xsc_A(i)=res_xsc_A(i)/sc_num
                  res_ysc_A(i)=res_ysc_A(i)/sc_num
                  res_zsc_A(i)=res_zsc_A(i)/sc_num
               else
                  res_xsc_A(i)=res_x_A(i)
                  res_ysc_A(i)=res_y_A(i)
                  res_zsc_A(i)=res_z_A(i)
               endif
            endif
         enddo
         do i=1,nres_B
            res_xsc_B(i)=0
            res_ysc_B(i)=0
            res_zsc_B(i)=0
         enddo
         do i=1,nres_B
            if(res_na_B(i).eq.'GLY')then
               res_xsc_B(i)=res_x_B(i)
               res_ysc_B(i)=res_y_B(i)
               res_zsc_B(i)=res_z_B(i)
            else
               sc_num=0
               if(res_na_B(i).eq.'ARG')then
                  do j=1,natom_B
                     if((resseq_B(j).eq.res_seq_B(i)).and.
     &                    (atom_chainid_B(j).eq.res_chainid_B(i)))then
                        if((atna_B(j).eq.' NH1').OR.
     &                       (atna_B(j).eq.' NH2'))then
                           sc_num=sc_num+1
                           res_xsc_B(i)=res_xsc_B(i)+x_B(j)
                           res_ysc_B(i)=res_ysc_B(i)+y_B(j)
                           res_zsc_B(i)=res_zsc_B(i)+z_B(j)
                        endif
                     endif
                  enddo              
               elseif(res_na_B(i).eq.'LYS')then
                  do j=1,natom_B
                     if((resseq_B(j).eq.res_seq_B(i)).and.
     &                    (atom_chainid_B(j).eq.res_chainid_B(i)))then
                        if(atna_B(j).eq.' NZ ')then
                           sc_num=sc_num+1
                           res_xsc_B(i)=res_xsc_B(i)+x_B(j)
                           res_ysc_B(i)=res_ysc_B(i)+y_B(j)
                           res_zsc_B(i)=res_zsc_B(i)+z_B(j)
                        endif
                     endif
                  enddo          
               elseif(res_na_B(i).eq.'ASP')then
                  do j=1,natom_B
                     if((resseq_B(j).eq.res_seq_B(i)).and.
     &                    (atom_chainid_B(j).eq.res_chainid_B(i)))then
                        if((atna_B(j).eq.' OD1').OR.
     &                       (atna_B(j).eq.' OD2'))then
                           sc_num=sc_num+1
                           res_xsc_B(i)=res_xsc_B(i)+x_B(j)
                           res_ysc_B(i)=res_ysc_B(i)+y_B(j)
                           res_zsc_B(i)=res_zsc_B(i)+z_B(j)
                        endif
                     endif
                  enddo              
               elseif(res_na_B(i).eq.'GLU')then
                  do j=1,natom_B
                     if((resseq_B(j).eq.res_seq_B(i)).and.
     &                    (atom_chainid_B(j).eq.res_chainid_B(i)))then
                        if((atna_B(j).eq.' OE1').OR.
     &                       (atna_B(j).eq.' OE2'))then
                           sc_num=sc_num+1
                           res_xsc_B(i)=res_xsc_B(i)+x_B(j)
                           res_ysc_B(i)=res_ysc_B(i)+y_B(j)
                           res_zsc_B(i)=res_zsc_B(i)+z_B(j)
                        endif
                     endif
                  enddo             
               elseif(res_na_B(i).eq.'HIS')then
                  do j=1,natom_B
                     if((resseq_B(j).eq.res_seq_B(i)).and.
     &                    (atom_chainid_B(j).eq.res_chainid_B(i)))then
                        if(atna_B(j).eq.' NE2')then
                           sc_num=sc_num+1
                           res_xsc_B(i)=res_xsc_B(i)+x_B(j)
                           res_ysc_B(i)=res_ysc_B(i)+y_B(j)
                           res_zsc_B(i)=res_zsc_B(i)+z_B(j)
                        endif
                     endif
                  enddo        
               else
                  do j=1,natom_B
                     if((resseq_B(j).eq.res_seq_B(i)).and.
     &                    (atom_chainid_B(j).eq.res_chainid_B(i)))then
                        if((atna_B(j).ne.' CA ')
     &                       .AND.(atna_B(j).ne.' N  ')
     &                       .AND.(atna_B(j).ne.' O  ')
     &                       .AND.(atna_B(j).ne.' C  ')
     &                       .AND.(atna_B(j).ne.' H  '))then
                           sc_num=sc_num+1
                           res_xsc_B(i)=res_xsc_B(i)+x_B(j)
                           res_ysc_B(i)=res_ysc_B(i)+y_B(j)
                           res_zsc_B(i)=res_zsc_B(i)+z_B(j)
                        endif
                     endif
                  enddo
               endif
               if(sc_num.ne.0)then
                  res_xsc_B(i)=res_xsc_B(i)/sc_num
                  res_ysc_B(i)=res_ysc_B(i)/sc_num
                  res_zsc_B(i)=res_zsc_B(i)/sc_num
               else
                  res_xsc_B(i)=res_x_B(i)
                  res_ysc_B(i)=res_y_B(i)
                  res_zsc_B(i)=res_z_B(i)
               endif
            endif
         enddo
        
         do i=1,nres_B
            do j=1,2
               respair_list(i,j)=0
            enddo
            respair_dist(i)=0
         enddo
         respair_num=0
         do i=1,nres_A
            do j=1,nres_B
               rescontact_flag=0
               do i2=1,natom_A
                  if((resseq_A(i2).eq.res_seq_A(i)).and.
     &                 (atom_chainid_A(i2).eq.res_chainid_A(i)))then
                     if((atna_A(i2).ne.' CA ').AND.
     &                    (atna_A(i2).ne.' N  ')
     &                    .AND.(atna_A(i2).ne.' O  ')
     &                    .AND.(atna_A(i2).ne.' C  ')
     &                    .AND.(atna_A(i2).ne.' H  ')
     &                    .AND.(atna_A(i2)(2:2).ne.'H'))then
                        do j2=1,natom_B
                           if((resseq_B(j2).eq.res_seq_B(j)).and.
     &                          (atom_chainid_B(j2).eq.res_chainid_B(j))
     &                          )then
                              if((atna_B(j2).ne.' CA ').AND.
     &                             (atna_B(j2).ne.' N  ')
     &                             .AND.(atna_B(j2).ne.' O  ')
     &                             .AND.(atna_B(j2).ne.' C  ')
     &                             .AND.(atna_B(j2).ne.' H  ')
     &                             .AND.(atna_B(j2)(2:2).ne.'H'))then
                                 dist=sqrt((x_A(i2)-x_B(j2))**2+
     &                                (y_A(i2)-y_B(j2))**2+
     &                                (z_A(i2)-z_B(j2))**2)
                                 if(dist.le.dist_cutoff)then
                                    rescontact_flag=1
                                 endif
                              endif
                           endif
                        enddo
                     endif
                  endif
               enddo
               if(rescontact_flag.eq.1)then
                  dist=sqrt((res_x_A(i)-res_x_B(j))**2+
     &                 (res_y_A(i)-res_y_B(j))**2+
     &                 (res_z_A(i)-res_z_B(j))**2)
                  respair_num=respair_num+1
                  respair_list(respair_num,1)=res_seq_A(i)
                  respair_list(respair_num,2)=res_seq_B(j)
                  respair_dist(respair_num)=dist
               endif
            enddo
         enddo
         
         do i=1,natom_B
            do j=1,2
               pair_list(i,j)=0
            enddo
            pair_dist(i)=0
         enddo
         pair_num=0
         do i=1,natom_A
            if((atna_A(i).ne.' CA ').AND.(atna_A(i).ne.' N  ')
     &           .AND.(atna_A(i).ne.' O  ')
     &           .AND.(atna_A(i).ne.' C  ')
     &           .AND.(atna_A(i).ne.' H  ')
     &           .AND.(atna_A(i)(2:2).ne.'H'))then
               do j=1,natom_B
                  if((atna_B(j).ne.' CA ').AND.(atna_B(j).ne.' N  ')
     &                 .AND.(atna_B(j).ne.' O  ')
     &                 .AND.(atna_B(j).ne.' C  ')
     &                 .AND.(atna_B(j).ne.' H  ')
     &                 .AND.(atna_B(j)(2:2).ne.'H'))then
                     dist=sqrt((x_A(i)-x_B(j))**2+
     &                    (y_A(i)-y_B(j))**2+
     &                    (z_A(i)-z_B(j))**2)
                     if(dist.le.dist_cutoff)then
                        pair_num=pair_num+1
                        pair_list(pair_num,1)=i
                        pair_list(pair_num,2)=j
                        pair_dist(pair_num)=dist
                     endif
                  endif
               enddo
            endif
         enddo
         
         index=0
         do i=1,nres_A
            buribility_x_A(i)=0
            buribility_y_A(i)=0
            buribility_z_A(i)=0
            buribility_r_A(i)=0
            buribility_A(i)=0
            do j=1,nres_A
               if(abs(i-j).gt.1)then
                  rij=sqrt((res_x_A(i)-res_x_A(j))**2+
     &                 (res_y_A(i)-res_y_A(j))**2+
     &                 (res_z_A(i)-res_z_A(j))**2)
                  if(rij.le.rc)then
                     buribility_x_A(i)=buribility_x_A(i)+
     &                    (res_x_A(j)-res_x_A(i))
                     buribility_y_A(i)=buribility_y_A(i)+
     &                    (res_y_A(j)-res_y_A(i))
                     buribility_z_A(i)=buribility_z_A(i)+
     &                    (res_z_A(j)-res_z_A(i))
                     buribility_r_A(i)=buribility_r_A(i)+rij
                  endif
               endif
            enddo
            buribility_A(i)=(sqrt(buribility_x_A(i)**2+
     &           buribility_y_A(i)**2+
     &           buribility_z_A(i)**2))/buribility_r_A(i)
         enddo
         
         index=0
         do i=1,nres_B
            buribility_x_B(i)=0
            buribility_y_B(i)=0
            buribility_z_B(i)=0
            buribility_r_B(i)=0
            buribility_B(i)=0
            do j=1,nres_B
               if(abs(i-j).gt.1)then
                  rij=sqrt((res_x_B(i)-res_x_B(j))**2+
     &                 (res_y_B(i)-res_y_B(j))**2+
     &                 (res_z_B(i)-res_z_B(j))**2)
                  if(rij.le.rc)then
                     buribility_x_B(i)=buribility_x_B(i)+
     &                    (res_x_B(j)-res_x_B(i))
                     buribility_y_B(i)=buribility_y_B(i)+
     &                    (res_y_B(j)-res_y_B(i))
                     buribility_z_B(i)=buribility_z_B(i)+
     &                    (res_z_B(j)-res_z_B(i))
                     buribility_r_B(i)=buribility_r_B(i)+rij
                  endif
               endif
            enddo
            buribility_B(i)=(sqrt(buribility_x_B(i)**2+
     &           buribility_y_B(i)**2+
     &           buribility_z_B(i)**2))/buribility_r_B(i)
         enddo
         
         SolvantEffect=0
         index=0
         do i=1,nres_A
            do j=1,nres_B
               pair_flag=0
               pair_index=0
               do i2=1,respair_num
                  if((respair_list(i2,1).eq.res_seq_A(i)).AND.
     &                 (respair_list(i2,2).eq.res_seq_B(j)))then
                     pair_flag=1
                     pair_index=i2
                  endif
               enddo                  
               if(pair_flag.eq.1)then                       
                  dist=sqrt((res_xsc_A(i)-res_xsc_B(j))**2+
     &                 (res_ysc_A(i)-res_ysc_B(j))**2+
     &                 (res_zsc_A(i)-res_zsc_B(j))**2)
                  
                  index=index+1

                  if(dist.le.7.5)then
                     

                     if(res_na_A(i).eq.'ALA')then
                        SolvantEffect=SolvantEffect-1.8
                     elseif(res_na_A(i).eq.'ILE')then
                        SolvantEffect=SolvantEffect-4.5
                     elseif(res_na_A(i).eq.'LEU')then
                        SolvantEffect=SolvantEffect-3.8
                     elseif(res_na_A(i).eq.'PHE')then
                        SolvantEffect=SolvantEffect-2.8
                     elseif(res_na_A(i).eq.'VAL')then
                        SolvantEffect=SolvantEffect-4.2
                     elseif(res_na_A(i).eq.'MET')then
                        SolvantEffect=SolvantEffect-1.9
                     elseif(res_na_A(i).eq.'CYS')then
                        SolvantEffect=SolvantEffect-2.5
                     endif
                     
                     if(res_na_B(j).eq.'ALA')then
                        SolvantEffect=SolvantEffect-1.8
                     elseif(res_na_B(j).eq.'ILE')then
                        SolvantEffect=SolvantEffect-4.5
                     elseif(res_na_B(j).eq.'LEU')then
                        SolvantEffect=SolvantEffect-3.8
                     elseif(res_na_B(j).eq.'PHE')then
                        SolvantEffect=SolvantEffect-2.8
                     elseif(res_na_B(j).eq.'VAL')then
                        SolvantEffect=SolvantEffect-4.2
                     elseif(res_na_B(j).eq.'MET')then
                        SolvantEffect=SolvantEffect-1.9
                     elseif(res_na_B(j).eq.'CYS')then
                        SolvantEffect=SolvantEffect-2.5
                     endif

                  endif
               endif
            enddo
         enddo
         
         
         ElectricPotential=0
         charge_A=0
         charge_B=0
         do i=1,nres_A
            if((res_na_A(i).eq.'ARG').OR.(res_na_A(i).eq.'LYS'))then
               charge_A=1.0
            elseif((res_na_A(i).eq.'GLU').OR.(res_na_A(i).eq.'ASP'))then
               charge_A=-1.0
            elseif(res_na_A(i).eq.'HIS')then
               charge_A=0.5
            else
               charge_A=0.0
            endif
            do j=1,nres_B
               if((res_na_B(j).eq.'ARG').OR.(res_na_B(j).eq.'LYS'))then
                  charge_B=1.0
               elseif((res_na_B(j).eq.'GLU').OR.(res_na_B(j).eq.'ASP')
     &                 )then
                  charge_B=-1.0
               elseif(res_na_B(j).eq.'HIS')then
                  charge_B=0.5
               else
                  charge_B=0.0
               endif
               pair_flag=0
               pair_index=0
               do i2=1,respair_num
                  if((respair_list(i2,1).eq.res_seq_A(i)).AND.
     &                 (respair_list(i2,2).eq.res_seq_B(j)))then
                     pair_flag=1
                     pair_index=i2
                  endif
               enddo                  
               if(pair_flag.eq.1)then        
                  
                  dist=sqrt((res_xsc_A(i)-res_xsc_B(j))**2+
     &                 (res_ysc_A(i)-res_ysc_B(j))**2+
     &                 (res_zsc_A(i)-res_zsc_B(j))**2)
                  ElectricPotential=ElectricPotential+
     &                 ElectricConstant*charge_A*charge_B/
     &                 (10.0*dist*DEXP(dist/DebyeLength))
               endif
            enddo
         enddo


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c>>  Multiple trajectories
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c>>
         do n_t=1,num_trajec

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c   construct the initial position of molecules in 3D
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c>>          center of mass of A and B
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            cm1_a_x=0
            cm1_a_y=0
            cm1_a_z=0
            cm1_b_x=0
            cm1_b_y=0
            cm1_b_z=0
            do i=1,nres_A
               cm1_a_x=cm1_a_x+res_x_A(i)
               cm1_a_y=cm1_a_y+res_y_A(i)
               cm1_a_z=cm1_a_z+res_z_A(i)
            enddo
            do i=1,nres_B
               cm1_b_x=cm1_b_x+res_x_B(i)
               cm1_b_y=cm1_b_y+res_y_B(i)
               cm1_b_z=cm1_b_z+res_z_B(i)
            enddo
            cm1_a_x=cm1_a_x/(real(nres_A))
            cm1_a_y=cm1_a_y/(real(nres_A))
            cm1_a_z=cm1_a_z/(real(nres_A))
            cm1_b_x=cm1_b_x/(real(nres_B))
            cm1_b_y=cm1_b_y/(real(nres_B))
            cm1_b_z=cm1_b_z/(real(nres_B))
            
            dist_cm1_IJ=sqrt(
     &           (cm1_b_x-cm1_a_x)**2+
     &           (cm1_b_y-cm1_a_y)**2+
     &           (cm1_b_z-cm1_a_z)**2)      
            
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c>>          calculate the X, Y and Z axes of proteins A and B
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
            RB_B_x(1)=cm1_b_x
            RB_B_y(1)=cm1_b_y
            RB_B_z(1)=cm1_b_z
            
            RB_A_x(1)=cm1_a_x
            RB_A_y(1)=cm1_a_y
            RB_A_z(1)=cm1_a_z
            
            RB_B_x(4)=(cm1_b_x+cm1_a_x)/2
            RB_B_y(4)=(cm1_b_y+cm1_a_y)/2
            RB_B_z(4)=(cm1_b_z+cm1_a_z)/2
            
            RB_A_x(4)=(cm1_b_x+cm1_a_x)/2
            RB_A_y(4)=(cm1_b_y+cm1_a_y)/2
            RB_A_z(4)=(cm1_b_z+cm1_a_z)/2
            
            dist_temp=sqrt((RB_B_x(4)-RB_B_x(1))**2+
     &           (RB_B_y(4)-RB_B_y(1))**2)
            RB_B_x(2)=cm1_b_x+(RB_B_y(4)-RB_B_y(1))
     &           *dist_cm1_IJ/(2*dist_temp)
            RB_B_y(2)=cm1_b_y+(RB_B_x(1)-RB_B_x(4))
     &           *dist_cm1_IJ/(2*dist_temp)
            RB_B_z(2)=cm1_b_z
            
            RB_A_x(3)=cm1_a_x-(RB_B_y(4)-RB_B_y(1))
     &           *dist_cm1_IJ/(2*dist_temp)
            RB_A_y(3)=cm1_a_y-(RB_B_x(1)-RB_B_x(4))
     &           *dist_cm1_IJ/(2*dist_temp)
            RB_A_z(3)=cm1_a_z
            
            vect_x=(RB_B_y(2)-cm1_b_y)*(RB_B_z(4)-cm1_b_z)
     &           -(RB_B_z(2)-cm1_b_z)*(RB_B_y(4)-cm1_b_y)
            vect_y=(RB_B_z(2)-cm1_b_z)*(RB_B_x(4)-cm1_b_x)
     &           -(RB_B_x(2)-cm1_b_x)*(RB_B_z(4)-cm1_b_z)
            vect_z=(RB_B_x(2)-cm1_b_x)*(RB_B_y(4)-cm1_b_y)
     &           -(RB_B_y(2)-cm1_b_y)*(RB_B_x(4)-cm1_b_x)
            dist_vect=sqrt(vect_x**2+vect_y**2+vect_z**2)
            
            RB_B_x(3)=cm1_b_x+vect_x
     &        *dist_cm1_IJ/(2*dist_vect)      
            RB_B_y(3)=cm1_b_y+vect_y
     &           *dist_cm1_IJ/(2*dist_vect)      
            RB_B_z(3)=cm1_b_z+vect_z
     &           *dist_cm1_IJ/(2*dist_vect)      
            
            RB_A_x(2)=cm1_a_x+vect_x
     &           *dist_cm1_IJ/(2*dist_vect)      
            RB_A_y(2)=cm1_a_y+vect_y
     &           *dist_cm1_IJ/(2*dist_vect)      
            RB_A_z(2)=cm1_a_z+vect_z
     &           *dist_cm1_IJ/(2*dist_vect)      
            
 100        continue
            
            do i=1,RB_res_num
               RB_B_x_01(i)=RB_B_x(i)           
               RB_B_y_01(i)=RB_B_y(i)           
               RB_B_z_01(i)=RB_B_z(i)            
            enddo
            do i=1,nres_B
               res_x_B0(i)=res_x_B(i)    
               res_y_B0(i)=res_y_B(i)    
               res_z_B0(i)=res_z_B(i)        
            enddo
            do i=1,nres_B
               res_xsc_B0(i)=res_xsc_B(i)    
               res_ysc_B0(i)=res_ysc_B(i)    
               res_zsc_B0(i)=res_zsc_B(i)        
            enddo         
            
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c>>         structural deviation of protein A
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c>>>>>   random position
            
            rand_x=2*rand3(r3)*dist_disturb-dist_disturb
            rand_y=2*rand3(r3)*dist_disturb-dist_disturb
            rand_z=2*rand3(r3)*dist_disturb-dist_disturb
            
            do i=1,RB_res_num
               RB_A_x_00(i)=RB_A_x(i)+rand_x           
               RB_A_y_00(i)=RB_A_y(i)+rand_y           
               RB_A_z_00(i)=RB_A_z(i)+rand_z            
            enddo
            do i=1,nres_A
               res_x_A1(i)=res_x_A(i)+rand_x    
               res_y_A1(i)=res_y_A(i)+rand_y    
               res_z_A1(i)=res_z_A(i)+rand_z        
            enddo
            do i=1,nres_A
               res_xsc_A1(i)=res_xsc_A(i)+rand_x    
               res_ysc_A1(i)=res_ysc_A(i)+rand_y    
               res_zsc_A1(i)=res_zsc_A(i)+rand_z        
            enddo
            
c>>>>>   random orientation

            theta=(2*rand3(r3)-1)*pai*angle_disturb/180.0
            phi=(2*rand3(r3)-1)*pai*angle_disturb/180.0
            psai=(2*rand3(r3)-1)*pai*angle_disturb/180.0
            
            t(1,1)=cos(psai)*cos(phi)-cos(theta)*sin(phi)*sin(psai)
            t(1,2)=-sin(psai)*cos(phi)-cos(theta)*sin(phi)*cos(psai)
            t(1,3)=sin(theta)*sin(phi)
            
            t(2,1)=cos(psai)*sin(phi)+cos(theta)*cos(phi)*sin(psai)
            t(2,2)=-sin(psai)*sin(phi)+cos(theta)*cos(phi)*cos(psai)
            t(2,3)=-sin(theta)*cos(phi)
            
            t(3,1)=sin(psai)*sin(theta)
            t(3,2)=cos(psai)*sin(theta)
            t(3,3)=cos(theta)
            
            RB_A_x_01(1)=RB_A_x_00(1) 
            RB_A_y_01(1)=RB_A_y_00(1)
            RB_A_z_01(1)=RB_A_z_00(1)
            do k=2,RB_res_num
               RB_A_x_01(k)=t(1,1)*(RB_A_x_00(k)-RB_A_x_01(1))
     &              +t(1,2)*(RB_A_y_00(k)-RB_A_y_01(1))
     &              +t(1,3)*(RB_A_z_00(k)-RB_A_z_01(1))
     &              +RB_A_x_01(1)
               RB_A_y_01(k)=t(2,1)*(RB_A_x_00(k)-RB_A_x_01(1))
     &              +t(2,2)*(RB_A_y_00(k)-RB_A_y_01(1))
     &              +t(2,3)*(RB_A_z_00(k)-RB_A_z_01(1))
     &              +RB_A_y_01(1)
               RB_A_z_01(k)=t(3,1)*(RB_A_x_00(k)-RB_A_x_01(1))
     &              +t(3,2)*(RB_A_y_00(k)-RB_A_y_01(1))
     &              +t(3,3)*(RB_A_z_00(k)-RB_A_z_01(1))
     &              +RB_A_z_01(1)
            enddo
            do k=1,nres_A
               res_x_A0(k)=t(1,1)*(res_x_A1(k)-RB_A_x_01(1))
     &              +t(1,2)*(res_y_A1(k)-RB_A_y_01(1))
     &              +t(1,3)*(res_z_A1(k)-RB_A_z_01(1))
     &              +RB_A_x_01(1)
               res_y_A0(k)=t(2,1)*(res_x_A1(k)-RB_A_x_01(1))
     &              +t(2,2)*(res_y_A1(k)-RB_A_y_01(1))
     &              +t(2,3)*(res_z_A1(k)-RB_A_z_01(1))
     &              +RB_A_y_01(1)
               res_z_A0(k)=t(3,1)*(res_x_A1(k)-RB_A_x_01(1))
     &              +t(3,2)*(res_y_A1(k)-RB_A_y_01(1))
     &              +t(3,3)*(res_z_A1(k)-RB_A_z_01(1))
     &              +RB_A_z_01(1)
            enddo
            do k=1,nres_A
               res_xsc_A0(k)=t(1,1)*(res_xsc_A1(k)-RB_A_x_01(1))
     &              +t(1,2)*(res_ysc_A1(k)-RB_A_y_01(1))
     &              +t(1,3)*(res_zsc_A1(k)-RB_A_z_01(1))
     &              +RB_A_x_01(1)
               res_ysc_A0(k)=t(2,1)*(res_xsc_A1(k)-RB_A_x_01(1))
     &              +t(2,2)*(res_ysc_A1(k)-RB_A_y_01(1))
     &              +t(2,3)*(res_zsc_A1(k)-RB_A_z_01(1))
     &              +RB_A_y_01(1)
               res_zsc_A0(k)=t(3,1)*(res_xsc_A1(k)-RB_A_x_01(1))
     &              +t(3,2)*(res_ysc_A1(k)-RB_A_y_01(1))
     &              +t(3,3)*(res_zsc_A1(k)-RB_A_z_01(1))
     &              +RB_A_z_01(1)
            enddo


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c>>    the initial configuration has to satisfy constrain
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c>>   1) no clash
        
            clash_flag=0
            
            do i=1,nres_A
               do j=1,nres_B
                  dist=sqrt((res_x_A0(i)-res_x_B0(j))**2+
     &                 (res_y_A0(i)-res_y_B0(j))**2+
     &                 (res_z_A0(i)-res_z_B0(j))**2)
                  if(dist.le.initial_clash)then
                     clash_flag=1
                  endif
               enddo
            enddo
            
            do i=1,nres_A
               do j=1,nres_B
                  dist=sqrt((res_xsc_A0(i)-res_xsc_B0(j))**2+
     &                 (res_ysc_A0(i)-res_ysc_B0(j))**2+
     &                 (res_zsc_A0(i)-res_zsc_B0(j))**2)
                  if(dist.le.initial_clash)then
                     clash_flag=1
                  endif
               enddo
            enddo
            
            do i=1,nres_A
               do j=1,nres_B
                  dist=sqrt((res_x_A0(i)-res_xsc_B0(j))**2+
     &                 (res_y_A0(i)-res_ysc_B0(j))**2+
     &                 (res_z_A0(i)-res_zsc_B0(j))**2)
                  if(dist.le.initial_clash)then
                     clash_flag=1
                  endif
               enddo
            enddo
            
            do i=1,nres_A
               do j=1,nres_B
                  dist=sqrt((res_xsc_A0(i)-res_x_B0(j))**2+
     &                 (res_ysc_A0(i)-res_y_B0(j))**2+
     &                 (res_zsc_A0(i)-res_z_B0(j))**2)
                  if(dist.le.initial_clash)then
                     clash_flag=1
                  endif
               enddo
            enddo
            

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c>>   2) with certain distance and angular constrain

c>>>>>>>>>>>>>   calculate distance of two reaction sites
            dist=sqrt((RB_B_x_01(4)-
     &           RB_A_x_01(4))**2+
     &           (RB_B_y_01(4)-
     &           RB_A_y_01(4))**2+
     &           (RB_B_z_01(4)-
     &           RB_A_z_01(4))**2)
c>>>>>>>>>>>> calculate two angle
c>> perpendicular angle
            theta_pd=0
            point_x(1)=RB_A_x_01(1)-RB_A_x_01(4)
            point_y(1)=RB_A_y_01(1)-RB_A_y_01(4)
            point_z(1)=RB_A_z_01(1)-RB_A_z_01(4)
            point_x(2)=0
            point_y(2)=0
            point_z(2)=0
            point_x(3)=RB_B_x_01(1)-RB_B_x_01(4)
            point_y(3)=RB_B_y_01(1)-RB_B_y_01(4)
            point_z(3)=RB_B_z_01(1)-RB_B_z_01(4)
            call gettheta(point_x,point_y,point_z,
     &           theta_pd)
c>> orientational angle
            theta_ot=0
            point_x(1)=RB_A_x_01(1)-RB_A_x_01(2)
            point_y(1)=RB_A_y_01(1)-RB_A_y_01(2)
            point_z(1)=RB_A_z_01(1)-RB_A_z_01(2)
            point_x(2)=0
            point_y(2)=0
            point_z(2)=0
            point_x(3)=RB_B_x_01(1)-RB_B_x_01(2)
            point_y(3)=RB_B_y_01(1)-RB_B_y_01(2)
            point_z(3)=RB_B_z_01(1)-RB_B_z_01(2)
            call gettheta(point_x,point_y,point_z,
     &           theta_ot)
            if((dist.gt.dist_constrain).OR.
     &           (clash_flag.eq.1))then
               goto 100
            endif

c>>  calculate how many native contact formed

            native_contact=0
            do i=1,nres_A
               do j=1,nres_B
                  dist=sqrt((res_x_A0(i)-res_x_B0(j))**2+
     &                 (res_y_A0(i)-res_y_B0(j))**2+
     &                 (res_z_A0(i)-res_z_B0(j))**2)
                  pair_flag=0
                  pair_index=0
                  do k=1,pair_num
                     if((pair_list(k,1).eq.i).AND.
     &                    (pair_list(k,2).eq.j))then
                        pair_flag=1
                        pair_index=k
                     endif
                  enddo                  
                  if(pair_flag.eq.1)then                     
                     if(dist.lt.(respair_dist(pair_index)
     &                    +well_width))then
                        native_contact=native_contact+1          
                     endif                      
                  endif                  
               enddo
            enddo
            native_contact_ratio=real(native_contact)/real(respair_num)      

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  end of construction of initial conformation
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            
 2100       format(A6,I5,1x,A4,1x,A3,1x,A1,I4,4x,3F8.3)
 3000       format(A3)     
            
            do j=1,nres_A
               res_xsc_A_old(j)=res_xsc_A0(j)
               res_ysc_A_old(j)=res_ysc_A0(j)
               res_zsc_A_old(j)=res_zsc_A0(j)
            enddo
            do j=1,nres_B
               res_xsc_B_old(j)=res_xsc_B0(j)
               res_ysc_B_old(j)=res_ysc_B0(j)
               res_zsc_B_old(j)=res_zsc_B0(j)
            enddo
            
            do j=1,nres_A
               res_x_A_old(j)=res_x_A0(j)
               res_y_A_old(j)=res_y_A0(j)
               res_z_A_old(j)=res_z_A0(j)
            enddo
            do j=1,nres_B
               res_x_B_old(j)=res_x_B0(j)
               res_y_B_old(j)=res_y_B0(j)
               res_z_B_old(j)=res_z_B0(j)
            enddo
            
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc  calculate the total energy of current conformation
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c>>  long-range energy
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c>>  repulsive iteraction as Van Der Waals potential

            ene_vdw=0
            
            do i=1,nres_A
               do j=1,nres_B                   
                  rij2=sqrt((res_x_A_old(i)-res_x_B_old(j))**2+
     &                 (res_y_A_old(i)-res_y_B_old(j))**2+
     &                 (res_z_A_old(i)-res_z_B_old(j))**2)
                  if(rij2.lt.3.0)then
                     inverser=1/rij2
                     inverser6=inverser**6
                     ene_vdw=ene_vdw+
     &                    2*125*125*
     &                    ((125*125/2)*(inverser6**2)-inverser6)
                  endif
               enddo
            enddo
            
            do i=1,nres_A
               do j=1,nres_B                   
                  rij2=sqrt((res_xsc_A_old(i)-res_xsc_B_old(j))**2+
     &                 (res_ysc_A_old(i)-res_ysc_B_old(j))**2+
     &                 (res_zsc_A_old(i)-res_zsc_B_old(j))**2)
                  if(rij2.lt.3.0)then
                     inverser=1/rij2
                     inverser6=inverser**6
                     ene_vdw=ene_vdw+
     &                    2*125*125*
     &                    ((125*125/2)*(inverser6**2)-inverser6)
                  endif
               enddo
            enddo
            
            do i=1,nres_A
               do j=1,nres_B                   
                  rij2=sqrt((res_x_A_old(i)-res_xsc_B_old(j))**2+
     &                 (res_y_A_old(i)-res_ysc_B_old(j))**2+
     &                 (res_z_A_old(i)-res_zsc_B_old(j))**2)
                  if(rij2.lt.3.0)then
                     inverser=1/rij2
                     inverser6=inverser**6
                     ene_vdw=ene_vdw+
     &                    2*125*125*
     &                    ((125*125/2)*(inverser6**2)-inverser6)
                  endif
               enddo
            enddo
            
            do i=1,nres_A
               do j=1,nres_B                   
                  rij2=sqrt((res_xsc_A_old(i)-res_x_B_old(j))**2+
     &                 (res_ysc_A_old(i)-res_y_B_old(j))**2+
     &                 (res_zsc_A_old(i)-res_z_B_old(j))**2)
                  if(rij2.lt.3.0)then
                     inverser=1/rij2
                     inverser6=inverser**6
                     ene_vdw=ene_vdw+
     &                    2*125*125*
     &                    ((125*125/2)*(inverser6**2)-inverser6)
                  endif
               enddo
            enddo
            


c>> the Solvant effect

            SolvantEffect=0

            do i=1,nres_A
               do j=1,nres_B
                  
                  dist=sqrt((res_xsc_A_old(i)-res_xsc_B_old(j))**2+
     &                 (res_ysc_A_old(i)-res_ysc_B_old(j))**2+
     &                 (res_zsc_A_old(i)-res_zsc_B_old(j))**2)
                  
                  if(dist.le.7.5)then
                     
                     if(res_na_A(i).eq.'ALA')then
                        SolvantEffect=SolvantEffect-1.8
                     elseif(res_na_A(i).eq.'ILE')then
                        SolvantEffect=SolvantEffect-4.5
                     elseif(res_na_A(i).eq.'LEU')then
                        SolvantEffect=SolvantEffect-3.8
                     elseif(res_na_A(i).eq.'PHE')then
                        SolvantEffect=SolvantEffect-2.8
                     elseif(res_na_A(i).eq.'VAL')then
                        SolvantEffect=SolvantEffect-4.2
                     elseif(res_na_A(i).eq.'MET')then
                        SolvantEffect=SolvantEffect-1.9
                     elseif(res_na_A(i).eq.'CYS')then
                        SolvantEffect=SolvantEffect-2.5
                     endif
                     
                     if(res_na_B(j).eq.'ALA')then
                        SolvantEffect=SolvantEffect-1.8
                     elseif(res_na_B(j).eq.'ILE')then
                        SolvantEffect=SolvantEffect-4.5
                     elseif(res_na_B(j).eq.'LEU')then
                        SolvantEffect=SolvantEffect-3.8
                     elseif(res_na_B(j).eq.'PHE')then
                        SolvantEffect=SolvantEffect-2.8
                     elseif(res_na_B(j).eq.'VAL')then
                        SolvantEffect=SolvantEffect-4.2
                     elseif(res_na_B(j).eq.'MET')then
                        SolvantEffect=SolvantEffect-1.9
                     elseif(res_na_B(j).eq.'CYS')then
                        SolvantEffect=SolvantEffect-2.5
                     endif
                     
                  endif
               enddo
            enddo

c>>   the electrostatic potential

            ElectricPotential=0
            charge_A=0
            charge_B=0
            do i=1,nres_A
               if((res_na_A(i).eq.'ARG')
     &              .OR.(res_na_A(i).eq.'LYS'))then
                  charge_A=1.0
               elseif((res_na_A(i).eq.'GLU')
     &                 .OR.(res_na_A(i).eq.'ASP'))then
                  charge_A=-1.0
               elseif(res_na_A(i).eq.'HIS')then
                  charge_A=0.5
               else
                  charge_A=0.0
               endif
               do j=1,nres_B
                  if((res_na_B(j).eq.'ARG')
     &                 .OR.(res_na_B(j).eq.'LYS'))then
                     charge_B=1.0
                  elseif((res_na_B(j).eq.'GLU')
     &                    .OR.(res_na_B(j).eq.'ASP'))then
                     charge_B=-1.0
                  elseif(res_na_B(j).eq.'HIS')then
                     charge_B=0.5
                  else
                     charge_B=0.0
                  endif
    
                  
                  dist=sqrt((res_xsc_A_old(i)-res_xsc_B_old(j))**2+
     &                 (res_ysc_A_old(i)-res_ysc_B_old(j))**2+
     &                 (res_zsc_A_old(i)-res_zsc_B_old(j))**2)
                  ElectricPotential=ElectricPotential+
     &                 ElectricConstant*charge_A*charge_B/
     &                 (10.0*dist*DEXP(dist/DebyeLength))


               enddo
            enddo

            sum_diff=ene_vdw+SolvantEffect*W_S+ElectricPotential*W_E
            
            
            formDimer_flag=0
            formDimer_time=simu_step*time_step
            record_formDimer_flag(n_t)=0
            record_formDimer_time(n_t)=simu_step*time_step

cccccccccccccccccccccccccccccccccccccccccccccc

            do mc_time_step=1,simu_step
            
               simu_temp=temp_ini-
     &              (temp_ini-temp_fin)/real(simu_step)
     &              *real(mc_time_step)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc        begin main simulation of association in current traj
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c>>>    diffusion of molecule A:
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc  global translational diffusion as RB
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            
               distance_step=2*sqrt(diffusion_constant_A
     &              *time_step/6.0)*rand3(r3)
               theta=rand4(r4)*pai
               phai=rand3(r3)*2*pai
               do j=1,nres_A
                  res_x_A_lcts(j)=
     &                 res_x_A_old(j)+
     &                 distance_step*sin(theta)*cos(phai)
                  res_y_A_lcts(j)=
     &                 res_y_A_old(j)+
     &                 distance_step*sin(theta)*sin(phai)
                  res_z_A_lcts(j)=
     &                 res_z_A_old(j)+
     &                 distance_step*cos(theta)
               enddo
               do j=1,nres_A
                  res_xsc_A_lcts(j)=
     &                 res_xsc_A_old(j)+
     &                 distance_step*sin(theta)*cos(phai)
                  res_ysc_A_lcts(j)=
     &                 res_ysc_A_old(j)+
     &                 distance_step*sin(theta)*sin(phai)
                  res_zsc_A_lcts(j)=
     &                 res_zsc_A_old(j)+
     &                 distance_step*cos(theta)
               enddo
               
     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc   rotational diffusion as RB
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               cm1_a_x=0
               cm1_a_y=0
               cm1_a_z=0
               do j=1,nres_A
                  cm1_a_x=cm1_a_x+res_x_A_lcts(j)
                  cm1_a_y=cm1_a_y+res_y_A_lcts(j)
                  cm1_a_z=cm1_a_z+res_z_A_lcts(j)
               enddo
               cm1_a_x=cm1_a_x/(real(nres_A))
               cm1_a_y=cm1_a_y/(real(nres_A))
               cm1_a_z=cm1_a_z/(real(nres_A))
               
               theta=(2*rand3(r3)-1)*2*rot_D_A
     &              *time_step*pai/180.0
               phi=(2*rand3(r3)-1)*2*rot_D_A
     &              *time_step*pai/180.0 
               psai=(2*rand4(r4)-1)*2*rot_D_A
     &              *time_step*pai/180.0
               
               t(1,1)=cos(psai)*cos(phi)-
     &              cos(theta)*sin(phi)*sin(psai)
               t(1,2)=-sin(psai)*cos(phi)-
     &              cos(theta)*sin(phi)*cos(psai)
               t(1,3)=sin(theta)*sin(phi)
               
               t(2,1)=cos(psai)*sin(phi)+
     &              cos(theta)*cos(phi)*sin(psai)
               t(2,2)=-sin(psai)*sin(phi)+
     &              cos(theta)*cos(phi)*cos(psai)
               t(2,3)=-sin(theta)*cos(phi)
               
               t(3,1)=sin(psai)*sin(theta)
               t(3,2)=cos(psai)*sin(theta)
               t(3,3)=cos(theta)
               
               do k=1,nres_A
                  res_x_A_new(k)=
     &                 t(1,1)*(res_x_A_lcts(k)-cm1_a_x)+
     &                 t(1,2)*(res_y_A_lcts(k)-cm1_a_y)+
     &                 t(1,3)*(res_z_A_lcts(k)-cm1_a_z)+
     &                 cm1_a_x
                  res_y_A_new(k)=
     &                 t(2,1)*(res_x_A_lcts(k)-cm1_a_x)+
     &                 t(2,2)*(res_y_A_lcts(k)-cm1_a_y)+
     &                 t(2,3)*(res_z_A_lcts(k)-cm1_a_z)+
     &                 cm1_a_y
                  res_z_A_new(k)=
     &                 t(3,1)*(res_x_A_lcts(k)-cm1_a_x)+
     &                 t(3,2)*(res_y_A_lcts(k)-cm1_a_y)+
     &                 t(3,3)*(res_z_A_lcts(k)-cm1_a_z)+
     &                 cm1_a_z
               enddo
               do k=1,nres_A
                  res_xsc_A_new(k)=
     &                 t(1,1)*(res_xsc_A_lcts(k)-cm1_a_x)+
     &                 t(1,2)*(res_ysc_A_lcts(k)-cm1_a_y)+
     &                 t(1,3)*(res_zsc_A_lcts(k)-cm1_a_z)+
     &                 cm1_a_x
                  res_ysc_A_new(k)=
     &                 t(2,1)*(res_xsc_A_lcts(k)-cm1_a_x)+
     &                 t(2,2)*(res_ysc_A_lcts(k)-cm1_a_y)+
     &                 t(2,3)*(res_zsc_A_lcts(k)-cm1_a_z)+
     &                 cm1_a_y
                  res_zsc_A_new(k)=
     &                 t(3,1)*(res_xsc_A_lcts(k)-cm1_a_x)+
     &                 t(3,2)*(res_ysc_A_lcts(k)-cm1_a_y)+
     &                 t(3,3)*(res_zsc_A_lcts(k)-cm1_a_z)+
     &                 cm1_a_z
               enddo
            
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c>>>    diffusion of molecule B:
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc  global translational diffusion as RB
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            
               distance_step=2*sqrt(diffusion_constant_B
     &              *time_step/6.0)*rand3(r3)
               theta=rand3(r3)*pai
               phai=rand3(r3)*2*pai
               do j=1,nres_B
                  res_x_B_lcts(j)=
     &                 res_x_B_old(j)+
     &                 distance_step*sin(theta)*cos(phai)
                  res_y_B_lcts(j)=
     &                 res_y_B_old(j)+
     &                 distance_step*sin(theta)*sin(phai)
                  res_z_B_lcts(j)=
     &                 res_z_B_old(j)+
     &                 distance_step*cos(theta)
               enddo
               do j=1,nres_B
                  res_xsc_B_lcts(j)=
     &                 res_xsc_B_old(j)+
     &                 distance_step*sin(theta)*cos(phai)
                  res_ysc_B_lcts(j)=
     &                 res_ysc_B_old(j)+
     &                 distance_step*sin(theta)*sin(phai)
                  res_zsc_B_lcts(j)=
     &                 res_zsc_B_old(j)+
     &                 distance_step*cos(theta)
               enddo


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc   rotational diffusion as RB
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               cm1_b_x=0
               cm1_b_y=0
               cm1_b_z=0
               do j=1,nres_B
                  cm1_b_x=cm1_b_x+res_x_B_lcts(j)
                  cm1_b_y=cm1_b_y+res_y_B_lcts(j)
                  cm1_b_z=cm1_b_z+res_z_B_lcts(j)
               enddo
               cm1_b_x=cm1_b_x/(real(nres_B))
               cm1_b_y=cm1_b_y/(real(nres_B))
               cm1_b_z=cm1_b_z/(real(nres_B))
               
               theta=(2*rand3(r3)-1)*2*rot_D_B
     &              *time_step*pai/180.0
               phi=(2*rand4(r4)-1)*2*rot_D_B
     &              *time_step*pai/180.0 
               psai=(2*rand3(r3)-1)*2*rot_D_B
     &              *time_step*pai/180.0
               
               t(1,1)=cos(psai)*cos(phi)-
     &              cos(theta)*sin(phi)*sin(psai)
               t(1,2)=-sin(psai)*cos(phi)-
     &              cos(theta)*sin(phi)*cos(psai)
               t(1,3)=sin(theta)*sin(phi)
               
               t(2,1)=cos(psai)*sin(phi)+
     &              cos(theta)*cos(phi)*sin(psai)
               t(2,2)=-sin(psai)*sin(phi)+
     &              cos(theta)*cos(phi)*cos(psai)
               t(2,3)=-sin(theta)*cos(phi)
               
               t(3,1)=sin(psai)*sin(theta)
               t(3,2)=cos(psai)*sin(theta)
               t(3,3)=cos(theta)
               
               do k=1,nres_B
                  res_x_B_new(k)=
     &                 t(1,1)*(res_x_B_lcts(k)-cm1_b_x)+
     &                 t(1,2)*(res_y_B_lcts(k)-cm1_b_y)+
     &                 t(1,3)*(res_z_B_lcts(k)-cm1_b_z)+
     &                 cm1_b_x
                  res_y_B_new(k)=
     &                 t(2,1)*(res_x_B_lcts(k)-cm1_b_x)+
     &                 t(2,2)*(res_y_B_lcts(k)-cm1_b_y)+
     &                 t(2,3)*(res_z_B_lcts(k)-cm1_b_z)+
     &                 cm1_b_y
                  res_z_B_new(k)=
     &                 t(3,1)*(res_x_B_lcts(k)-cm1_b_x)+
     &                 t(3,2)*(res_y_B_lcts(k)-cm1_b_y)+
     &                 t(3,3)*(res_z_B_lcts(k)-cm1_b_z)+
     &                 cm1_b_z
               enddo
               
               do k=1,nres_B
                  res_xsc_B_new(k)=
     &                 t(1,1)*(res_xsc_B_lcts(k)-cm1_b_x)+
     &                 t(1,2)*(res_ysc_B_lcts(k)-cm1_b_y)+
     &                 t(1,3)*(res_zsc_B_lcts(k)-cm1_b_z)+
     &                 cm1_b_x
                  res_ysc_B_new(k)=
     &                 t(2,1)*(res_xsc_B_lcts(k)-cm1_b_x)+
     &                 t(2,2)*(res_ysc_B_lcts(k)-cm1_b_y)+
     &                 t(2,3)*(res_zsc_B_lcts(k)-cm1_b_z)+
     &                 cm1_b_y
                  res_zsc_B_new(k)=
     &                 t(3,1)*(res_xsc_B_lcts(k)-cm1_b_x)+
     &                 t(3,2)*(res_ysc_B_lcts(k)-cm1_b_y)+
     &                 t(3,3)*(res_zsc_B_lcts(k)-cm1_b_z)+
     &                 cm1_b_z
               enddo
               
c            print*,n_t,mc_time_step,cm1_a_x,cm1_a_y,cm1_a_z,
c     &           cm1_b_x,cm1_b_y,cm1_b_z
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc  calculate the total energy of current conformation
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c>>  long-range energy
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c>>  repulsive iteraction as Van Der Waals potential

               ene_vdw=0
             
               do i=1,nres_A
                  do j=1,nres_B                   
                     rij2=sqrt((res_x_A_new(i)-res_x_B_new(j))**2+
     &                    (res_y_A_new(i)-res_y_B_new(j))**2+
     &                    (res_z_A_new(i)-res_z_B_new(j))**2)
                     if(rij2.lt.3.0)then
                        inverser=1/rij2
                        inverser6=inverser**6
                        ene_vdw=ene_vdw+
     &                       2*125*125*
     &                       ((125*125/2)*(inverser6**2)-inverser6)
                     endif
                  enddo
               enddo
               
               do i=1,nres_A
                  do j=1,nres_B                   
                     rij2=sqrt((res_xsc_A_new(i)-res_xsc_B_new(j))**2+
     &                    (res_ysc_A_new(i)-res_ysc_B_new(j))**2+
     &                    (res_zsc_A_new(i)-res_zsc_B_new(j))**2)
                     if(rij2.lt.3.0)then
                        inverser=1/rij2
                        inverser6=inverser**6
                        ene_vdw=ene_vdw+
     &                       2*125*125*
     &                       ((125*125/2)*(inverser6**2)-inverser6)
                     endif
                  enddo
               enddo
               
               do i=1,nres_A
                  do j=1,nres_B                   
                     rij2=sqrt((res_x_A_new(i)-res_xsc_B_new(j))**2+
     &                    (res_y_A_new(i)-res_ysc_B_new(j))**2+
     &                    (res_z_A_new(i)-res_zsc_B_new(j))**2)
                     if(rij2.lt.3.0)then
                        inverser=1/rij2
                        inverser6=inverser**6
                        ene_vdw=ene_vdw+
     &                       2*125*125*
     &                       ((125*125/2)*(inverser6**2)-inverser6)
                     endif
                  enddo
               enddo
               
               do i=1,nres_A
                  do j=1,nres_B                   
                     rij2=sqrt((res_xsc_A_new(i)-res_x_B_new(j))**2+
     &                    (res_ysc_A_new(i)-res_y_B_new(j))**2+
     &                    (res_zsc_A_new(i)-res_z_B_new(j))**2)
                     if(rij2.lt.3.0)then
                        inverser=1/rij2
                        inverser6=inverser**6
                        ene_vdw=ene_vdw+
     &                       2*125*125*
     &                       ((125*125/2)*(inverser6**2)-inverser6)
                     endif
                  enddo
               enddo
              

c>> the Solvant effect


               SolvantEffect=0
               
               do i=1,nres_A
                  do j=1,nres_B
                     
                     dist=sqrt((res_xsc_A_new(i)-res_xsc_B_new(j))**2+
     &                    (res_ysc_A_new(i)-res_ysc_B_new(j))**2+
     &                    (res_zsc_A_new(i)-res_zsc_B_new(j))**2)
                     
                     if(dist.le.7.5)then
                        
                        if(res_na_A(i).eq.'ALA')then
                           SolvantEffect=SolvantEffect-1.8
                        elseif(res_na_A(i).eq.'ILE')then
                           SolvantEffect=SolvantEffect-4.5
                        elseif(res_na_A(i).eq.'LEU')then
                           SolvantEffect=SolvantEffect-3.8
                        elseif(res_na_A(i).eq.'PHE')then
                           SolvantEffect=SolvantEffect-2.8
                        elseif(res_na_A(i).eq.'VAL')then
                           SolvantEffect=SolvantEffect-4.2
                        elseif(res_na_A(i).eq.'MET')then
                           SolvantEffect=SolvantEffect-1.9
                        elseif(res_na_A(i).eq.'CYS')then
                           SolvantEffect=SolvantEffect-2.5
                        endif
                        
                        if(res_na_B(j).eq.'ALA')then
                           SolvantEffect=SolvantEffect-1.8
                        elseif(res_na_B(j).eq.'ILE')then
                           SolvantEffect=SolvantEffect-4.5
                        elseif(res_na_B(j).eq.'LEU')then
                           SolvantEffect=SolvantEffect-3.8
                        elseif(res_na_B(j).eq.'PHE')then
                           SolvantEffect=SolvantEffect-2.8
                        elseif(res_na_B(j).eq.'VAL')then
                           SolvantEffect=SolvantEffect-4.2
                        elseif(res_na_B(j).eq.'MET')then
                           SolvantEffect=SolvantEffect-1.9
                        elseif(res_na_B(j).eq.'CYS')then
                           SolvantEffect=SolvantEffect-2.5
                        endif
                        
                     endif
c              
                  enddo
               enddo

c>>   the electrostatic potential

               ElectricPotential=0
               charge_A=0
               charge_B=0
               do i=1,nres_A
                  if((res_na_A(i).eq.'ARG')
     &                 .OR.(res_na_A(i).eq.'LYS'))then
                     charge_A=1.0
                  elseif((res_na_A(i).eq.'GLU')
     &                    .OR.(res_na_A(i).eq.'ASP'))then
                     charge_A=-1.0
                  elseif(res_na_A(i).eq.'HIS')then
                     charge_A=0.5
                  else
                     charge_A=0.0
                  endif
                  do j=1,nres_B
                     if((res_na_B(j).eq.'ARG')
     &                    .OR.(res_na_B(j).eq.'LYS'))then
                        charge_B=1.0
                     elseif((res_na_B(j).eq.'GLU')
     &                       .OR.(res_na_B(j).eq.'ASP'))then
                        charge_B=-1.0
                     elseif(res_na_B(j).eq.'HIS')then
                        charge_B=0.5
                     else
                        charge_B=0.0
                     endif
                     
                     
                     dist=sqrt((res_xsc_A_new(i)-res_xsc_B_new(j))**2+
     &                    (res_ysc_A_new(i)-res_ysc_B_new(j))**2+
     &                    (res_zsc_A_new(i)-res_zsc_B_new(j))**2)
                     ElectricPotential=ElectricPotential+
     &                    ElectricConstant*charge_A*charge_B/
     &                    (10.0*dist*DEXP(dist/DebyeLength))
                     
c             
                  enddo
               enddo
            
               sum_diff2=ene_vdw+SolvantEffect*W_S+ElectricPotential*W_E


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Monte-carlo metropolis algorithm
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               if((sum_diff2.le.sum_diff).OR.
     &              ((sum_diff2.gt.sum_diff).AND.
     &              (exp(-(sum_diff2-sum_diff)/simu_temp)
     &              .gt.rand5(r5))))then
                  sum_diff=sum_diff2
                  
                  do i=1,nres_A
                     res_x_A_old(i)=res_x_A_new(i)
                     res_y_A_old(i)=res_y_A_new(i)
                     res_z_A_old(i)=res_z_A_new(i)
                  enddo
                  do i=1,nres_B
                     res_x_B_old(i)=res_x_B_new(i)
                     res_y_B_old(i)=res_y_B_new(i)
                     res_z_B_old(i)=res_z_B_new(i)
                  enddo
                  
                  do i=1,nres_A
                     res_xsc_A_old(i)=res_xsc_A_new(i)
                     res_ysc_A_old(i)=res_ysc_A_new(i)
                     res_zsc_A_old(i)=res_zsc_A_new(i)
                  enddo
                  do i=1,nres_B
                     res_xsc_B_old(i)=res_xsc_B_new(i)
                     res_ysc_B_old(i)=res_ysc_B_new(i)
                     res_zsc_B_old(i)=res_zsc_B_new(i)
                  enddo
                  
               endif   

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc  calculate how many native contact formed
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               native_contact=0
               do i=1,nres_A
                  do j=1,nres_B
                     dist=sqrt((res_x_A_old(i)-res_x_B_old(j))**2+
     &                    (res_y_A_old(i)-res_y_B_old(j))**2+
     &                    (res_z_A_old(i)-res_z_B_old(j))**2)
                     pair_flag=0
                     pair_index=0
                     do k=1,respair_num
                        if((respair_list(k,1).eq.res_seq_A(i)).AND.
     &                       (respair_list(k,2).eq.res_seq_B(j)))then
                           pair_flag=1
                           pair_index=k
                        endif
                     enddo                  
                     if(pair_flag.eq.1)then                     
                        if(dist.lt.(respair_dist(pair_index)
     &                       +well_width))then
                           native_contact=native_contact+1          
                        endif                     
                     endif                  
                  enddo
               enddo
               native_contact_ratio=real(native_contact)
     &              /real(respair_num)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc  calculate the RMSD between A and B
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               do j=1,1000
                  xa(j)=0
                  ya(j)=0
                  za(j)=0
                  xb(j)=0
                  yb(j)=0
                  zb(j)=0
                  atompr(1,j)=0
                  atompr(2,j)=0
               enddo
               index=0
               do i=1,nres_A
                  index=index+1
                  xa(index)=res_x_A(i)
                  ya(index)=res_y_A(i)
                  za(index)=res_z_A(i)
               enddo
               do i=1,nres_B
                  index=index+1
                  xa(index)=res_x_B(i)
                  ya(index)=res_y_B(i)
                  za(index)=res_z_B(i)
               enddo
               index=0
               do i=1,nres_A
                  index=index+1
                  xb(index)=res_x_A_old(i)
                  yb(index)=res_y_A_old(i)
                  zb(index)=res_z_A_old(i)
               enddo
               do i=1,nres_B
                  index=index+1
                  xb(index)=res_x_B_old(i)
                  yb(index)=res_y_B_old(i)
                  zb(index)=res_z_B_old(i)
               enddo
               do i=1,nres_A+nres_B
                  atompr(1,i)=i
                  atompr(2,i)=i
               enddo
               npair=nres_A+nres_B       
               npair_a=nres_A+nres_B
               npair_b=nres_A+nres_B
               call ROTLSQ
     &              (xa,ya,za,npair_a,xb,yb,zb,npair_b,atompr,npair)
               RMSD=0
               do i=1,nres_A+nres_B
                  RMSD=RMSD+(xa(i)-xb(i))**2
     &                 +(ya(i)-yb(i))**2
     &                 +(za(i)-zb(i))**2
               enddo
               RMSD=sqrt(RMSD/(real(nres_A+nres_B)))
            

ccccccccccccccccccccccccccccccccccccccccccccc
cc  distance between two proteins
cccccccccccccccccccccccccccccccccccccccccccc

               cm1_a_x=0
               cm1_a_y=0
               cm1_a_z=0
               do j=1,nres_A
                  cm1_a_x=cm1_a_x+res_x_A_old(j)
                  cm1_a_y=cm1_a_y+res_y_A_old(j)
                  cm1_a_z=cm1_a_z+res_z_A_old(j)
               enddo
               cm1_a_x=cm1_a_x/(real(nres_A))
               cm1_a_y=cm1_a_y/(real(nres_A))
               cm1_a_z=cm1_a_z/(real(nres_A))
               
               cm1_b_x=0
               cm1_b_y=0
               cm1_b_z=0
               do j=1,nres_B
                  cm1_b_x=cm1_b_x+res_x_B_old(j)
                  cm1_b_y=cm1_b_y+res_y_B_old(j)
                  cm1_b_z=cm1_b_z+res_z_B_old(j)
               enddo
               cm1_b_x=cm1_b_x/(real(nres_B))
               cm1_b_y=cm1_b_y/(real(nres_B))
               cm1_b_z=cm1_b_z/(real(nres_B))
               
               dist_AB_COM=sqrt((cm1_a_x-cm1_b_x)**2
     &              +(cm1_a_y-cm1_b_y)**2
     &              +(cm1_a_z-cm1_b_z)**2)
               
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc  output
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c            print*,n_t,mc_time_step,sum_diff,native_contact_ratio,RMSD

               if((output_flag_str.eq.1).AND.
     &              (MOD(mc_time_step,trj_output_freq).eq.0))then            
                  open (unit=10,file=
     &                 'Result_ResNP_trj_'//pdbid(id)
     &                 //'_'//
     &                 chainid_A2(1)//chainid_A2(2)//chainid_A2(3)
     &                 //'_'//
     &                 chainid_B2(1)//chainid_B2(2)//chainid_B2(3)
     &                 //'_d'//
     &                 char(INT(dist_constrain_matr(id)/10.0)+48)//
     &                 char(int(dist_constrain_matr(id)
     &                 -int(dist_constrain_matr(id)/10.0)*10.0)+48)
     &                 //'.pdb',
     &                 status='unknown',access='append')
                  write(10,5000) 'index',n_t,mc_time_step,sum_diff
                  write(10,2100)('ATOM  ',j,' CA ',res_na_A(j), 
     &                 'A',j,res_x_A_old(j),res_y_A_old(j),
     &                 res_z_A_old(j),j=1,nres_A)
                  write(10,3000) 'TER'
                  write(10,2100)('ATOM  ',j,' CA ',res_na_B(j), 
     &                 'B',j,res_x_B_old(j),res_y_B_old(j),
     &                 res_z_B_old(j),j=1,nres_B)
                  write(10,3000) 'TER'
                  write(10,3000) 'END'      
                  close(10)         
               endif
               
               if((output_flag_trjlog.eq.1).AND.
     &              (MOD(mc_time_step,trj_output_freq).eq.0))then       
                  open (unit=10,file=
     &                 'Result_ResNP_trjlog_'//pdbid(id)
     &                 //'_'//
     &                 chainid_A2(1)//chainid_A2(2)//chainid_A2(3)
     &                 //'_'//
     &                 chainid_B2(1)//chainid_B2(2)//chainid_B2(3)
     &                 //'_d'//
     &                 char(INT(dist_constrain_matr(id)/10.0)+48)//
     &                 char(int(dist_constrain_matr(id)
     &                 -int(dist_constrain_matr(id)/10.0)*10.0)+48)
     &                 //'.dat',
     &                 status='unknown',access='append')
                  write(10,1979) 'index',n_t,mc_time_step,
     &                 dist_AB_COM,native_contact_ratio,RMSD
                  close(10)
               endif
 1979          format(A5,I3,I10,F12.5,1x,F12.5,1x,F12.5)
               
 5000          format(A5,I3,I10,F12.5)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    judge if native dimer formed
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               if(native_contact.ge.3)then

                  formDimer_flag=1
                  formDimer_time=mc_time_step*time_step
                  do k=1,nres_A
                     res_x_A_fin(k,n_t)=res_x_A_old(k)
                     res_y_A_fin(k,n_t)=res_y_A_old(k)
                     res_z_A_fin(k,n_t)=res_z_A_old(k)
                  enddo
                  do k=1,nres_B
                     res_x_B_fin(k,n_t)=res_x_B_old(k)
                     res_y_B_fin(k,n_t)=res_y_B_old(k)
                     res_z_B_fin(k,n_t)=res_z_B_old(k)
                  enddo
                  goto 1980

               endif
               
               if((mc_time_step.eq.simu_step)
     &              .AND.
     &              (native_contact.lt.3)
     &              )then
                  formDimer_flag=0
c               formDimer_time=mc_time_step*time_step
                  do k=1,nres_A
                     res_x_A_fin(k,n_t)=res_x_A_old(k)
                     res_y_A_fin(k,n_t)=res_y_A_old(k)
                     res_z_A_fin(k,n_t)=res_z_A_old(k)
                  enddo
                  do k=1,nres_B
                     res_x_B_fin(k,n_t)=res_x_B_old(k)
                     res_y_B_fin(k,n_t)=res_y_B_old(k)
                     res_z_B_fin(k,n_t)=res_z_B_old(k)
                  enddo
               endif

cccccccccccccccccccccccccccccccccccccccccccc
cc   end main simulation of accociation
cccccccccccccccccccccccccccccccccccccccccccc
 1977          continue
            enddo
            
 1980       continue
            
            if(formDimer_flag.eq.1)then
               record_formDimer_flag(n_t)=1
               record_formDimer_time(n_t)=formDimer_time
            endif
            
            if(output_flag_rec.eq.1)then
               open (unit=10,file=
     &              'output/Result_ResNP_rec_'//pdbid(id)
     &              //'_'//
     &              chainid_A2(1)//chainid_A2(2)//chainid_A2(3)
     &              //'_'//
     &              chainid_B2(1)//chainid_B2(2)//chainid_B2(3)
     &              //'_d'//
     &              char(INT(dist_constrain_matr(id)/10.0)+48)//
     &              char(int(dist_constrain_matr(id)
     &              -int(dist_constrain_matr(id)/10.0)*10.0)+48)
     &              //'.dat',
     &              status='unknown',access='append')
               write(10,1978) n_t,record_formDimer_flag(n_t),
     &              record_formDimer_time(n_t),native_contact_ratio,RMSD
               close(10)
            endif
            if(output_flag_trjrec.eq.1)then
               open (unit=10,file=
     &              'output/Result_ResNP_recstr_'//pdbid(id)
     &              //'_'//
     &              chainid_A2(1)//chainid_A2(2)//chainid_A2(3)
     &              //'_'//
     &              chainid_B2(1)//chainid_B2(2)//chainid_B2(3)
     &              //'_d'//
     &              char(INT(dist_constrain_matr(id)/10.0)+48)//
     &              char(int(dist_constrain_matr(id)
     &              -int(dist_constrain_matr(id)/10.0)*10.0)+48)
     &              //'.pdb',
     &              status='unknown',access='append')
               write(10,5000) 'index',n_t,mc_time_step,sum_diff
               write(10,2100)('ATOM  ',j,' CA ',res_na_A(j), 
     &              'A',j,res_x_A_fin(j,n_t),res_y_A_fin(j,n_t),
     &              res_z_A_fin(j,n_t),j=1,nres_A)
               write(10,3000) 'TER'
               write(10,2100)('ATOM  ',j,' CA ',res_na_B(j), 
     &              'B',j,res_x_B_fin(j,n_t),res_y_B_fin(j,n_t),
     &              res_z_B_fin(j,n_t),j=1,nres_B)
               write(10,3000) 'TER'
               write(10,3000) 'END'      
               close(10)         
            endif
 1978       format(I5,I3,1x,3F15.5)


ccccccccccccccccccccccccccccccccccc
c   end of multiple traj
ccccccccccccccccccccccccccccccccccc
            
         enddo

cccccccccccccccccccccccccccc
c  end of multiple pdb
cccccccccccccccccccccccccccc

      enddo
      
      stop
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc

      real  function rand3(r3)
      double precision s,u,v,r3
      s=65536.0
      u=2053.0
      v=13849.0
      m=r3/s
      r3=r3-m*s
      r3=u*r3+v
      m=r3/s
      r3=r3-m*s
      rand3=r3/s
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccc

      real  function rand4(r4)
      double precision s,u,v,r4
      s=65536.0
      u=2053.0
      v=13849.0
      m=r4/s
      r4=r4-m*s
      r4=u*r4+v
      m=r4/s
      r4=r4-m*s
      rand4=r4/s
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc                      
ccccccccccccccccccccccccccccccccccccccccccccccccc

      real  function rand5(r5)
      double precision s,u,v,r5
      s=65536.0
      u=2053.0
      v=13849.0
      m=r5/s
      r5=r5-m*s
      r5=u*r5+v
      m=r5/s
      r5=r5-m*s
      rand5=r5/s
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc    
ccccccccccccccccccccccccccccccccccccccccccccccccc            
            
      subroutine gettheta(point_x,point_y,point_z,test_theta)

      implicit none
      real*8 point_x(3),point_y(3),point_z(3),test_theta
ccc
      real*8 lx(2),ly(2),lz(2),lr(2)
      real*8 conv,doth1,doth2,conv2
      integer in

ccccccccccccccccccc
c>>  bond lengt
ccccccccccccccccccc

      do in=1,2
         lx(in)=0
         ly(in)=0
         lz(in)=0
         lr(in)=0
      enddo

      lx(1)=point_x(2)-point_x(1)
      ly(1)=point_y(2)-point_y(1)
      lz(1)=point_z(2)-point_z(1)
      lr(1)=sqrt(lx(1)**2+ly(1)**2+lz(1)**2)
      
      lx(2)=point_x(3)-point_x(2)
      ly(2)=point_y(3)-point_y(2)
      lz(2)=point_z(3)-point_z(2)
      lr(2)=sqrt(lx(2)**2+ly(2)**2+lz(2)**2)
 
cccccccccccccccccccc
c>>  theta value
cccccccccccccccccccc
      
      test_theta=0

      conv=180/3.14159 
      do in=1,1
         doth1=-1*(lx(in+1)*lx(in)+ly(in+1)*ly(in)+lz(in+1)*lz(in))
         doth2=doth1/(lr(in+1)*lr(in))
         if(doth2.gt.1)then
            doth2=1
         endif
         if(doth2.lt.-1)then
            doth2=-1
         endif
         test_theta=acos(doth2)*conv
      enddo

cccccccccccccccc

      return
      end

cccccccccccccccccccccccccccccccccccc
                        
ccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine readpdb(pdbid,numchain,
     &     chainid,
     &     usedchainlength,usedchain_restype_ca,
     &     usedchain_x_ca,usedchain_y_ca,usedchain_z_ca,
     &     usedchain_resseq_ca,usedchainatomlength,usedchain_atno,
     &     usedchain_restype,usedchain_atomtype,
     &     usedchain_x,usedchain_y,usedchain_z,
     &     usedchain_resseq,
     &     usedchain_reschainid,usedchain_atomchainid,
     &     flag)

      implicit none
      character*4 pdbid
      integer numchain
      character*1 chainid(10)
      integer chain_bg
      character*1 chain_bgicode
      integer chain_ed
      character*1 chain_edicode
      character*100 direction
      integer usedchainlength
      integer usedchainatomlength
      real*8 usedchain_x_ca(5000)
      real*8 usedchain_y_ca(5000)
      real*8 usedchain_z_ca(5000)
      character*3 usedchain_restype_ca(5000)
      integer usedchain_resseq_ca(5000)
      real*8 usedchain_x(50000)
      real*8 usedchain_y(50000)
      real*8 usedchain_z(50000)
      character*3 usedchain_restype(50000)
      character*4 usedchain_atomtype(50000)
      integer usedchain_resseq(50000)
      integer usedchain_atno(50000)
      character*1 usedchain_reschainid(5000)
      character*1 usedchain_atomchainid(50000)
      integer flag

      character*100 file4
      character*11 file
      integer*4 ij
      character*1 chainidused(10)
      integer wholechainlength
      integer in,i,im,j,in2
      character*6 atomserial
      character*6 recna
      character*4 atomname
      integer index_atom(100000)
      character*6 recna_atom(100000)
      integer atno_atom(100000)
      character*4 atna_atom(100000)
      character*1 altloc_atom(100000)
      character*3 resna_atom(100000)
      character*1 chainid_atom(100000)
      integer resseq_atom(100000)
      character*1 icode_atom(100000)
      real*8 x_atom(100000),y_atom(100000),z_atom(100000)
      real*8 occ_atom(100000),temf_atom(100000)
      character*4 seqid_atom(100000)
      character*1 altloc2
      character*3 resna2 
      character*2 chainid2
      character*4 resseq2
      character*1 icode2
      character*3 expmethod
      integer exp_flag

cccccccccccccccccccccccc

      usedchainlength=0
      do i=1,5000
         usedchain_restype_ca(i)='   '
         usedchain_x_ca(i)=0
         usedchain_y_ca(i)=0
         usedchain_z_ca(i)=0
         usedchain_reschainid(i)=' '
      enddo
      usedchainatomlength=0
      do i=1,50000
         usedchain_restype(i)='   '
         usedchain_atomtype(i)='    '
         usedchain_resseq(i)=0
         usedchain_x(i)=0
         usedchain_y(i)=0
         usedchain_z(i)=0
         usedchain_atno(i)=0
         usedchain_atomchainid(i)=' '
      enddo

cccccccccccccccccccccc
c  open file
cccccccccccccccccccccc

      file=pdbid//'.pdb'
     
ccccccccccccccccccccc
c read the chain id
ccccccccccccccccccccc

      do i=1,10
         chainidused(i)=chainid(i)
         if(chainidused(i).eq.'0')then
            chainidused(i)=' '
         endif
      enddo

 60   format (A1)

cccccccccccc

      direction='pdb/'

      ij=index(direction,' ')

      file4=direction(1:ij-1)//file

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c search the information about protein structure
ccccccccccccccccccccccccccccccccccccccccccccccccccccc

      open (unit=1,file=file4,status='old',
     &     access='sequential')

      im=0
      in=0

 900  read(1,5000) recna,atomserial,atomname
     &     ,altloc2,resna2, 
     &     chainid2,resseq2,icode2
            
      im=im+1

      if((recna.eq.'ENDMDL').OR.(recna.eq.'END   '))then
         goto 910
      endif

      if((recna.eq.'ATOM  ')
     &     .AND.((altloc2.eq.' ')
     &     .or.(altloc2.eq.'A')))then
         
         in=in+1
         index_atom(in)=im
     
      endif

      goto 900

 910  continue
            
      close(1)

 3000 format(A6)
 5000 format(A6,A6,A4,A1,A3,A2,A4,A1)

ccccccccccccccccccccccccc
c      print*,'in',in,'im',im
ccccccccccccccccccccccccc

      open (unit=1,file=file4,status='old',
     &     access='sequential')
         
      do i=1,index_atom(1)-1
         read(1,3000) recna
      enddo
           
      do i=1,in-1
               
         read(1,6000) recna_atom(i),atno_atom(i),atna_atom(i)
     &        ,altloc_atom(i),resna_atom(i), 
     &        chainid_atom(i),resseq_atom(i),icode_atom(i)
     &        ,x_atom(i),y_atom(i),z_atom(i),occ_atom(i),
     &        temf_atom(i),seqid_atom(i)
 6000    format(A6,I5,1x,A4,A1,A3,1X,A1,I4,A1,3x,3F8.3,
     &        F6.2,F6.2,6x,A4)
         
         do j=index_atom(i)+1,index_atom(i+1)-1
            read(1,3000) recna
         enddo
         
      enddo
      
      wholechainlength=in-1
         
      close(1)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     extract the information of chosen protein chain
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      in=0
      in2=0
      do j=1,numchain
         do i=1,wholechainlength
            if(chainid_atom(i).eq.chainidused(j))then
               in=in+1
               if(atna_atom(i).eq.' CA ')then
                  in2=in2+1
               endif
            endif
         enddo
      enddo
  
      if((in.gt.50000).or.(in2.gt.5000))then
         flag=1
         return
      else
         in=0
         in2=0

         do j=1,numchain
            do i=1,wholechainlength

               if(chainid_atom(i).eq.chainidused(j))then
                  
                  in=in+1
                  usedchain_atno(in)=atno_atom(i)
                  usedchain_x(in)=x_atom(i)
                  usedchain_y(in)=y_atom(i)
                  usedchain_z(in)=z_atom(i)
                  usedchain_restype(in)=resna_atom(i)
                  usedchain_resseq(in)=resseq_atom(i)
                  usedchain_atomtype(in)=atna_atom(i) 
                  usedchain_atomchainid(in)=chainid_atom(i)
                  
                  if(atna_atom(i).eq.' CA ')then
                     
                     in2=in2+1
                     usedchain_x_ca(in2)=x_atom(i)
                     usedchain_y_ca(in2)=y_atom(i)
                     usedchain_z_ca(in2)=z_atom(i)
                     usedchain_restype_ca(in2)=resna_atom(i)
                     usedchain_resseq_ca(in2)=resseq_atom(i)
                     usedchain_reschainid(in2)=chainid_atom(i)
                  endif               
                  
               endif
            enddo
         enddo
      
         usedchainatomlength=in
         usedchainlength=in2
                           
         return

      endif

ccccccccc

      end

ccccccccccccccccccccccccccccccccccccccccccccccccc            
            
