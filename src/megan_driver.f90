program main
   ! program:        MEGAN v3.3
   ! description:    biogenic VOCs emission model (by Alex Guenther)
   ! authors:        Alex Guenther, Ling Huang, Xuemei Wang, Jeff Willison, Hui Wang among others.
   ! programmed by:  Ramiro A. Espada (from Lakes Environmental Software), Hui Wang (UC Irvine)
   use mpi
   use netcdf   
   use datetime_module, only: datetime, timedelta, strptime!, secondsSinceEpoch
   use voc_mod   !megan module: (megan_voc)
   use co2_reader
   use, intrinsic :: iso_fortran_env, only: real64, int16
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   !use nox_mod   !megan module: (megan_nox)
   !use bdsnp    !megan module: (bdsnp_nox)
   
   implicit none

   INCLUDE 'tables/SPC_NOCONVER.EXT'
   !Variables:
   real, parameter :: FillValue = -999.0 
   integer :: iostat
   integer :: t,i,s,k,j!,k
   integer :: lai_num,laic_idx,laip_idx
   integer :: t_24,t_240,t_total

   !date-time vars:
   character(4) ::  yyyy,current_year                           
   character(3) ::  ddd,current_jday                            
   character(2) ::   mm,dd,hh,current_day="99",current_month="99"
   character(len=19) :: current_date
   type(datetime)    :: current_date_s, end_date_s, next_date_s
   type(datetime)    :: output_start_s, output_month_s
   integer(kind=8)   :: cur_sec_epoch

   !input variables:
   !==================soil variables========================
   !integer, allocatable, dimension(:,:)     :: stype                         !(x,y)   <- from wrfout
   !integer, allocatable, dimension(:,:)     :: arid,non_arid,landtype        !(x,y)   <- from prep_megan
   !real,    allocatable, dimension(:,:)     :: mapfac                !(x,y)   <- from wrfout
   !real,    allocatable, dimension(:,:)     :: smois,stemp
   !real   , allocatable, dimension(:,:)     :: cell_area                     !(x,y)   <- from prep_megan
   !real,    allocatable, dimension(:,:)     :: ndep,fert                 !(x,y,t) <- from prep_megan
   !==================above ground variables================
   integer(kind=8), allocatable, dimension(:) :: times
   real,    allocatable, dimension(:)       :: lon,lat
   real,    allocatable, dimension(:)       :: lon_local
   real,    allocatable, dimension(:,:)     :: temp,dtemp,ppfd
   real,    allocatable, dimension(:,:)     :: u10,v10,pres,rh
   real,    allocatable, dimension(:,:)     :: wind
   real,    allocatable, dimension(:,:)     :: rain
   real,    allocatable, dimension(:,:,:)   :: ctf,ef,ldf_in,lai             !(x,y,*) <- from prep_megan

   !intermediate vars:
   character(256)                      :: met_file,pmet_file
   character(256)                      :: out_file
   character(256)                      :: lai_file      !lai file
   real, allocatable, dimension(:,:)   :: temp_min,temp_max
   real, allocatable, dimension(:,:)   :: wind_max,temp24_avg,ppfd24_avg
   real, allocatable, dimension(:,:)   :: temp240_avg,ppfd240_avg
   real, allocatable, dimension(:,:,:) :: temp24, ppfd24, wind24
   real, allocatable, dimension(:,:,:) :: temp240,ppfd240

   !output vars:
   real, allocatable, dimension(:,:,:)   :: out_buffer
   real, allocatable, dimension(:,:,:,:) :: out_buffer_emis   !(x,y,1,nclass)
   real(real64), allocatable :: output_sum(:,:,:)
   integer, allocatable :: output_count(:,:,:)
   integer :: output_hours, output_samples=0, output_record=0
   character(len=7) :: output_month=''

   !megan namelist variables:
   character(len=19) :: start_date, end_date
   integer(kind=8)   :: cur_secd

   character(256)    :: met_file_path !path to met file
   character(256)    :: lai_file_path !path to lai file
   character(256)    :: pft_file_path !path to pft file
   character(256)    :: ef_file_path  !path to ef file
   character(256)    :: output_path   !path to output file
   character(len=16) :: output_frequency='24h' !1h, 3h, 6h, 24h or monthly
   

   character(256)    :: met_files     !global meteo files
   character(256)    :: lai_files     !global lai files
   character(256)    :: pft_files     !global pft file
   character(256)    :: ef_file       !emission factor file
   character(256)    :: co2_file       !emission factor file
   !flower and litter emission flag; Hui Wang
   logical           :: run_flower=.false., run_litter=.false.
   logical           :: run_co2   =.false.
   logical           :: real_spinup_met=.false., output_ef_file=.false.
   logical           :: diagnose=.false.,output_gamma=.false.
   integer           :: run_flower_flag,run_litter_flag,run_co2_flag
   integer           :: real_spinup_met_flag,output_ef_file_flag
   integer           :: diagnose_flag,output_gamma_flag

   character(3)      :: nlai='12'
   real              :: lai_scale_factor=0.1
   integer           :: ilen,nlat,nlon,ntime
   integer, allocatable :: co2_year(:), co2_month(:)
   real,    allocatable :: co2_avg(:)
   integer              :: co2_nrows
   real                 :: co2_value=420.
   !region defined parameters
   !mpi related
   integer :: ierr, rank, nprocs
   integer :: istart,iend

   !---read namelist variables and parameters
   namelist/megan_nl/start_date,end_date,&
                     met_file_path,lai_file_path,pft_file_path,&
                     ef_file_path,output_path,output_frequency,&
                     met_files, pft_files, ef_file, lai_files,&
                     nlai,lai_scale_factor,run_co2,co2_file,&
                     run_flower, run_litter,real_spinup_met,output_ef_file,output_gamma,diagnose

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)


   if (rank == 0)then
   !------------------------------------------------------------------
   print '(" =================================" )' 
   print '("     ╔╦╗╔═╗╔═╗╔═╗╔╗╔              " )'
   print '("     ║║║║╣ ║ ╦╠═╣║║║              " )'
   print '("     ╩ ╩╚═╝╚═╝╩ ╩╝╚╝ GLOBAL (v3.3)" )'
   print '(" =================================" )'
   !--------------------------------------------------------------------
   !reading megan namelist
   read(*,nml=megan_nl, iostat=iostat)
   !reading prep_megan namelist
   !read(*,nml=windowdefs, iostat=iostat)
   run_flower_flag      = flag2int(run_flower)
   run_litter_flag      = flag2int(run_litter)
   real_spinup_met_flag = flag2int(real_spinup_met)
   output_ef_file_flag  = flag2int(output_ef_file)
   run_co2_flag         = flag2int(run_co2)
   output_gamma_flag    = flag2int(output_gamma)
   diagnose_flag        = flag2int(diagnose)

   if( iostat /= 0 ) then
     call safe_mpi_exit("megan: failed to read namelist",iostat)
   end if
   call parse_output_frequency(output_frequency, output_hours)
   current_date_s = strptime(start_date,'%Y-%m-%d %H:%M:%S')
   end_date_s = strptime(end_date,'%Y-%m-%d %H:%M:%S')
   if (.not. current_date_s%isValid() .or. .not. end_date_s%isValid()) then
      call safe_mpi_exit('megan: invalid start_date or end_date', 1)
   end if
   if (end_date_s < current_date_s) then
      call safe_mpi_exit('megan: end_date precedes start_date', 1)
   end if
   if (current_date_s%getMinute() /= 0 .or. current_date_s%getSecond() /= 0 .or. &
       end_date_s%getMinute() /= 0 .or. end_date_s%getSecond() /= 0) then
      call safe_mpi_exit('megan: start_date and end_date must be on whole hours', 1)
   end if
   print '(A,A)', ' Output averaging frequency: ', trim(output_frequency)
   !prepare variables
   lai_num = atoi(nlai)
   
   print '(/" Read meteorology data information.")'
   current_date_s = strptime(start_date,'%Y-%m-%d %H:%M:%S')
   current_date=current_date_s%strftime("%Y %m %d %j %H")
   read(current_date ,*) yyyy,mm,dd,ddd,hh                !
   
   pmet_file=previous_filename(met_files,met_file_path, yyyy, mm )
   met_file =update_filename(met_files,met_file_path, yyyy, mm )
   write(*,'(A, A)') "Reading: ", trim(met_file)
   call get_grid(met_file,nlat,nlon,lat,lon)
       !read co2 file
       if(run_co2_flag == 1)then
       call read_co2_csv(co2_file, co2_year, co2_month,&
                         co2_avg, co2_nrows)
       end if
   end if!rank ==0

   !broadcast the information from the inputs
   call MPI_Bcast(met_file,     len(met_file),   MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(pmet_file,    len(pmet_file),  MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(start_date,   len(start_date), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(end_date,     len(end_date),   MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(output_path, len(output_path), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(output_frequency, len(output_frequency), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(output_hours, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(nlon,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(nlat,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)   
   call MPI_Bcast(lai_num, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)   
   call MPI_Bcast(run_flower_flag, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)   
   call MPI_Bcast(run_litter_flag, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)   
   call MPI_Bcast(run_co2_flag, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)   
   call MPI_Bcast(real_spinup_met_flag, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)   
   call MPI_Bcast(output_ef_file_flag,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)   
   call MPI_Bcast(output_gamma_flag,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)   
   call MPI_Bcast(diagnose_flag,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)   
   if (.not. allocated(lat)  ) allocate(lat(nlat))
   if (.not. allocated(lon)  ) allocate(lon(nlon))
   call MPI_Bcast(lat, nlat, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(lon, nlon, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
   if(run_co2_flag == 1)then
      call MPI_Bcast(co2_nrows,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (.not. allocated(co2_year)  )  allocate(co2_year(co2_nrows))
      if (.not. allocated(co2_month)  ) allocate(co2_month(co2_nrows))
      if (.not. allocated(co2_avg)  )   allocate(co2_avg(co2_nrows))
      call MPI_Bcast(co2_year, co2_nrows, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(co2_month, co2_nrows, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(co2_avg, co2_nrows, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
   end if 


   !================================================================
   !==========================Domain divided =======================
   !================================================================
   call divide_domain(nlon,rank,nprocs,istart,iend)
   ilen = iend - istart + 1
   write(*,'(A,I4,A,I4,A,I4)') "Rank:", rank, " istart:", istart, " iend:", iend
   call get_times_parallel(met_file,times)
   !--- Allocate variables
   allocate(lon_local(ilen))
   allocate(out_buffer_emis(ilen,nlat,1,nclass))
   allocate(out_buffer(ilen,nlat,nclass))
   allocate(output_sum(ilen,nlat,nclass), output_count(ilen,nlat,nclass))
   allocate( temp(ilen, nlat))
   allocate(dtemp(ilen, nlat))
   allocate(  u10(ilen, nlat))
   allocate(  v10(ilen, nlat))
   allocate( wind(ilen, nlat))
   allocate( pres(ilen, nlat))
   allocate( ppfd(ilen, nlat))
   allocate(   rh(ilen, nlat))
   !--- Allocate acclimation variables
   allocate( temp24(ilen, nlat, 24))
   allocate( ppfd24(ilen, nlat, 24))
   allocate( wind24(ilen, nlat, 24))
   allocate( temp240(ilen, nlat, 240))
   allocate( ppfd240(ilen, nlat, 240))
   !--- Min/Max/Averaged variables
   allocate( temp_min(ilen, nlat))
   allocate( temp_max(ilen, nlat))
   allocate( wind_max(ilen, nlat))
   allocate( temp24_avg(ilen, nlat))
   allocate( ppfd24_avg(ilen, nlat))
   allocate( temp240_avg(ilen, nlat))
   allocate( ppfd240_avg(ilen, nlat))
  
   !---vegetation data-------
   allocate(lai(ilen,nlat,lai_num+1))
   allocate(ctf(ilen,nlat,7))
   allocate(ef(ilen,nlat,19))
   allocate(ldf_in(ilen,nlat,4))
 
   out_buffer_emis=0.0
   output_sum=0.0_real64
   output_count=0
   temp=0.0  ;dtemp=0.0 ;u10=0.0; v10=0.0
   wind=0.0  ;pres=0.0  ;ppfd=0.0;rh=0.0
   temp24=0.0;ppfd24=0.0;wind24=0.0
   temp240=0.0;ppfd240=0.0
   temp_min=0.0;temp_max=0.0
   wind_max=0.0;temp24_avg =0.0;temp240_avg=0.0
   ppfd24_avg=0.0;ppfd240_avg=0.0
   lai = 0.0;ctf=0.0
   lon_local = lon(istart:iend) 
 
   if (rank .eq. 0) then
      print '(/" Init. temporal loop.. ")'
   end if
   !initiate time variables
   !====================================================
   current_date_s = strptime(start_date,'%Y-%m-%d %H:%M:%S')
   end_date_s     = strptime(  end_date,'%Y-%m-%d %H:%M:%S')
   current_date=current_date_s%strftime("%Y %m %d %j %H")
   read(current_date ,*) yyyy,mm,dd,ddd,hh                !
   current_year =YYYY
   current_month=MM
   current_day  =DD
   current_jday =DDD
   !====================================================
   !====================================================
   !initialize variables
   t_24           =1
   t_240          =1
   t_total        =1!total time step
   !====================================================
   if(real_spinup_met_flag .eq. 1) then
   call initiate_averaged_data(met_file,pmet_file,&
                                ilen,nlat,nlon,current_date_s,&
                                temp24,ppfd24,wind24,&
                                temp240,ppfd240,&
                                temp_max,temp_min,wind_max,&
                                temp24_avg,temp240_avg,&
                                ppfd24_avg,ppfd240_avg,FillValue)
   else
         temp24(:,:,:)  =288.0 !15deg Celsius (!CHECK VALUES!)
         ppfd24(:,:,:)  =400.  !              (!CHECK VALUES!)
         temp240(:,:,:) =288.0 !15deg Celsius (!CHECK VALUES!)
         ppfd240(:,:,:) =400.  !              (!CHECK VALUES!)
         wind24(:,:,:)  =2.0   !              (!CHECK VALUES!)
         do i=1,ilen
         do j=1,nlat
         temp_min(i,j)      = minval(temp24(i,j,:))
         temp_max(i,j)      = maxval(temp24(i,j,:))
         wind_max(i,j)      = maxval(wind24(i,j,:))
         temp24_avg(i,j)    = sum(temp24(i,j,:))/24.
         ppfd24_avg(i,j)    = sum(ppfd24(i,j,:))/24.
         temp240_avg(i,j)   = sum(temp240(i,j,:))/240. !time_len
         ppfd240_avg(i,j)   = sum(ppfd240(i,j,:))/240. !time_len
         end do
         end do

   end if
   call get_pft_data(pft_file_path,pft_files,&
                     nlat,nlon,ilen,ctf,yyyy)
   call get_ef_data(ef_file_path,ef_file,&
                     nlat,nlon,ilen,ctf(:,:,4:7),ef,ldf_in)
   if(output_ef_file_flag .eq. 1)then
         if(rank == 0) then
            out_file = trim(output_path)//"/"//"Gridded_MEGAN_EF_" // yyyy // ".nc"
            call create_output_ef_file(out_file,nlat,nlon,lat,lon)
         end if
         do k = 1, nclass
         call gather_output_2d(out_file,"EF_"//trim(mgn_spc(k)),ef(:,:,k),&
                            ilen,  nlat,  &
                            nlon, rank, nprocs)
         end do
         do k = 3, 6
         call gather_output_2d(out_file,"LDF_"//trim(mgn_spc(k)),ldf_in(:,:,k-2),&
                            ilen,  nlat,  &
                            nlon, rank, nprocs)
         end do
   end if
   if(run_co2_flag == 1)then
   co2_value=get_co2(co2_year, co2_month, co2_avg, atoi(yyyy), atoi(mm))
       if(rank ==0) then
       print*,rank,"CO2 concentration is ",co2_value
       end if
   end if
   !====================================================
   call get_lai_data(lai_file_path,lai_files,&
                     nlat,nlon,ilen,lai,&
                     lai_num,lai_scale_factor,yyyy)
   do while ( current_date_s <= end_date_s ) !temporal loop

      current_date=YYYY//"-"//MM//"-"//DD//" "//HH//":00:00"
      if (rank .eq. 0) print '(/"Current date ",a/)', trim(current_date)

      !update the meteorology file
      if (current_month .ne. mm) then
            pmet_file=met_file
            met_file=update_filename(met_files,met_file_path,yyyy,mm)
            call get_times_parallel(met_file,times)
            !update_co2_value
            if(run_co2_flag == 1)then
            co2_value=get_co2(co2_year, co2_month, co2_avg, atoi(yyyy), atoi(mm))
            if(rank ==0) then
            print*,"CO2 concentration is ",co2_value
            end if
            end if
      end if
      !update the lai file
      if (current_year .ne. yyyy) then
      call get_pft_data(pft_file_path,pft_files,&
                        nlat,nlon,ilen,ctf,yyyy)
      call get_ef_data(ef_file_path,ef_file,&
                     nlat,nlon,ilen,ctf(:,:,4:7),ef,ldf_in)
      call get_lai_data(lai_file_path,lai_files,&
                        nlat,nlon,ilen,lai,&
                        lai_num,lai_scale_factor,yyyy)
      end if

      !update time
      current_year =YYYY
      current_month=MM
      current_day  =DD
      current_jday =DDD
      cur_sec_epoch=current_date_s%secondsSinceEpoch()
      
      t=findloc(times == cur_sec_epoch, .true., 1)   !get index in time dimension of current date-time
      
      if (t == 0) then
         call safe_mpi_exit( "Error: current time not found!",t)
      end if

      call get_hourly_data(pmet_file,met_file, nlat, nlon, ilen,&
                              temp,dtemp,ppfd,pres,&
                              u10,v10,wind,rh,&
                              temp24,ppfd24,wind24,&
                              temp240,ppfd240,&
                              t, t_24, t_240, hh)

      !current/previous lai time index
      laip_idx = laiidx(ddd,lai_num)
      laic_idx = laip_idx + 1 
      !----------------------                                                             !run megan_voc
      out_buffer = 0.
      
      if(rank .eq. 0) print*,"   > Exec. megan_voc"

      call megan_voc(atoi(yyyy),atoi(ddd),atoi(hh),      & !date: year, julian day, hour.
             ilen,nlat,lon_local,lat,                    & !dimensions (ncols,nrows) & coordinates
             temp,ppfd,                                  & !Tmp.[ºK], PPFD [umol m-2 s-1]
             wind,pres,rh,                               & !Wind spd.[m/s], Press.[Pa], Humdty.[%]
             lai(:,:,laip_idx), lai(:,:,laic_idx),       & !LAI (past) [1], LAI (current) [1]
             ctf(:,:,1:6), ldf_in,                       & !Canopy type frac. [1],
             temp_max,temp_min,wind_max,                 & !max temp, min temp, max wind
             temp24_avg,temp240_avg,ppfd24_avg,ppfd240_avg,& !daily avg of temp & ppfd
             out_buffer,                                   &
             run_flower_flag,run_litter_flag,              &
             run_co2_flag, co2_value,diagnose_flag,FillValue                  ) 

      ! Average emissions after applying EF; gamma and diagnostics are averaged directly.
      if (diagnose_flag /= 1 .and. output_gamma_flag /= 1) then
         where (out_buffer /= FillValue .and. ef /= FillValue)
            out_buffer = out_buffer * ef
         elsewhere
            out_buffer = FillValue
         end where
      end if
      if (output_samples == 0) output_start_s = current_date_s
      call accumulate_output(out_buffer, output_sum, output_count)
      output_samples = output_samples + 1
      call get_averaged_data(rank,nlat,ilen,t_total,&
                                temp24,ppfd24,wind24,&
                                temp240,ppfd240,&
                                temp_max,temp_min,wind_max,&
                                temp24_avg,temp240_avg,&
                                ppfd24_avg,ppfd240_avg)
      !----------------------                                                              ! run megan_nox 
      !soil NO model:
      !if ( run_bdsnp ) then
      !    !call bdsnp_nox()                                
      !else
      !    call megan_nox(atoi(yyyy),atoi(ddd),atoi(hh),  & !date: year, julian day, hour.
      !           grid%nx,grid%ny,                        & !dimensions: (ncols nrows)
      !           lat,                                    & !latitude coordinates
      !           tmp,rain,                               & !temperature [ºK], precipitation rate [mm]
      !           lsm,stype,stemp,smois,                  & !land-surface-model, soil_type_clasification, soil temperature [ºK], soil mositure [m3/m3]
      !           ctf, laic,                               & !canopy type fraction [1], leaf-area-index [1]
      !           out_buffer(:,:,i_NO,atoi(HH))           ) !emision flux array [mole m-2 s-1]
      !endif
      

      next_date_s = current_date_s + timedelta(hours=1)
      if (output_interval_complete(current_date_s, end_date_s, output_hours)) then
         ! Flush at the interval boundary, month boundary, or final simulated hour.
         ! A partial first/last interval contains only the hours actually simulated.
         out_buffer_emis(:,:,1,:) = FillValue
         where (output_count > 0)
            out_buffer_emis(:,:,1,:) = real(output_sum / max(output_count, 1))
         end where
         if (output_month /= yyyy//'-'//mm) then
            output_month = yyyy//'-'//mm
            output_month_s = datetime(atoi(yyyy), atoi(mm), 1)
            output_record = 0
            out_file = trim(output_path)//'/Global_emis_MEGAN_'//output_month//'.nc'
            if (rank == 0) then
               call create_output_file(out_file,output_month//'-01',nlat,nlon,lat,lon,&
                                       output_gamma_flag,diagnose_flag,output_frequency)
            end if
         end if
         output_record = output_record + 1
         if (rank == 0) then
            call write_output_time(out_file, output_record, output_start_s, next_date_s, &
                                   output_month_s, output_samples)
         end if
         do k = 1, nclass
         if( diagnose_flag .eq. 1)then
         call gather_output_3d(out_file,trim(mgn_diag_var(k)),out_buffer_emis(:,:,:,k),&
                            ilen,  nlat, 1, &
                            nlon, rank, nprocs,0,output_record)
         else
         call gather_output_3d(out_file,trim(mgn_spc(k)),out_buffer_emis(:,:,:,k),&
                            ilen,  nlat, 1, &
                            nlon, rank, nprocs,output_gamma_flag,output_record)
         end if
         end do 
         output_sum = 0.0_real64
         output_count = 0
         output_samples = 0
      endif    

      current_date_s = next_date_s
      !Define next expected date
      current_date=current_date_s%strftime("%Y %m %d %j %H")
      read(current_date ,*) yyyy,mm,dd,ddd,hh                !
      !get current date in format: %Y-%m-%d %H:%M:%S
      t_total = t_total + 1
      !long-term time_loop
      t_24  = t_24 + 1
      t_240 = t_240 + 1

      if(t_24 .eq. 25) then
      t_24 = 1
      end if

      if(t_240 .eq. 241) then
      t_240 = 1
      end if
   end do!time loop
   
   if (allocated(out_buffer_emis)) deallocate(out_buffer_emis)
   if (allocated(out_buffer)) deallocate(out_buffer)
   if (allocated(output_sum)) deallocate(output_sum)
   if (allocated(output_count)) deallocate(output_count)
   if (allocated(temp)) deallocate( temp)
   if (allocated(dtemp))deallocate(dtemp)
   if (allocated(u10))  deallocate(  u10)
   if (allocated(v10))  deallocate(  v10)
   if (allocated(wind)) deallocate( wind)
   if (allocated(pres)) deallocate( pres)
   if (allocated(ppfd)) deallocate( ppfd)
   if (allocated(rh))   deallocate(   rh)
   if (allocated(lai)) deallocate( lai)
   if (allocated(ctf)) deallocate( ctf)
   if (allocated(ef)) deallocate( ef)
   if (allocated(ldf_in)) deallocate( ldf_in)
   !if (allocated(temp24)) deallocate( temp24)
   !if (allocated(ppfd24)) deallocate( ppfd24)
   !if (allocated(wind24)) deallocate( wind24)
   !if (allocated(temp240)) deallocate( temp240)
   !if (allocated(ppfd240)) deallocate( ppfd240)
   !if (allocated(temp_min)) deallocate( temp_min)
   !if (allocated(temp_max)) deallocate( temp_max)
   !if (allocated(wind_max)) deallocate( wind_max)
   !if (allocated(temp24_avg)) deallocate( temp24_avg)
   !if (allocated(temp240_avg)) deallocate( temp240_avg)
   !if (allocated(ppfd24_avg)) deallocate( ppfd24_avg)
   !if (allocated(ppfd240_avg)) deallocate( ppfd240_avg)
   
   if (rank==0)then
   print*, 'MEGAN finished succesfully.'
   end if
   call MPI_Finalize(ierr)
   !call safe_mpi_exit("megan: failed to read namelist",iostat)

contains

!UTILITIES -------------------------------------------------------------
   subroutine check(status)            !netcdf error-check function
     integer, intent(in) :: status
     if (status /= nf90_noerr) then
       write(*,*) nf90_strerror(status);
       call safe_mpi_exit('netcdf error', status)
     end if
   end subroutine check
   integer function laiidx(ddd,nlai)
    implicit none
    character(len=3), intent(in) :: ddd
    integer,          intent(in) :: nlai
    integer                      :: doy
    doy = atoi(ddd)

    laiidx = ceiling(real(doy)/real(365/nlai))

    !make sure the range of LAI
    if (laiidx > nlai) laiidx = nlai
    if (laiidx < 1) laiidx = 1 
   end function

   integer function atoi(str)
    implicit none
    character(len=*), intent(in) :: str
    integer :: iostat
    read(str,*,iostat=iostat) atoi
    if (iostat /= 0) then
        print *, "Error: Invalid integer format in string"
        atoi = 0
    end if
   end function atoi 
   integer function flag2int(logic_flag)
    implicit none
    logical, intent(in) :: logic_flag
    if (logic_flag) then
        flag2int = 1
    else
        flag2int = 0
    end if
    
   end function flag2int 
   character(len=20) function itoa(i)       !int -> string
      implicit none
      integer, intent(in) :: i
      !character(len=20) :: itoa
      write(itoa, '(i0)') i
      itoa = adjustl(itoa)
   end function
   character(len=16) function rtoa(r)       !real -> string
      implicit none
      real, intent(in) :: r
      !character(len=16) :: rtoa
      write(rtoa, '(F16.3)') r
      rtoa = adjustl(rtoa)
   end function
   
   pure function replace(string, s1, s2) result(str)
       !replace substring "s1" with "s2" on string
       implicit none
       character(*), intent(in)       :: string
       character(*), intent(in)       :: s1,s2
       character(len(string)+len(s2)) :: str          !not very elegant
       integer :: i,j,n,n1,n2!,dif
       str=string;n =len(str); n1=len(s1); n2=len(s2)
       if ( n1 == n2 ) then
          do i=1,n-n1
             if ( str(i:i+n1) == s1 ) str(i:i+n1) = s2
          end do
       else
           i=1
           do while ( i < len(trim(str)))
              if ( str(i:i+n1-1) == s1 ) then
                   str(i+n2:n)=str(i+n1:n)            !make space on str for replacement
                   str(i:i+n2-1) = s2                 !replace! 
              endif
             i=i+1
           enddo
       endif
   end function
   
   !input -----------------------------------------------------------------
   subroutine get_grid(meteo_file,nlat,nlon,lat,lon)
      implicit none
      character(256) ,   intent(in)    :: meteo_file
      integer,           intent(out)   :: nlat,nlon
      real, allocatable, intent(inout) :: lat(:),lon(:)
      integer :: ncid,dimid
      integer :: varid
   
      !get lat/lon from meteo file:
      call check(nf90_open(trim(meteo_file), nf90_nowrite, ncid ))
         !grid dimensions
         call check (nf90_inq_dimid(ncid,'latitude'  ,  dimid   ))
         call check (nf90_inquire_dimension(ncid, dimid,len=nlat ))
         call check (nf90_inq_dimid(ncid,'longitude',  dimid   ))
         call check (nf90_inquire_dimension(ncid, dimid,len=nlon ))
         nlat = nlat - 2 
         !!lat lon coordinates
         if (.not. allocated(lat)  ) allocate(lat(nlat))
         if (.not. allocated(lon)  ) allocate(lon(nlon))
         call check( nf90_inq_varid(ncid,'latitude' , varid))
         call check( nf90_get_var(ncId, varId, lat, start=[2], count=[nlat]))
         call check( nf90_inq_varid(ncid,'longitude', varid))
         call check( nf90_get_var(ncId, varId, lon) )!, start=[g%gx0,g%gy0], count=[g%nx,g%ny]))
      call check(nf90_close(ncid))
      !get rid of the polar points
   end subroutine
   
   subroutine get_times_parallel(meteo_file,times)
      implicit none
      character(250) ,intent(in)  :: meteo_file
      integer(kind=8), allocatable, intent(inout) :: times(:)
      integer :: ncid,dimid,varid,nt
      integer :: rank,ierr
  
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

      if (rank == 0)then
         call check(nf90_open(trim(meteo_file), nf90_nowrite, ncid ))
         call check(nf90_inq_dimid(ncid,'valid_time' ,  dimid   ))
         call check(nf90_inquire_dimension(ncid, dimId,len=nt ))
         
         if ( allocated(Times) ) deallocate(times)
         allocate(times(nt))
         call check(nf90_inq_varid(ncid,'valid_time', varid))
         call check(nf90_get_var(ncid, varid, times ))
         call check(nf90_close(ncid))
      end if
      call MPI_Bcast(nt, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (rank /= 0) then
         if (allocated(times)) deallocate(times)
         allocate(times(nt))
      endif

      call MPI_Bcast(times, nt, MPI_INTEGER8, 0, MPI_COMM_WORLD, ierr)
   end subroutine get_times_parallel
   
   character(250) function update_filename(file_names,file_path, yyyy, mm)
       !update file name according to flags. if no flags do nothing
       implicit none
       character(len=*),intent(in) :: file_names
       character(len=*),intent(in) :: file_path
       character(len=4),intent(in) :: YYYY
       character(len=2),intent(in) :: MM
       update_filename=trim(file_path)//'/'//trim(file_names)
       update_filename=replace(update_filename,"<time>", YYYY//"_"//MM)
   end function
   character(250) function previous_filename(file_names,file_path, yyyy, mm)
       !update file name according to flags. if no flags do nothing
       implicit none
       character(len=*),intent(in) :: file_names
       character(len=*),intent(in) :: file_path
       character(len=4),intent(in) :: YYYY
       character(len=2),intent(in) :: MM
       character(len=4) :: pyyyy
       character(len=2) :: pmm
       integer  :: yyyy_num, mm_num

       read(yyyy, *) yyyy_num
       read(mm, *) mm_num
       
       if(mm_num .eq. 1) then
          mm_num = 12
          yyyy_num = yyyy_num - 1
       else
          mm_num = mm_num - 1
       end if

       write(pyyyy, '(I4)') yyyy_num
       write(pmm, '(I2.2)') mm_num

       previous_filename=trim(file_path)//'/'//trim(file_names)
       previous_filename=replace(previous_filename,"<time>", pyyyy//"_"//pmm)
   end function
   
   !subroutine get_static_data(g) 
   !  implicit none
   !  type(grid_type) :: g
   !  integer         :: i,j,k
   !  integer         :: ncid, var_id
   !  character(len=10),dimension(19) :: ef_vars=["EF_ISOP   ", "EF_MBO    ", "EF_MT_PINE",&
   !                                              "EF_MT_ACYC", "EF_MT_CAMP", "EF_MT_SABI",&
   !                                              "EF_MT_AROM", "EF_NO     ", "EF_SQT_HR ",&
   !                                              "EF_SQT_LR ", "EF_MEOH   ", "EF_ACTO   ",&
   !                                              "EF_ETOH   ", "EF_ACID   ", "EF_LVOC   ",&
   !                                              "EF_OXPROD ", "EF_STRESS ", "EF_OTHER  ",&
   !                                              "EF_CO     "]
   !  character(len=5),dimension(4)   :: ldf_vars=["LDF03","LDF04","LDF05","LDF06"]
   !
   !  !Allocation of variables to use
   !  allocate(      cell_area(g%nx,g%ny)        )
   !  allocate(      ef(g%nx,g%ny,size( ef_vars)))
   !  allocate(  ldf_in(g%nx,g%ny,size(ldf_vars)))
   !  allocate(     ctf(g%nx,g%ny,NRTYP         ))
   !
   !  if ( run_bdsnp ) then
   !     allocate(    arid(g%nx,g%ny            ))
   !     allocate(non_arid(g%nx,g%ny            ))
   !     allocate(landtype(g%nx,g%ny            ))
   !     print '("   Reading: ",A50)',trim(static_file)//":LAND" !static_file !land_file !debug
   !     !LAND                                                                               
   !     call check(nf90_open(trim(static_file), nf90_write, ncid ))
   !        call check( nf90_inq_varid(ncid,'LANDTYPE', var_id ))
   !        call check( nf90_get_var(ncid, var_id, LANDTYPE ))

   !        call check( nf90_inq_varid(ncid,'ARID'    , var_id ))
   !        call check( nf90_get_var(ncid, var_id, ARID     ))

   !        call check( nf90_inq_varid(ncid,'NONARID' , var_id ))
   !        call check( nf90_get_var(ncid, var_id, NON_ARID ))
   !     call check(nf90_close(ncid))
   !  endif
   !  !CTS, EFS, LDF 
   !   print '("   Reading: ",A50)',trim(static_file) !debug
   !  call check(nf90_open(trim(static_file), nf90_write, ncid ))
   !      call check( nf90_inq_varid(ncid,'cell_area', var_id ))
   !      call check( nf90_get_var(ncid,var_id,cell_area))

   !      call check( nf90_inq_varid(ncid,'CTF', var_id ))
   !      call check( nf90_get_var(ncid,var_id,CTF,[1,1,1],[g%nx,g%ny,NRTYP]))

   !      call check( nf90_inq_varid(ncid,"EFS", var_id ))
   !      call check( nf90_get_var(ncid, var_id , ef ))   !new v3.3
   !    
   !      call check( nf90_inq_varid(ncid,"LDF", var_id ))
   !      call check( nf90_get_var(ncid, var_id , ldf_in ))  !new v3.3
   !  call check(nf90_close(ncid))
   !
   !  !From meteo:
   !  print '("   Reading: ",A50)',wrf_static_file !debug
   !  if (.not. allocated(stype))   allocate(  stype(g%nx,g%ny))
   !  if (.not. allocated(mapfac))  allocate( mapfac(g%nx,g%ny))
   !  call check(nf90_open(trim(wrf_static_file), nf90_write, ncid ))
   !      call check(nf90_inq_varid(ncid,'ISLTYP'  , var_id))
   !      call check(nf90_get_var(ncid, var_id,  stype, [1,1,1], [g%nx,g%ny,1]  ))
   !      call check(nf90_inq_varid(ncid,'MAPFAC_M', var_id)) 
   !      call check(nf90_get_var(ncid, var_id, mapfac, [1,1,1], [g%nx,g%ny,1]  ))
   !  call check(nf90_close(ncid))
   !end subroutine

   subroutine get_hourly_data(pmet_file,met_file, nlat, nlon, ilen,&
                              temp,dtemp,ppfd,pres,&
                              u10,v10,wind,rh,&
                              temp24,ppfd24,wind24,&
                              temp240,ppfd240,&
                               t, t_24, t_240, hh)
     use netcdf
     use mpi
     implicit none
   
     character(len=256), intent(in) :: met_file,pmet_file
     character(len=2),   intent(in) :: hh
     integer, intent(in)            :: nlat, nlon, ilen, t
     integer, intent(inout)         :: t_24, t_240
     real, intent(inout) :: temp(:,:),dtemp(:,:),u10(:,:),v10(:,:)
     real, intent(inout) :: wind(:,:),pres(:,:),ppfd(:,:),rh(:,:)
     real, intent(inout) :: temp24(:,:,:), ppfd24(:,:,:),wind24(:,:,:)
     real, intent(inout) :: temp240(:,:,:),ppfd240(:,:,:)
   
     !real,    allocatable :: temp_all(:,:), dtemp_all(:,:)
     !real,    allocatable :: u10_all(:,:), v10_all(:,:)
     !real,    allocatable :: pres_all(:,:)
     real,    allocatable :: ppfd_t1(:,:),ppfd_t2(:,:),data_all(:,:)
     real,    allocatable :: buf(:),bufrecv(:)
     integer, allocatable :: sendcounts(:), displs(:)
     integer :: ncid, var_id, ierr,  rank, nprocs
     integer :: ncid_p, var_id_p, ierr_p,dimid,nt_p
     integer :: i,j,p,is,ir,expected,istart,iend
     double precision :: t_start, t_end, t_elapsed
   
     ! Parallel info
     call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
     call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
   
     allocate(sendcounts(nprocs), displs(nprocs))
     do i = 0, nprocs - 1
       call compute_ilen(i, nlon, nprocs, sendcounts(i+1))
       sendcounts(i+1) = sendcounts(i+1) * nlat
     end do
   
     displs(1) = 0
     do i = 2, nprocs
       displs(i) = displs(i-1) + sendcounts(i-1)
     end do
   
     ! Local output (each process)
     if (rank == 0) t_start = MPI_Wtime()
     if (rank == 0) then
       allocate( data_all(nlon, nlat))
       write(*,'(A, A, A, I3 )') "Reading: ", trim(met_file)," at t=",t
       call check(nf90_open(trim(met_file), nf90_nowrite, ncid))
       call check(nf90_inq_varid(ncid, 't2m', var_id))
       call check(nf90_get_var(ncid, var_id, data_all, start=[1,2,t], count=[nlon,nlat,1]))
     end if
     ! temp
     call scatter_data(data_all, temp, nlon, nlat, ilen, rank,&
                      nprocs, sendcounts, displs, FillValue, ierr)
   
     if (rank == 0) then
       call check(nf90_inq_varid(ncid, 'd2m', var_id))
       call check(nf90_get_var(ncid, var_id, data_all, start=[1,2,t], count=[nlon,nlat,1]))
     end if
     ! dtemp
     call scatter_data(data_all, dtemp, nlon, nlat, ilen, rank,&
                      nprocs, sendcounts, displs, FillValue, ierr)
   
     if (rank == 0) then
       call check(nf90_inq_varid(ncid, 'u10', var_id))
       call check(nf90_get_var(ncid, var_id, data_all, start=[1,2,t], count=[nlon,nlat,1]))
     end if
     ! u10
     call scatter_data(data_all, u10, nlon, nlat, ilen, rank,&
                      nprocs, sendcounts, displs, FillValue, ierr)
   
     if (rank == 0) then
       call check(nf90_inq_varid(ncid, 'v10', var_id))
       call check(nf90_get_var(ncid, var_id, data_all, start=[1,2,t], count=[nlon,nlat,1]))
     end if
     ! v10
     call scatter_data(data_all, v10, nlon, nlat, ilen, rank,&
                      nprocs, sendcounts, displs, FillValue, ierr)
   
     if (rank == 0) then
       call check(nf90_inq_varid(ncid, 'sp', var_id))
       call check(nf90_get_var(ncid, var_id, data_all, start=[1,2,t], count=[nlon,nlat,1]))
     end if
     ! pressure
     call scatter_data(data_all, pres, nlon, nlat, ilen, rank,&
                      nprocs, sendcounts, displs, FillValue, ierr)
      
      if (rank == 0) then
           t_end = MPI_Wtime()
           print *, "Other var reading wall time:", t_end - t_start, "seconds"
           t_start = MPI_Wtime()
      end if
  
     if (rank == 0) then
       allocate( ppfd_t1(nlon, nlat))
       allocate( ppfd_t2(nlon, nlat))
       call check(nf90_inq_varid(ncid, 'ssrd', var_id))
       call check(nf90_get_var(ncid, var_id, ppfd_t2, start=[1,2,t], count=[nlon,nlat,1]))
       if(t .eq. 1) then
       write(*,'(A, A)') "Reading ssrd from the previous met file:", trim(pmet_file)
       call check(nf90_open(trim(pmet_file), nf90_nowrite, ncid_p))
       call check(nf90_inq_dimid(ncid_p,'valid_time' ,  dimid))
       call check(nf90_inquire_dimension(ncid_p, dimid,len=nt_p ))
       
       call check(nf90_inq_varid(ncid_p, 'ssrd', var_id_p))
       call check(nf90_get_var(ncid_p, var_id_p, ppfd_t1, start=[1,2,nt_p], count=[nlon,nlat,1]))
       call check(nf90_close(ncid_p))
       else
       call check(nf90_inq_varid(ncid, 'ssrd', var_id))
       call check(nf90_get_var(ncid, var_id, ppfd_t1, start=[1,2,t-1], count=[nlon,nlat,1]))
       end if
       
       if( hh .eq. "01") then
       data_all =ppfd_t2
       else
       data_all =ppfd_t2 - ppfd_t1
       end if
       call check(nf90_close(ncid))
       deallocate( ppfd_t1)
       deallocate( ppfd_t2)
     end if
   
     ! ppfd
     call scatter_data(data_all, ppfd, nlon, nlat, ilen, rank,&
                      nprocs, sendcounts, displs, FillValue, ierr)
      if (rank == 0) then
           t_end = MPI_Wtime()
           print *, "PPFD reading wall time:", t_end - t_start, "seconds"
      end if
     
     ! Derived quantities
     wind = sqrt(u10**2 + v10**2)
     
     !Ground Incident Radiation [J m-2] to PPFD (Photosynthetic Photon Flux Density [W m-2])
     ! ppfd = par   * 4.5   !par to ppfd
     ! par  = rgrnd * 0.5   !total rad to Photosyntetic Active Radiation (PAR)
     ppfd=ppfd*4.5*0.5/3600.
   
     ! water mixing ratio (Kg/Kg)
     rh = 0.622 * (6.112 * exp((17.625*(dtemp - 273.15))/(243.04 + (dtemp - 273.15)))) / &
           ((pres/100.0) - 6.112 * exp((17.625*(dtemp - 273.15))/(243.04 + (dtemp - 273.15)))) 
     !rh = 100*exp((17.625*(dtemp - 273.15))/(243.04 + (dtemp - 273.15))) / &
     !           exp((17.625*(temp  - 273.15))/(243.04 + (temp  - 273.15)))
     !calculate relative humidity based on T2 and DT2
     temp24(:,:,t_24) = temp
     ppfd24(:,:,t_24) = ppfd
     wind24(:,:,t_24) = wind
     
     temp240(:,:,t_240) = temp
     ppfd240(:,:,t_240) = ppfd
     !deallocate(buf) 
     if (rank == 0) then
     deallocate( data_all)
     !deallocate(dtemp_all)
     !deallocate(  u10_all)
     !deallocate(  v10_all)
     !deallocate( pres_all)
     !deallocate( ppfd_all)
     end if
   end subroutine get_hourly_data
!-----------------------------------------------------------------
   subroutine initiate_averaged_data(cmet_file,pmet_file,&
                                ilen,nlat,nlon,current_date_s,&
                                temp24,ppfd24,wind24,&
                                temp240,ppfd240,&
                                temp_max,temp_min,wind_max,&
                                temp24_avg,temp240_avg,&
                                ppfd24_avg,ppfd240_avg, FillValue)
     use datetime_module, only: datetime, timedelta, strptime!, secondsSinceEpoch
     implicit none
     character(len=256), intent(in) :: cmet_file,pmet_file
     integer,            intent(in) :: ilen,nlat,nlon
     type(datetime),     intent(in) :: current_date_s
     real, intent(in)            :: FillValue
     real, intent(inout) :: temp24(:,:,:), ppfd24(:,:,:),wind24(:,:,:)
     real, intent(inout) :: temp240(:,:,:),ppfd240(:,:,:)
     real, intent(inout) :: temp_max(:,:),temp_min(:,:),wind_max(:,:)
     real, intent(inout) :: temp24_avg(:,:),temp240_avg(:,:)
     real, intent(inout) :: ppfd24_avg(:,:),ppfd240_avg(:,:)
     
     
     real,  allocatable  :: ppfd_t1(:,:),ppfd_t2(:,:),data_all(:,:)
     type(datetime)      :: p240_date_s
     character(len=19)   :: p240_time
     integer             :: ncid, var_id, dimid, ierr, rank, nprocs
     integer             :: n,t,nt,f_flag,i,j
     integer(kind=8)     :: cur_sec_local_epoch
     integer(kind=8), allocatable, dimension(:) :: ptimes,ctimes
     integer, allocatable :: sendcounts(:), displs(:)
     character(4) ::  yyyy
     character(3) ::  ddd
     character(2) ::   mm,dd,hh
     character(len=19) :: current_date
     
     !get the time dimensions
     call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
     call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
     
     call get_times_parallel(pmet_file,ptimes)
     call get_times_parallel(cmet_file,ctimes)

     p240_date_s  = current_date_s - timedelta(hours=240)
    
     allocate(sendcounts(nprocs), displs(nprocs))
     do i = 0, nprocs - 1
       call compute_ilen(i, nlon, nprocs, sendcounts(i+1))
       sendcounts(i+1) = sendcounts(i+1) * nlat
     end do
   
     displs(1) = 0
     do i = 2, nprocs
       displs(i) = displs(i-1) + sendcounts(i-1)
     end do
     n = 1
     f_flag=0 

     do while ( p240_date_s < current_date_s ) !temporal loop
        current_date=p240_date_s%strftime("%Y %m %d %j %H")
        read(current_date ,*) yyyy,mm,dd,ddd,hh                !
        cur_sec_local_epoch=p240_date_s%secondsSinceEpoch()

        !go through the previous file
        t=findloc(ptimes == cur_sec_local_epoch, .true., 1)
        f_flag = 1
        
        !go through the current file
        if (t == 0) then
        t=findloc(ctimes == cur_sec_local_epoch, .true., 1)
        f_flag = 2
        end if
        
        if (t == 0) then
        call safe_mpi_exit( "Error: averaged time not found!",t)
        end if
     
        if (rank == 0) then
          print*,yyyy//"-"//mm//"-"//dd//":"//hh
          print*,"Reading n=",n 
          if(.not.allocated(data_all)) allocate( data_all(nlon, nlat))
          if(.not. allocated(ppfd_t1)) allocate(  ppfd_t1(nlon, nlat))
          if(.not. allocated(ppfd_t2)) allocate(  ppfd_t2(nlon, nlat))

          if (f_flag == 1) then
            write(*,'(A, A, A, I3 )') "Reading: ", trim(pmet_file)," at t=",t
            call check(nf90_open(trim(pmet_file), nf90_nowrite, ncid))
          else if (f_flag == 2) then
            write(*,'(A, A, A, I3 )') "Reading: ", trim(cmet_file)," at t=",t
            call check(nf90_open(trim(cmet_file), nf90_nowrite, ncid))
          else
            call safe_mpi_exit( "Error: averaged time not found!",t)
          end if
          call check(nf90_inq_varid(ncid, 't2m', var_id))
          call check(nf90_get_var(ncid, var_id, data_all, start=[1,2,t], count=[nlon,nlat,1]))
          call check(nf90_inq_varid(ncid, 'ssrd', var_id))
          call check(nf90_get_var(ncid, var_id, ppfd_t2, start=[1,2,t], count=[nlon,nlat,1]))
          call check(nf90_close(ncid))
        end if
        call scatter_data(data_all, temp240(:,:,n), nlon, nlat, ilen, rank,&
                      nprocs, sendcounts, displs, FillValue, ierr)
       
        if (rank == 0) then
          
          if (f_flag == 2) then!from the current file
              if(t .eq. 1) then
                 write(*,'(A, A)') "Reading ssrd from the previous met file:", trim(pmet_file)
                 call check(nf90_open(trim(pmet_file), nf90_nowrite, ncid))
                 call check(nf90_inq_dimid(ncid,'valid_time' ,  dimid))
                 call check(nf90_inquire_dimension(ncid, dimid,len=nt ))
                 
                 call check(nf90_inq_varid(ncid, 'ssrd', var_id))
                 call check(nf90_get_var(ncid, var_id, ppfd_t1, start=[1,2,nt], count=[nlon,nlat,1]))
                 call check(nf90_close(ncid))
              else
                 write(*,'(A, A)') "Reading ssrd from the previous met file:", trim(pmet_file)
                 call check(nf90_open(trim(cmet_file), nf90_nowrite, ncid))
                 call check(nf90_inq_varid(ncid, 'ssrd', var_id))
                 call check(nf90_get_var(ncid, var_id, ppfd_t1, start=[1,2,t-1], count=[nlon,nlat,1]))
                 call check(nf90_close(ncid))
              end if
          else if (f_flag ==1) then
                 write(*,'(A, A)') "Reading ssrd from the current met file:", trim(pmet_file)
                 call check(nf90_open(trim(pmet_file), nf90_nowrite, ncid))
                 call check(nf90_inq_varid(ncid, 'ssrd', var_id))
                 call check(nf90_get_var(ncid, var_id, ppfd_t1, start=[1,2,t-1], count=[nlon,nlat,1]))
                 call check(nf90_close(ncid))
          end if         

 
          if( hh .eq. "01") then
              !ppfd=ppfd*4.5*0.5/3600.
              data_all =(ppfd_t2)*2.25/3600.
          else
              data_all =(ppfd_t2 - ppfd_t1)*2.25/3600.
          end if
        end if!rank ==0
       call scatter_data(data_all, ppfd240(:,:,n), nlon, nlat, ilen, rank,&
                      nprocs, sendcounts, displs, FillValue, ierr)
   
       p240_date_s   = p240_date_s + timedelta(hours=1)
       n =n +1
     end do
     ! ppfd
     if( rank == 0)then
     deallocate( data_all)
     deallocate( ppfd_t2)
     deallocate( ppfd_t1)
     end if 
     temp24 = temp240(:,:,217:240)
     ppfd24 = ppfd240(:,:,217:240)
 
     do i=1,ilen
     do j=1,nlat
         if (temp24(i,j,1) /= FillValue)then
         temp_min(i,j)      = minval(temp24(i,j,:))
         temp_max(i,j)      = maxval(temp24(i,j,:))
         wind_max(i,j)      = maxval(wind24(i,j,:))
         temp24_avg(i,j)    = sum(temp24(i,j,:))/24.
         ppfd24_avg(i,j)    = sum(ppfd24(i,j,:))/24.
         temp240_avg(i,j)   = sum(temp240(i,j,:))/240. !time_len
         ppfd240_avg(i,j)   = sum(ppfd240(i,j,:))/240. !time_len
         else
         temp_min(i,j)      = FillValue!minval(temp24(i,j,:))
         temp_max(i,j)      = FillValue!maxval(temp24(i,j,:))
         wind_max(i,j)      = FillValue!maxval(wind24(i,j,:))
         temp24_avg(i,j)    = FillValue!sum(temp24(i,j,:))/24.
         ppfd24_avg(i,j)    = FillValue!sum(ppfd24(i,j,:))/24.
         temp240_avg(i,j)   = FillValue!sum(temp24(i,j,:))/24.
         ppfd240_avg(i,j)   = FillValue!sum(ppfd24(i,j,:))/24.
         end if
     end do
     end do
   end subroutine
      
!-----------------------------------------------------------------
   subroutine get_averaged_data(rank,nlat,nlon,t_total,&
                                temp24,ppfd24,wind24,&
                                temp240,ppfd240,&
                                temp_max,temp_min,wind_max,&
                                temp24_avg,temp240_avg,&
                                ppfd24_avg,ppfd240_avg)
      implicit none
      !character(len=3),intent(in) :: ddd                      
      !integer ::ncid,var_id
      integer, intent(in) :: rank
      integer, intent(in) :: nlat,nlon
      integer, intent(in) :: t_total
      real, intent(in) :: temp24(:,:,:),ppfd24(:,:,:),wind24(:,:,:)
      real, intent(in) :: temp240(:,:,:),ppfd240(:,:,:)
      real, intent(inout) :: temp_max(:,:),temp_min(:,:),wind_max(:,:)
      real, intent(inout) :: temp24_avg(:,:),temp240_avg(:,:)
      real, intent(inout) :: ppfd24_avg(:,:),ppfd240_avg(:,:)
      !integer ::t_24
      integer :: i,j 
      
      !if(rank .eq. 0) print '(" writing out file: ")'
      if(rank .eq. 0) print '(" Calculating 1-day/10-day averaged variables ")'
   
      ! 1day averaged data
      !if (t_total < 24) then 
      !    !initialize default variable values:
      !     temp24_avg(:,:) =288.0 !15deg Celsius (!CHECK VALUES!)
      !     ppfd24_avg(:,:) =600.  !              (!CHECK VALUES!)
      !     temp_min(:,:) =283.0 !10deg Celsius (!CHECK VALUES!)
      !     temp_max(:,:) =293.0 !20deg Celsius (!CHECK VALUES!)
      !     wind_max(:,:) =  2.0 !              (!CHECK VALUES!)
      !else
           do i=1,nlon
           do j=1,nlat
               if (temp24(i,j,1) /= FillValue)then
               temp_min(i,j)      = minval(temp24(i,j,:))
               temp_max(i,j)      = maxval(temp24(i,j,:))
               wind_max(i,j)      = maxval(wind24(i,j,:))
               temp24_avg(i,j)    = sum(temp24(i,j,:))/24.
               ppfd24_avg(i,j)    = sum(ppfd24(i,j,:))/24.
               temp240_avg(i,j)  = sum(temp240(i,j,:))/240. !time_len
               ppfd240_avg(i,j)  = sum(ppfd240(i,j,:))/240. !time_len
               else
               temp_min(i,j)      = FillValue!minval(temp24(i,j,:))
               temp_max(i,j)      = FillValue!maxval(temp24(i,j,:))
               wind_max(i,j)      = FillValue!maxval(wind24(i,j,:))
               temp24_avg(i,j)    = FillValue!sum(temp24(i,j,:))/24.
               ppfd24_avg(i,j)    = FillValue!sum(ppfd24(i,j,:))/24.
               temp240_avg(i,j)    = FillValue!sum(temp24(i,j,:))/24.
               ppfd240_avg(i,j)    = FillValue!sum(ppfd24(i,j,:))/24.
               end if
           end do
           end do
      !end if
      !if (t_total < 240) then
      !     temp240_avg =283.0; !10deg Celsius (!CHECK VALUES!)
      !     ppfd240_avg =400. ; !              (!CHECK VALUES!)
      !else
      !     do i=1,nlon
      !     do j=1,nlat
      !         if (temp24(i,j,1) /= FillValue)then
      !         else
      !         end if
      !     end do
      !     end do
      !end if  
 
      !if ( run_bdsnp ) then
      !   if (.not. allocated(fert)) then; allocate( fert(g%nx,g%ny));endif
      !   call check(nf90_open(trim(dynamic_file), nf90_write, ncid ))
      !   call check(   nf90_inq_varid(ncid,'FERT'//DDD, var_id ))
      !   call check( nf90_get_var(ncid, var_id , FERT ))
      !   call check(nf90_close(ncid))
      !endif
   
   end subroutine
!======================Plant Functional Type=======================
   subroutine get_pft_data(pft_path,pft_file_prefix,&
                           nlat,nlon,ilen,pft,yyyy)
     implicit none
     character(256),  intent(in)    :: pft_path,pft_file_prefix
     integer,         intent(in)    :: nlat,nlon,ilen
     real,         intent(inout)    :: pft(:,:,:)
     character(len=4),intent(in)    :: yyyy

     character(256)       :: pft_file
     integer              :: ncid,var_id
     integer              :: rank,nprocs
     logical              :: file_exists
     integer, allocatable :: sendcounts(:), displs(:)
     real,    allocatable :: data_all(:,:)!, btr_all(:,:),trop_all(:,:),ntr_all(:,:)
     !real,    allocatable :: shrub_all(:,:),crop_all(:,:),grass_all(:,:)
    
     
     
     ! Scatter each 2D var using 1D reshape
     call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
     call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
     allocate(sendcounts(nprocs), displs(nprocs))
     do i = 0, nprocs - 1
       call compute_ilen(i, nlon, nprocs, sendcounts(i+1))
       sendcounts(i+1) = sendcounts(i+1) * nlat
     end do
   
     displs(1) = 0
     do i = 2, nprocs
       displs(i) = displs(i-1) + sendcounts(i-1)
     end do


     if (rank == 0) then
        !file name with prefix string
        pft_file=trim(pft_path)//'/'//trim(pft_file_prefix)
        ! create the file name
        if ( index(pft_file,"<time>") /= 0 ) then
            pft_file=replace(pft_file,"<time>", yyyy)
        end if

        !check if the file existed  
        inquire(file=pft_file, exist=file_exists)
        if (.not. file_exists) then
           write(*, '(A, A)') "File missing: ", trim(pft_file)
           call safe_mpi_exit("PFT file error", 0)
        else
           write(*, '(A, A)') "Reading: ", trim(pft_file)
        end if

        if(.not. allocated (data_all))   allocate( data_all(nlon, nlat))
        call check(nf90_open(trim(pft_file), nf90_nowrite, ncid))  
        call check(nf90_inq_varid(ncid, 'NEEDTR', var_id))
        call check(nf90_get_var(ncid, var_id, data_all, start=[1,2], count=[nlon,nlat]))
     end if
     !needleleaf tree
     call scatter_data(data_all, pft(:,:,1), nlon, nlat, ilen, rank,&
                      nprocs, sendcounts, displs, FillValue, ierr)

     if (rank == 0) then
        call check(nf90_inq_varid(ncid, 'TROPBE', var_id))
        call check(nf90_get_var(ncid, var_id, data_all, start=[1,2], count=[nlon,nlat]))
     end if
     !tropical tree
     call scatter_data(data_all, pft(:,:,2), nlon, nlat, ilen, rank,&
                      nprocs, sendcounts, displs, FillValue, ierr)
        
     if (rank == 0) then
        call check(nf90_inq_varid(ncid, 'BROADTR', var_id))
        call check(nf90_get_var(ncid, var_id, data_all, start=[1,2], count=[nlon,nlat]))
     end if
     !broadleaf tree
     call scatter_data(data_all, pft(:,:,3), nlon, nlat, ilen, rank,&
                      nprocs, sendcounts, displs, FillValue, ierr)

     if (rank == 0) then
        call check(nf90_inq_varid(ncid, 'SHRUB', var_id))
        call check(nf90_get_var(ncid, var_id, data_all, start=[1,2], count=[nlon,nlat]))
     end if
     !shrub
     call scatter_data(data_all, pft(:,:,4), nlon, nlat, ilen, rank,&
                      nprocs, sendcounts, displs, FillValue, ierr)

     if (rank == 0) then
        call check(nf90_inq_varid(ncid, 'GRASS', var_id))
        call check(nf90_get_var(ncid, var_id, data_all, start=[1,2], count=[nlon,nlat]))
     end if
     !grass
     call scatter_data(data_all, pft(:,:,5), nlon, nlat, ilen, rank,&
                      nprocs, sendcounts, displs, FillValue, ierr)

     if (rank == 0) then
        call check(nf90_inq_varid(ncid, 'CROP', var_id))
        call check(nf90_get_var(ncid, var_id, data_all, start=[1,2], count=[nlon,nlat]))
     end if
     !crop
     call scatter_data(data_all, pft(:,:,6), nlon, nlat, ilen, rank,&
                      nprocs, sendcounts, displs, FillValue, ierr)


     if (rank == 0) then
        call check(nf90_inq_varid(ncid, 'TREE', var_id))
        call check(nf90_get_var(ncid, var_id, data_all, start=[1,2], count=[nlon,nlat]))
        call check(nf90_close(ncid))
     end if
     !tree (needle + tropical + broadleaf)
     call scatter_data(data_all, pft(:,:,7), nlon, nlat, ilen, rank,&
                      nprocs, sendcounts, displs, FillValue, ierr)
     
     deallocate(sendcounts)
     deallocate(displs)
     if (rank == 0) then
     if(allocated (data_all))   deallocate( data_all )
     end if 
   
   end subroutine get_pft_data
!======================emission factor=======================
   subroutine get_ef_data(ef_path,ef_file_prefix,&
                           nlat,nlon,ilen,pft,ef_out,ldf_out)
     implicit none
     character(256),  intent(in)    :: ef_path,ef_file_prefix
     integer,         intent(in)    :: nlat,nlon,ilen
     real,         intent(inout)    :: ef_out(:,:,:),ldf_out(:,:,:)
     real,         intent(in)       :: pft(:,:,:)

     character(256)       :: ef_file
     character(20), dimension(4) :: var_ef,var_ldf
     integer              :: ncid,var_id
     integer              :: rank,nprocs
     integer              :: v,i
     logical              :: file_exists
     integer, allocatable :: sendcounts(:), displs(:)
     real,    allocatable :: data_all(:,:),data_local(:,:)
     real,    allocatable :: vcf(:,:)
    
    
     var_ef  = [character(len=8) :: "SHRUB_EF","GRASS_EF","CROP_EF","TREE_EF"]
     var_ldf = [character(len=9) :: "SHRUB_LDF","GRASS_LDF","CROP_LDF","TREE_LDF"]
     

     if(.not. allocated (data_local))   allocate( data_local(ilen, nlat))
     if(.not. allocated (vcf))   allocate( vcf(ilen, nlat))
    
     vcf = 0.
     ef_out = 0.
     ldf_out = 0.
     do i=1,4
       vcf = vcf + pft(:,:,i)
     end do  
     ! Scatter each 2D var using 1D reshape
     call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
     call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
     allocate(sendcounts(nprocs), displs(nprocs))
     do i = 0, nprocs - 1
       call compute_ilen(i, nlon, nprocs, sendcounts(i+1))
       sendcounts(i+1) = sendcounts(i+1) * nlat
     end do
   
     displs(1) = 0
     do i = 2, nprocs
       displs(i) = displs(i-1) + sendcounts(i-1)
     end do

     ef_file=trim(ef_path)//'/'//trim(ef_file_prefix)
     !write(*,'(A, A)') "Reading EF file: ", trim(ef_file)
     if (rank == 0) then
        inquire(file=ef_file, exist=file_exists)
        if (.not. file_exists) then
           write(*, '(A, A)') "File missing: ", trim(ef_file)
           call safe_mpi_exit("EF file error", 0)
        else
           write(*, '(A, A)') "Reading EF file: ", trim(ef_file)
        end if
     end if
     
     do i=1,19!MEGAN species
     do v=1,4 !GF tpyes
     if (rank == 0) then
        if(.not. allocated (data_all))   allocate( data_all(nlon, nlat))
        call check(nf90_open(trim(ef_file), nf90_nowrite, ncid))  
        call check(nf90_inq_varid(ncid,var_ef(v), var_id))
        call check(nf90_get_var(ncid, var_id, data_all, start=[1,2,i], count=[nlon,nlat]))
        call check(nf90_close(ncid))
     end if
     !spread data
     call scatter_data(data_all, data_local, nlon, nlat, ilen, rank,&
                      nprocs, sendcounts, displs, FillValue, ierr)

     where (vcf /= 0.0)
     ef_out(:,:,i) = ef_out(:,:,i) + data_local * pft(:,:,v) / vcf
     !elsewhere
     !ef_out(:,:,i) = ef_out(:,:,i) + 0.0  ! 或者你可以选择跳过或加别的默认值
     end where
     !ef_out(:,:,i) = ef_out(:,:,i) + data_local*pft(:,:,v)/vcf
     end do 
     end do 
     
     do i=1,4 !MEGAN LDF species
     do v=1,4 !GF types
     if (rank == 0) then
        if(.not. allocated (data_all))   allocate( data_all(nlon, nlat))
        call check(nf90_open(trim(ef_file), nf90_nowrite, ncid))  
        call check(nf90_inq_varid(ncid,var_ldf(v), var_id))
        call check(nf90_get_var(ncid, var_id, data_all, start=[1,2,i], count=[nlon,nlat]))
        call check(nf90_close(ncid))
     end if
     !spread data
     call scatter_data(data_all, data_local, nlon, nlat, ilen, rank,&
                      nprocs, sendcounts, displs, FillValue, ierr)

     where (vcf /= 0.0)
     ldf_out(:,:,i) = ldf_out(:,:,i) + data_local*pft(:,:,v)/vcf
     end where
     end do 
     end do 
     
     
     deallocate(sendcounts)
     deallocate(displs)
     if(allocated (data_all))   deallocate( data_all )
     if(allocated (data_local)) deallocate( data_local )
     if(allocated (vcf))        deallocate( vcf)
   
   end subroutine get_ef_data
   
!========================LAI================================= 
   subroutine get_lai_data(lai_path,lai_file_prefix,&
                           nlat,nlon,ilen,lai,&
                           nlai,lai_scalefactor,yyyy)
     implicit none
     character(256), intent(in)  :: lai_path,lai_file_prefix
     integer, intent(in) :: nlat,nlon,ilen
     integer, intent(in) :: nlai
     real,    intent(in) :: lai_scalefactor
     real,    intent(inout) :: lai(:,:,:)
     character(len=4),intent(in)    :: yyyy

     character(256)      :: laic_file,lai_file_input
     character(4)        :: cyear, pyear
     character(len=3)    :: ct,pt
     integer             :: ncid,var_id
     integer             :: rank,nprocs
     integer             :: cy, py
     integer             :: indx
     logical             :: file_exists
     integer, allocatable :: sendcounts(:), displs(:)
     real,    allocatable :: laic_all(:,:),laip_all(:,:)
     integer(kind=1),    allocatable :: laic_tmp(:,:),laip_tmp(:,:)
     
     !file name with prefix string
     laic_file=trim(lai_path)//'/'//trim(lai_file_prefix)
     
     cy  = atoi(yyyy)
     cyear = yyyy
     py = cy - 1
     write(pyear, '(I4)') py
     
     ! Scatter each 2D var using 1D reshape
     call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
     call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
     allocate(sendcounts(nprocs), displs(nprocs))
     do i = 0, nprocs - 1
       call compute_ilen(i, nlon, nprocs, sendcounts(i+1))
       sendcounts(i+1) = sendcounts(i+1) * nlat
     end do
   
     displs(1) = 0
     do i = 2, nprocs
       displs(i) = displs(i-1) + sendcounts(i-1)
     end do


     do i =1,nlai+1
     
         if (rank == 0) then
            ! create the file name
            if( i .eq. 1 )then
               write(cyear,'(I4)') py
               indx = nlai
               write(ct, '(I3.3)') indx
            else
               write(cyear,'(I4)') cy 
               indx = i-1
               write(ct, '(I3.3)') indx
            end if
            if ( index(laic_file,"<time>") /= 0 ) then
                lai_file_input=replace(laic_file,"<time>", cyear//ct)
            end if

            !check if the file existed  
            inquire(file=lai_file_input, exist=file_exists)
            if (.not. file_exists) then
               write(*, '(A, A)') "File missing: ", trim(lai_file_input)
               call safe_mpi_exit("LAI file error", 0)
            else
               write(*, '(A, A)') "Reading: ", trim(lai_file_input)
            end if

            if(.not. allocated (laic_tmp)) allocate( laic_tmp(nlon, nlat))
            if(.not. allocated (laic_all)) allocate( laic_all(nlon, nlat))
            call check(nf90_open(trim(lai_file_input), nf90_nowrite, ncid))  
            call check(nf90_inq_varid(ncid, 'LAI', var_id))
            call check(nf90_get_var(ncid, var_id, laic_tmp, start=[1,2], count=[nlon,nlat]))
            call check(nf90_close(ncid))
            laic_all = real(laic_tmp, kind=4)
            laic_all = laic_all*lai_scalefactor
         endif
         ! laic
         call scatter_data(laic_all, lai(:,:,i), nlon, nlat, ilen, rank,&
                          nprocs, sendcounts, displs, FillValue, ierr)
     end do
     deallocate(sendcounts)
     deallocate(displs)
     
     if (allocated(laic_tmp)) deallocate(laic_tmp)
     if (allocated(laic_all)) deallocate(laic_all)
 
   
   end subroutine get_lai_data
!---------------create output file----------------------------- 
   subroutine create_output_ef_file(out_file,nlat,nlon,lat,lon)
     implicit none
     character(len=*), intent(in) :: out_file
     integer,          intent(in) :: nlat, nlon
     real, allocatable,intent(in) :: lat(:),lon(:)
   
     integer :: ncid, t_dim_id, x_dim_id, y_dim_id, str_dim_id
     integer :: var_id_lat, var_id_lon, var_id_area, var_id_time, k, t
     integer, allocatable :: var_id_mech(:)
   
     allocate(var_id_mech(nclass))
   
     print '(" writing gridded out file: ",a)', trim(out_file)
   
     ! create and define netcdf file
     !call check(nf90_create(trim(out_file), nf90_clobber, ncid))
     call check(nf90_create(trim(out_file), IOR(nf90_clobber, nf90_netcdf4), ncid))
 
     call check(nf90_def_dim(ncid, "lon",      nlon, x_dim_id))
     call check(nf90_def_dim(ncid, "lat",      nlat, y_dim_id))
   
     ! define variables explicitly
     call check(nf90_def_var(ncid, "lat", nf90_float, [y_dim_id], var_id_lat))
     call check(nf90_put_att(ncid, var_id_lat, "units", "degrees_north"))
     call check(nf90_put_att(ncid, var_id_lat, "long_name", "latitude"))
   
     call check(nf90_def_var(ncid, "lon", nf90_float, [x_dim_id], var_id_lon))
     call check(nf90_put_att(ncid, var_id_lon, "units", "degrees_east"))
     call check(nf90_put_att(ncid, var_id_lon, "long_name", "longitude"))
   
     do k = 1, nclass
       call check(nf90_def_var(ncid, "EF_"//trim(mgn_spc(k)), nf90_float, &
                               [x_dim_id, y_dim_id], var_id_mech(k)))
       call check(nf90_put_att(ncid, var_id_mech(k), "units", "nmole m-2 s-1"))
       call check(nf90_put_att(ncid, var_id_mech(k), "_FillValue", FillValue))
       call check(nf90_put_att(ncid, var_id_mech(k), "missing_value", FillValue))
       call check(nf90_put_att(ncid, var_id_mech(k), "var_desc", &
                               trim(mgn_spc(k))//" emission factor"))
     end do
     do k=3,6  
       call check(nf90_def_var(ncid, "LDF_"//trim(mgn_spc(k)), nf90_float, &
                               [x_dim_id, y_dim_id], var_id_mech(k)))
       call check(nf90_put_att(ncid, var_id_mech(k), "units", "fraction"))
       call check(nf90_put_att(ncid, var_id_mech(k), "_FillValue", FillValue))
       call check(nf90_put_att(ncid, var_id_mech(k), "missing_value", FillValue))
       call check(nf90_put_att(ncid, var_id_mech(k), "var_desc", &
                               trim(mgn_spc(k))//" Light Dependent Factor"))
     end do
 
     call check(nf90_enddef(ncid))
   
     ! write data directly (no reopening)
     call check(nf90_put_var(ncid, var_id_lat,  lat))
     call check(nf90_put_var(ncid, var_id_lon,  lon))
     call check(nf90_close(ncid))
   
     deallocate(var_id_mech)
   end subroutine create_output_ef_file

   !OUTPUT ----------------------------------------------------------------
   subroutine parse_output_frequency(frequency, hours)
      character(len=*), intent(inout) :: frequency
      integer, intent(out) :: hours
      integer :: i, code

      frequency = adjustl(frequency)
      do i = 1, len_trim(frequency)
         code = iachar(frequency(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) frequency(i:i) = achar(code + 32)
      end do
      select case (trim(frequency))
      case ('1h')
         hours = 1
      case ('3h')
         hours = 3
      case ('6h')
         hours = 6
      case ('24h')
         hours = 24
      case ('monthly')
         hours = 0
      case default
         call safe_mpi_exit('megan: output_frequency must be 1h, 3h, 6h, 24h or monthly', 1)
      end select
   end subroutine parse_output_frequency

   logical function output_interval_complete(current, last, hours) result(complete)
      type(datetime), intent(in) :: current, last
      integer, intent(in) :: hours
      type(datetime) :: next

      next = current + timedelta(hours=1)
      complete = next > last .or. next%getMonth() /= current%getMonth() .or. &
                 next%getYear() /= current%getYear()
      ! Fixed-hour windows are aligned to midnight UTC, independently of start_date.
      if (hours > 0) complete = complete .or. mod(next%getHour(), hours) == 0
   end function output_interval_complete

   subroutine accumulate_output(values, sums, counts)
      real, intent(in) :: values(:,:,:)
      real(real64), intent(inout) :: sums(:,:,:)
      integer, intent(inout) :: counts(:,:,:)

      ! Do not average fill values or non-finite values into valid grid cells.
      where (values /= FillValue .and. ieee_is_finite(values))
         sums = sums + real(values, real64)
         counts = counts + 1
      end where
   end subroutine accumulate_output

   subroutine write_output_time(out_file, record, first, next, origin, samples)
      character(len=*), intent(in) :: out_file
      integer, intent(in) :: record, samples
      type(datetime), intent(in) :: first, next, origin
      real(real64) :: bounds(2,1)
      integer :: ncid, varid

      bounds(1,1) = real(first%secondsSinceEpoch() - origin%secondsSinceEpoch(), real64)
      bounds(2,1) = real(next%secondsSinceEpoch() - origin%secondsSinceEpoch(), real64)
      call check(nf90_open(trim(out_file), nf90_write, ncid))
      call check(nf90_inq_varid(ncid, 'time', varid))
      call check(nf90_put_var(ncid, varid, [sum(bounds(:,1))/2.0_real64], start=[record], count=[1]))
      call check(nf90_inq_varid(ncid, 'time_bnds', varid))
      call check(nf90_put_var(ncid, varid, bounds, start=[1,record], count=[2,1]))
      call check(nf90_inq_varid(ncid, 'sample_count', varid))
      call check(nf90_put_var(ncid, varid, [samples], start=[record], count=[1]))
      call check(nf90_close(ncid))
   end subroutine write_output_time

   subroutine create_output_file(out_file,current_date,&
                                 nlat,nlon,lat,lon,&
                                 output_gamma_flag,&
                                 diagnose_flag,frequency)
     implicit none
     character(len=*), intent(in) :: out_file
     character(len=*), intent(in) :: current_date
     character(len=*), intent(in) :: frequency
     integer,          intent(in) :: nlat, nlon
     real, allocatable,intent(in) :: lat(:),lon(:)
     integer,          intent(in) :: output_gamma_flag,diagnose_flag
   
     integer :: ncid, t_dim_id, x_dim_id, y_dim_id, b_dim_id
     integer :: var_id_lat, var_id_lon, var_id_time, var_id_bounds, var_id_samples, k
     integer, allocatable :: var_id_mech(:)
   
     allocate(var_id_mech(nclass))
   
     print '(" writing out file: ",a)', trim(out_file)
   
     ! create and define netcdf file
     !call check(nf90_create(trim(out_file), nf90_clobber, ncid))
     call check(nf90_create(trim(out_file), IOR(nf90_clobber, nf90_netcdf4), ncid))
 
     call check(nf90_def_dim(ncid, "time", nf90_unlimited, t_dim_id))
     call check(nf90_def_dim(ncid, "nv",          2, b_dim_id))
     call check(nf90_def_dim(ncid, "lon",      nlon, x_dim_id))
     call check(nf90_def_dim(ncid, "lat",      nlat, y_dim_id))
   
     ! define variables explicitly
     call check(nf90_def_var(ncid, "lat", nf90_float, [y_dim_id], var_id_lat))
     call check(nf90_put_att(ncid, var_id_lat, "units", "degrees_north"))
     call check(nf90_put_att(ncid, var_id_lat, "long_name", "latitude"))
   
     call check(nf90_def_var(ncid, "lon", nf90_float, [x_dim_id], var_id_lon))
     call check(nf90_put_att(ncid, var_id_lon, "units", "degrees_east"))
     call check(nf90_put_att(ncid, var_id_lon, "long_name", "longitude"))
   
   
     call check(nf90_put_att(ncid, nf90_global, 'output_frequency', trim(frequency)))
     call check(nf90_put_att(ncid, nf90_global, 'averaging', &
                'Mean of valid hourly samples; partial intervals use only simulated hours.'))
     call check(nf90_def_var(ncid, "time", nf90_double, [t_dim_id], var_id_time))
     call check(nf90_put_att(ncid, var_id_time, "units", &
                "seconds since "//current_date//" 00:00:00 utc"))
     call check(nf90_put_att(ncid, var_id_time, "long_name", "midpoint of averaging interval"))
     call check(nf90_put_att(ncid, var_id_time, "axis", "T"))
     call check(nf90_put_att(ncid, var_id_time, "calendar", "standard"))
     call check(nf90_put_att(ncid, var_id_time, "standard_name", "time"))
     call check(nf90_put_att(ncid, var_id_time, 'bounds', 'time_bnds'))
     call check(nf90_def_var(ncid, 'time_bnds', nf90_double, [b_dim_id,t_dim_id], var_id_bounds))
     call check(nf90_def_var(ncid, 'sample_count', nf90_int, [t_dim_id], var_id_samples))
     call check(nf90_put_att(ncid, var_id_samples, 'long_name', 'number of simulated hourly samples in interval'))
     call check(nf90_put_att(ncid, var_id_samples, 'units', '1'))
   
     do k = 1, nclass
       if(diagnose_flag == 1)then
          call check(nf90_def_var(ncid, trim(mgn_diag_var(k)), nf90_float, &
                                  [x_dim_id, y_dim_id, t_dim_id], var_id_mech(k)))
          call check(nf90_put_att(ncid, var_id_mech(k), "_FillValue", FillValue))
          call check(nf90_put_att(ncid, var_id_mech(k), "missing_value", FillValue))
          call check(nf90_put_att(ncid, var_id_mech(k), "var_desc", &
                                  trim(mgn_diag_var(k))//" for diagnose"))
       else
          if( output_gamma_flag == 1) then
          call check(nf90_def_var(ncid, trim(mgn_spc(k))//"_GAMMA", nf90_short, &
                                  [x_dim_id, y_dim_id, t_dim_id], var_id_mech(k)))
          call check(nf90_put_att(ncid, var_id_mech(k), "units", "unit_less"))
          call check(nf90_put_att(ncid, var_id_mech(k), "scale_factor", 0.01))
          call check(nf90_put_att(ncid, var_id_mech(k), "_FillValue", int(nf90_fill_short, int16)))
          call check(nf90_put_att(ncid, var_id_mech(k), "missing_value", int(nf90_fill_short, int16)))
          call check(nf90_put_att(ncid, var_id_mech(k), "var_desc", &
                                  trim(mgn_spc(k))//" gamma value"))
          else
          call check(nf90_def_var(ncid, trim(mgn_spc(k)), nf90_float, &
                                  [x_dim_id, y_dim_id, t_dim_id], var_id_mech(k)))
          call check(nf90_put_att(ncid, var_id_mech(k), "units", "nmole m-2 s-1"))
          call check(nf90_put_att(ncid, var_id_mech(k), "_FillValue", FillValue))
          call check(nf90_put_att(ncid, var_id_mech(k), "missing_value", FillValue))
          call check(nf90_put_att(ncid, var_id_mech(k), "var_desc", &
                                  trim(mgn_spc(k))//" emission rate"))
          end if
       end if
       call check(nf90_put_att(ncid, var_id_mech(k), 'cell_methods', 'time: mean'))
     end do
   
     call check(nf90_enddef(ncid))
   
     ! write data directly (no reopening)
     call check(nf90_put_var(ncid, var_id_lat,  lat))
     call check(nf90_put_var(ncid, var_id_lon,  lon))
   
     call check(nf90_close(ncid))
   
     deallocate(var_id_mech)
   end subroutine create_output_file


   subroutine divide_domain(nlon_total, rank, nprocs, istart, iend)
      integer, intent(in) :: nlon_total, rank, nprocs
      integer, intent(out) :: istart, iend
      integer :: nbase, extra

      nbase = nlon_total / nprocs
      extra = mod(nlon_total, nprocs)

      if (rank < extra) then
        istart = rank * (nbase+1) + 1
        iend = istart + nbase
      else
        istart = extra * (nbase+1) + (rank-extra)*nbase + 1
        iend = istart + nbase - 1
      end if
   end subroutine divide_domain
   
   subroutine scatter_data(temp_all, temp, nlon, nlat, ilen, rank,&
                                   nprocs, sendcounts, displs, FillValue, ierr)
        real, intent(in) :: temp_all(:,:)
        real, intent(inout) :: temp(:,:)
        integer, intent(in) :: nlon, nlat, ilen, rank, nprocs
        integer, intent(in) :: sendcounts(:), displs(:)
        real, intent(in) :: FillValue
        integer, intent(out) :: ierr

        real, allocatable :: buf(:), bufrecv(:)
        integer :: is, ir, i, j, istart, iend, p

        allocate(buf(nlon*nlat))
        allocate(bufrecv(ilen * nlat))
        if (rank == 0) then
            is = 1
            do p = 0, nprocs - 1
                call divide_domain(nlon, p, nprocs, istart, iend)
                do j = 1, nlat
                    do i = istart, iend
                        buf(is) = temp_all(i, j)
                        is = is + 1
                    end do
                end do
            end do
        end if

        call MPI_Scatterv(buf, sendcounts, displs, MPI_REAL, bufrecv, ilen * nlat, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        !if (rank == 0) then
        !    call MPI_Scatterv(buf, sendcounts, displs, MPI_REAL, MPI_IN_PLACE, 0, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        !else
        !    call MPI_Scatterv(buf, sendcounts, displs, MPI_REAL, bufrecv, ilen*nlat, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        !end if


        ir = 1
        do j = 1, nlat
            do i = 1, ilen
                temp(i, j) = bufrecv(ir)
                ir = ir + 1
            end do
        end do

        where (isnan(temp))
            temp = FillValue
        end where

        !if (rank == 0) then
        !    deallocate(buf)
        !end if
        deallocate(buf)
        deallocate(bufrecv)
    end subroutine scatter_data
!==========================================================
   subroutine gather_output_3d(out_file,var_name,temp_in,&
                            ilen, nlat, nt, &
                            nlon_total, rank, nprocs, output_gamma_flag, record)
      implicit none

      character(len=*), intent(in) :: out_file
      character(len=*), intent(in) :: var_name
      real,             intent(in) :: temp_in(:,:,:)
      integer,          intent(in) :: ilen!, istart
      integer,          intent(in) :: nlat, nt, nlon_total
      integer,          intent(in) :: rank, nprocs
      integer,          intent(in) :: output_gamma_flag, record

      ! gather buffers
      real,    allocatable :: sendbuf(:), recvbuf(:), temp_global(:,:,:)
      real,    allocatable :: temp_local(:,:,:),temp_reordered(:,:,:)
      integer, allocatable :: recvcounts(:), displs(:)
      integer(kind=2),  allocatable :: temp_global_int(:,:,:)
      integer :: ierr,sendcount
      integer :: ncid, varid
      integer :: total_recv, expected
      integer :: p,k,i,j
      integer :: ix,istart,iend
      character(len=256) :: msg

      sendcount = ilen*nlat*nt
      !print *, "Rank", rank,"ilen",ilen, "nlat:", nlat, "nt:",nt
      !print *, "Rank", rank,"nlon",nlon_total, "nlat:", nlat, "nt:",nt

      if ( .not. allocated( temp_local)) allocate(temp_local(ilen,nlat,nt))
      if ( .not. allocated( sendbuf))    allocate(sendbuf(sendcount))
      if ( .not. allocated( recvcounts)) allocate(recvcounts(nprocs))
      if ( .not. allocated( displs))     allocate(    displs(nprocs))
      temp_local = temp_in
      !sendbuf = reshape(temp_local, [sendcount])
      ix = 1
      do k = 1, nt
        do j = 1, nlat
        do i = 1, ilen
             sendbuf(ix) = temp_local(i, j, k)
             ix = ix + 1
        end do
        end do
      end do

      if (rank == 0) then
        print '("Writing variable ",a)', trim(var_name)
        if(.not. allocated(recvbuf))     allocate(    recvbuf(nlon_total*nlat*nt))
        if(.not. allocated(temp_global)) allocate(temp_global(nlon_total,nlat,nt))
        ! reconstruct recvcounts and displs
        do i = 0, nprocs - 1
          call compute_ilen(i, nlon_total, nprocs, recvcounts(i+1))
          recvcounts(i+1) = recvcounts(i+1) * nlat * nt
        end do

        displs(1) = 0
        do i = 2, nprocs
          displs(i) = displs(i-1) + recvcounts(i-1)
        end do
      else
        recvcounts = 0
        displs = 0
      end if

      ! gather data
      call MPI_Gatherv(sendbuf, sendcount, MPI_REAL, &
                       recvbuf, recvcounts, displs, MPI_REAL, &
                       0, MPI_COMM_WORLD, ierr)


      if (rank == 0) then
      end if
      if (rank == 0) then
        !temp_global = reshape(recvbuf, [nlon_total, nlat, nt])
        ix = 1
        do p = 0, nprocs - 1
          call divide_domain(nlon_total,p,nprocs,istart,iend)
              do k = 1, nt
              do j = 1, nlat
              do i = istart, iend
              !  temp_global(istart_p + i - 1, j, k) = recvbuf(ix)
                temp_global(i, j, k) = recvbuf(ix)
                ix = ix + 1
              end do
              end do
              end do
        end do
        if(output_gamma_flag == 1)then
        if(.not. allocated(temp_global_int)) allocate(temp_global_int(nlon_total,nlat,nt))
        temp_global_int = int(nf90_fill_short, int16)
        where (temp_global /= FillValue .and. ieee_is_finite(temp_global))
           temp_global_int = int(temp_global*100., kind=int16)
        end where
        call check(nf90_open(trim(out_file), nf90_write, ncid ))
        call check(nf90_inq_varid(ncid,var_name//"_GAMMA",varid))
        call check(nf90_put_var(ncid, varid, temp_global_int, start=[1, 1, record], count=[nlon_total, nlat, nt]))
        call check(nf90_close(ncid))
        deallocate(temp_global_int)
        else
        call check(nf90_open(trim(out_file), nf90_write, ncid ))
        call check(nf90_inq_varid(ncid,var_name,varid))
        call check(nf90_put_var(ncid, varid, temp_global, start=[1, 1, record], count=[nlon_total, nlat, nt]))
        call check(nf90_close(ncid))
        end if
        deallocate(    recvbuf)
        deallocate(temp_global)
      end if
      deallocate(temp_local)
      deallocate(sendbuf)
      deallocate(recvcounts)
      deallocate(displs)
   end subroutine

!=======================================================
   subroutine gather_output_2d(out_file,var_name,temp_in,&
                            ilen, nlat, &
                            nlon_total, rank, nprocs)
      implicit none

      character(len=*), intent(in) :: out_file
      character(len=*), intent(in) :: var_name
      real,             intent(in) :: temp_in(:,:)
      integer,          intent(in) :: ilen!, istart
      integer,          intent(in) :: nlat, nlon_total
      integer,          intent(in) :: rank, nprocs

      ! gather buffers
      real,    allocatable :: sendbuf(:), recvbuf(:), temp_global(:,:)
      real,    allocatable :: temp_local(:,:),temp_reordered(:,:)
      integer, allocatable :: recvcounts(:), displs(:)
      integer :: ierr,sendcount
      integer :: ncid, varid
      integer :: total_recv, expected
      integer :: p,k,i,j
      integer :: ix,istart,iend
      character(len=256) :: msg

      sendcount = ilen*nlat
      !print *, "Rank", rank,"ilen",ilen, "nlat:", nlat, "nt:",nt
      !print *, "Rank", rank,"nlon",nlon_total, "nlat:", nlat, "nt:",nt

      if ( .not. allocated( temp_local)) allocate(temp_local(ilen,nlat))
      if ( .not. allocated( sendbuf))    allocate(sendbuf(sendcount))
      if ( .not. allocated( recvcounts)) allocate(recvcounts(nprocs))
      if ( .not. allocated( displs))     allocate(    displs(nprocs))
      temp_local = temp_in
      !sendbuf = reshape(temp_local, [sendcount])
      ix = 1
      do j = 1, nlat
      do i = 1, ilen
           sendbuf(ix) = temp_local(i, j)
           ix = ix + 1
      end do
      end do

      if (rank == 0) then
        print '("Writing variable ",a)', trim(var_name)
        if(.not. allocated(recvbuf))     allocate(    recvbuf(nlon_total*nlat))
        if(.not. allocated(temp_global)) allocate(temp_global(nlon_total,nlat))
        ! reconstruct recvcounts and displs
        do i = 0, nprocs - 1
          call compute_ilen(i, nlon_total, nprocs, recvcounts(i+1))
          recvcounts(i+1) = recvcounts(i+1) * nlat
        end do

        displs(1) = 0
        do i = 2, nprocs
          displs(i) = displs(i-1) + recvcounts(i-1)
        end do
      else
        recvcounts = 0
        displs = 0
      end if

      ! gather data
      call MPI_Gatherv(sendbuf, sendcount, MPI_REAL, &
                       recvbuf, recvcounts, displs, MPI_REAL, &
                       0, MPI_COMM_WORLD, ierr)


      if (rank == 0) then
      end if
      if (rank == 0) then
        !temp_global = reshape(recvbuf, [nlon_total, nlat, nt])
        ix = 1
        do p = 0, nprocs - 1
          call divide_domain(nlon_total,p,nprocs,istart,iend)
              do j = 1, nlat
              do i = istart, iend
              !  temp_global(istart_p + i - 1, j, k) = recvbuf(ix)
                temp_global(i, j) = recvbuf(ix)
                ix = ix + 1
              end do
              end do
        end do
        call check(nf90_open(trim(out_file), nf90_write, ncid ))
        call check(nf90_inq_varid(ncid,var_name,varid))
        call check(nf90_put_var(ncid, varid, temp_global, start=[1, 1], count=[nlon_total, nlat]))
        call check(nf90_close(ncid))
        deallocate(    recvbuf)
        deallocate(temp_global)
      end if
      deallocate(temp_local)
      deallocate(sendbuf)
      deallocate(recvcounts)
      deallocate(displs)
   end subroutine

!=======================================================



   subroutine compute_ilen(rk, nlon_total, nprocs, ilen_out)
      integer, intent(in) :: rk, nlon_total, nprocs
      integer, intent(out) :: ilen_out
      integer :: nbase, remain

      nbase = nlon_total / nprocs
      remain = mod(nlon_total, nprocs)
      if (rk < remain) then
        ilen_out = nbase + 1
      else
        ilen_out = nbase
      end if
   end subroutine compute_ilen

   subroutine safe_mpi_exit(msg, errcode)
      use mpi
      implicit none
      character(len=*), intent(in) :: msg
      integer, intent(in), optional :: errcode
      integer :: ierr, my_rank, code

      call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)

      if (my_rank == 0) then
        print *, 'FATAL ERROR: ', trim(msg)
      end if

      if (present(errcode)) then
        code = errcode
      else
        code = 1
      end if

      call MPI_Abort(MPI_COMM_WORLD, code, ierr)
   end subroutine safe_mpi_exit
   
end program main

