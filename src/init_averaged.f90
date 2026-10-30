
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
        print*,yyyy//"-"//mm//"-"//dd//":"//hh
        
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
         temp240_avg(i,j)    = FillValue!sum(temp24(i,j,:))/24.
         ppfd240_avg(i,j)    = FillValue!sum(ppfd24(i,j,:))/24.
         end if
     end do
     end do
   end subroutine
      
