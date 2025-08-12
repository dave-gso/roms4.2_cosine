#include "cppdefs.h"
      MODULE sedbio_output_mod

#if defined BIO_UMAINE15 && defined SEDBIO
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2025 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module defines/writes sediment biology variables               !
!  into output NetCDF files.                                           !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
      USE mod_sedbio
!
      USE def_var_mod,     ONLY : def_var
      USE nf_fwrite2d_mod, ONLY : nf_fwrite2d
# ifdef SOLVE3D
      USE nf_fwrite3d_mod, ONLY : nf_fwrite3d
# endif
# ifdef STATIONS
      USE extract_sta_mod, ONLY : extract_sta2d
#  ifdef SOLVE3D
      USE extract_sta_mod, ONLY : extract_sta3d
#  endif
# endif
      USE strings_mod,     ONLY : FoundError
!
      implicit none
!
      PUBLIC :: sedbio_def_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PUBLIC :: sedbio_def_pio
# endif
# ifdef STATIONS
      PUBLIC :: sedbio_def_station_nf90
      PUBLIC :: sedbio_wrt_station_nf90
#  if defined PIO_LIB && defined DISTRIBUTE
      PUBLIC :: sedbio_def_station_pio
      PUBLIC :: sedbio_wrt_station_pio
#  endif
# endif
      PUBLIC :: sedbio_wrt_nf90
# if defined PIO_LIB && defined DISTRIBUTE
      PUBLIC :: sedbio_wrt_pio
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE sedbio_def_nf90 (ng, model, ldef, VarOut, S,       &
     &                                bs3dgrd, t2dgrd)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: ldef, VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: bs3dgrd(:), t2dgrd(:)
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: Lrestart, got_var(NV)
!
      integer, parameter :: Natt = 25

      integer :: i, itrc, j, nvd3, nvd4, status, vtype
!
      real(r8) :: Aval(6)
!
# ifdef ADJOINT
      character (len=21)  :: Prefix
# else
      character (len=13)  :: Prefix
# endif
      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", sedbio_def_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Define sediment biology output variables.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=S(ng)%name
!
!  Set NetCDF variable external data type. The restart file is usually
!  in double precision.
!
      IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
        Lrestart=.TRUE.
        Vtype=NF_FRST
      ELSE
        Lrestart=.FALSE.
        Vtype=NF_FOUT
      END IF
!
!  Defined sedbio variables.
!
      DEFINE : IF (ldef) THEN
!
!  Set number of dimensions for output variables.
!
# if defined WRITE_WATER && defined MASKING
        nvd3=2
        nvd4=2
# else
        nvd3=3
        nvd4=4
# endif
!
!  Set long name prefix string.
!
# ifdef ADJOINT
!!      Prefix='time-averaged adjoint'
        Prefix='adjoint'
# else
!!      Prefix='time-averaged'
        Prefix=CHAR(32)                                   ! blank space
# endif
!
!  Initialize local information variable arrays.
!
        DO i=1,Natt
          DO j=1,LEN(Vinfo(1))
            Vinfo(i)(j:j)=' '
          END DO
        END DO
        DO i=1,6
          Aval(i)=0.0_r8
        END DO
!
!  Define sediment particulate matter concentrations.
!
        DO i=1,NPOM
	    IF (Lrestart.or.VarOut(idbPM(i),ng)) THEN
            Vinfo( 1)=Vname(1,idbPM(i))
            IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
              WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                      &
     &                                     TRIM(Vname(2,idbPM(i)))
            ELSE
              Vinfo( 2)=Vname(2,idbPM(i))
		END IF
            Vinfo( 3)=Vname(3,idbPM(i))
            Vinfo(14)=Vname(4,idbPM(i))
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idbPM(i),ng),r8)
!            Aval(5)=REAL(bs3dvar,r8)
            status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idbPM(i)),  &
     &                     vtype, nvd4, bs3dgrd, Aval, Vinfo, ncname,   &
     &                     SetFillVal = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
!  Define sediment particulate matter decay rates.
!
        DO i=1,NDR
          IF (Lrestart.or.VarOut(idbDR(i),ng)) THEN
            Vinfo( 1)=Vname(1,idbDR(i))
            IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
              WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                      &
     &                                     TRIM(Vname(2,idbDR(i)))
            ELSE
              Vinfo( 2)=Vname(2,idbDR(i))
		END IF
            Vinfo( 3)=Vname(3,idbDR(i))
            Vinfo(14)=Vname(4,idbDR(i))
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idbDR(i),ng),r8)
!            Aval(5)=REAL(bs3dvar,r8)
            status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idbDR(i)),  &
     &                     vtype, nvd4, bs3dgrd, Aval, Vinfo, ncname,   &
     &                     SetFillVal = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
!  Define sediment porewater concentrations.
!
        DO i=1,NPWC
          IF (Lrestart.or.VarOut(idbPW(i),ng)) THEN
            Vinfo( 1)=Vname(1,idbPW(i))
	      IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
              WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idbPW(i)))
            ELSE
              Vinfo( 2)=Vname(2,idbPW(i))
            END IF
            Vinfo( 3)=Vname(3,idbPW(i))
            Vinfo(14)=Vname(4,idbPW(i))
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(22)='coordinates'
!            Aval(5)=REAL(r2dvar,r8)
            Aval(5)=REAL(Iinfo(1,idbPW(i),ng),r8)
            status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idbPW(i)),      &
     &                   vtype, nvd3, t2dgrd, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
!  Define sediment-water fluxes.
!
        DO i=1,NSF
          IF (Lrestart.or.VarOut(idbSF(i),ng)) THEN
            Vinfo( 1)=Vname(1,idbSF(i))
            IF (S(ng)%ncid.eq.AVG(ng)%ncid) THEN
              WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idbSF(i)))
            ELSE
              Vinfo( 2)=Vname(2,idbSF(i))
            END IF
            Vinfo( 3)=Vname(3,idbSF(i))
            Vinfo(14)=Vname(4,idbSF(i))
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(22)='coordinates'
!            Aval(5)=REAL(r2dvar,r8)
            Aval(5)=REAL(Iinfo(1,idbSF(i),ng),r8)
            status=def_var(ng, model, S(ng)%ncid, S(ng)%Vid(idbSF(i)),      &
     &                   vtype, nvd3, t2dgrd, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO

      END IF DEFINE
!
!-----------------------------------------------------------------------
!  Otherwise, check existing output file and prepare for appending
!  data.
!-----------------------------------------------------------------------
!
      QUERY : IF (.not.ldef) THEN
!
!  Initialize local logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  sedbio variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            S(ng)%Vid(idtime)=var_id(i)
          END IF
          DO itrc=1,NPOM
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idbPM(itrc)))) THEN
              got_var(idbPM(itrc))=.TRUE.
              S(ng)%Vid(idbPM(itrc))=var_id(i)
            END IF
          END DO
          DO itrc=1,NPWC
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idbPW(itrc)))) THEN
              got_var(idbPW(itrc))=.TRUE.
              S(ng)%Vid(idbPW(itrc))=var_id(i)
            END IF
          END DO
          DO itrc=1,NDR
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idbDR(itrc)))) THEN
              got_var(idbDR(itrc))=.TRUE.
              S(ng)%Vid(idbDR(itrc))=var_id(i)
            END IF
          END DO
          DO itrc=1,NSF
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idbSF(itrc)))) THEN
              got_var(idbSF(itrc))=.TRUE.
              S(ng)%Vid(idbSF(itrc))=var_id(i)
            END IF
          END DO

        END DO
!
!  Check if output variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO i=1,NPOM
          IF (.not.got_var(idbPM(i)).and.                              &
     &        (Lrestart.or.VarOut(idbPM(i),ng))) THEN        
            IF (Master) WRITE (stdout,10) TRIM(Vname(1,idbPM(i))),     &
     &                                    TRIM(ncname)
                exit_flag=3
                RETURN
          END IF
        END DO
        DO i=1,NPWC
          IF (.not.got_var(idbPW(i)).and.                              &
     &        (Lrestart.or.VarOut(idbPW(i),ng))) THEN
            IF (Master) WRITE (stdout,10) TRIM(Vname(1,idbPW(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO i=1,NDR
          IF (.not.got_var(idbDR(i)).and.                              &
     &        (Lrestart.or.VarOut(idbDR(i),ng))) THEN
            IF (Master) WRITE (stdout,10) TRIM(Vname(1,idbDR(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO i=1,NSF
          IF (.not.got_var(idbSF(i)).and.                              &
     &        (Lrestart.or.VarOut(idbSF(i),ng))) THEN
            IF (Master) WRITE (stdout,10) TRIM(Vname(1,idbSF(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO

      END IF QUERY
!
  10  FORMAT (/,' SEDBIO_DEF_NF90 - unable to find variable: ',     &
     &        a,2x,' in output NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE sedbio_def_nf90
!
!***********************************************************************
      SUBROUTINE sedbio_wrt_nf90 (ng, model, tile,                  &
    &                                 LBi, UBi, LBj, UBj,               &
    &                                 VarOut, S)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: Lrestart
!
      integer :: gfactor, gtype, i, status
!
      real(dp) :: scale
!
      character (len=256) :: ncname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", sedbio_wrt_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out sediment biology output variables into
!  specified output NetCDF file.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=S(ng)%name
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
# if defined WRITE_WATER && defined MASKING
      gfactor=-1
# else
      gfactor=1
# endif
!
!  Set NetCDF variable external data type. The restart file is usually
!  in double precision.
!
      IF (S(ng)%ncid.eq.RST(ng)%ncid) THEN
        Lrestart=.TRUE.
      ELSE
        Lrestart=.FALSE.
      END IF
!
!  Write out sediment particulate matter concentrations.
!
      DO i=1,NPOM
        IF (Lrestart.or.VarOut(idbPM(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*bs3dvar
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idbPM(i),           &
     &                       S(ng)%Vid(idbPM(i)),                       &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, nspc, scale,        &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       SEDBIOL(ng) % sedPOM(:,:,:,i),             &
     &                       SetFillVal= .FALSE.)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idbPM(i))), S(ng)%Rindex,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out sediment porewater concentrations.
!
      DO i=1,NPWC
        IF (Lrestart.or.VarOut(idbPW(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idbPW(i),           &
     &                     S(ng)%Vid(idbPW(i)),                         &
     &                     S(ng)%Rindex, gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     SEDBIOL(ng) % sedPoreWaterCon(:,:,i))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idbPW(i))), S(ng)%Rindex,  &
     &                        TRIM(ncname)
		END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out sediment particulate matter decay rates.
!
      DO i=1,NDR
        IF (Lrestart.or.VarOut(idbDR(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*bs3dvar
          status=nf_fwrite3d(ng, model, S(ng)%ncid, idbDR(i),           &
     &                       S(ng)%Vid(idbDR(i)),                       &
     &                       S(ng)%Rindex, gtype,                       &
     &                       LBi, UBi, LBj, UBj, 1, nspc, scale,        &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       SEDBIOL(ng) % sedDecayRate(:,:,:,i),       &
     &                       SetFillVal= .FALSE.)
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idbDR(i))), S(ng)%Rindex,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out sediment-water fluxes.
!
      DO i=1,NSF
        IF (Lrestart.or.VarOut(idbSF(i),ng)) THEN
          scale=1.0_r8
          gtype=gfactor*r2dvar
          status=nf_fwrite2d(ng, model, S(ng)%ncid, idbSF(i),           &
     &                     S(ng)%Vid(idbSF(i)),                         &
     &                     S(ng)%Rindex, gtype,                         &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     SEDBIOL(ng) % sedFlux(:,:,i))
          IF (FoundError(status, nf90_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idbSF(i))), S(ng)%Rindex,  &
     &                        TRIM(ncname)
		END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO

!
  10  FORMAT (/," SEDBIO_WRT_NF90 - error while writing variable '",&
     &        a,"', time record = ",i0,/,11x,'into file: ',a)
!
      RETURN
      END SUBROUTINE sedbio_wrt_nf90

# if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE sedbio_def_pio (ng, model, ldef, VarOut, S,        &
     &                               bs3dgrd, t2dgrd)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: ldef, VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: bs3dgrd(:), t2dgrd(:)
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: Lrestart, got_var(NV)
!
      integer, parameter :: Natt = 25

      integer :: FileH, i, itrc, j, nvd3, nvd4, status, vtype
!
      real(r8) :: Aval(6)
!
#  ifdef ADJOINT
      character (len=21)  :: Prefix
#  else
      character (len=13)  :: Prefix
#  endif
      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", sedbio_def_pio"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Define sediment biology output variables.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=S(ng)%name
      FileH=ABS(S(ng)%pioFile%fh)
!
!  Set NetCDF variable external data type. The restart file is usually
!  in double precision.
!
      IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
        Lrestart=.TRUE.
        Vtype=PIO_FRST
      ELSE
        Lrestart=.FALSE.
        Vtype=PIO_FOUT
      END IF
!
!  Defined sedbio variables.
!
      DEFINE : IF (ldef) THEN
!
!  Set number of dimensions for output variables.
!
#  if defined WRITE_WATER && defined MASKING
        nvd3=2
        nvd4=2
#  else
        nvd3=3
        nvd4=4
#  endif
!
!  Set long name prefix string.
!
#  ifdef ADJOINT
!!      Prefix='time-averaged adjoint'
        Prefix='adjoint'
#  else
!!      Prefix='time-averaged'
        Prefix=CHAR(32)                                   ! blank space
#  endif
!
!  Initialize local information variable arrays.
!
        DO i=1,Natt
          DO j=1,LEN(Vinfo(1))
            Vinfo(i)(j:j)=' '
          END DO
        END DO
        DO i=1,6
          Aval(i)=0.0_r8
        END DO
!
!  Define sediment particulate matter concentrations.
!
        DO i=1,NPOM
	    IF (Lrestart.or.VarOut(idbPM(i),ng)) THEN
            Vinfo( 1)=Vname(1,idbPM(i))
            IF (FileH.eq.ABS(AVG(ng)%pioFile%fh)) THEN
              WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                      &
     &                                     TRIM(Vname(2,idbPM(i)))
            ELSE
              Vinfo( 2)=Vname(2,idbPM(i))
		END IF
            Vinfo( 3)=Vname(3,idbPM(i))
            Vinfo(14)=Vname(4,idbPM(i))
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idbPM(i),ng),r8)
            S(ng)%pioVar(idbPM(i))%dkind=PIO_FOUT
            S(ng)%pioVar(idbPM(i))%gtype=bs3dvar
!
            status=def_var(ng, model, S(ng)%pioFile,                    &
     &                     S(ng)%pioVar(idbPM(i))%vd,                   &
     &                     vtype, nvd4, bs3dgrd, Aval, Vinfo, ncname,   &
     &                     SetFillVal = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
!  Define sediment particulate matter decay rates.
!
        DO i=1,NDR
          IF (Lrestart.or.VarOut(idbDR(i),ng)) THEN
            Vinfo( 1)=Vname(1,idbDR(i))
            IF (FileH.eq.ABS(AVG(ng)%pioFile%fh)) THEN
              WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                      &
     &                                     TRIM(Vname(2,idbDR(i)))
            ELSE
              Vinfo( 2)=Vname(2,idbDR(i))
		END IF
            Vinfo( 3)=Vname(3,idbDR(i))
            Vinfo(14)=Vname(4,idbDR(i))
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idbDR(i),ng),r8)
            S(ng)%pioVar(idbDR(i))%dkind=PIO_FOUT
            S(ng)%pioVar(idbDR(i))%gtype=bs3dvar
            status=def_var(ng, model, S(ng)%pioFile,                    &  
                           S(ng)%pioVar(idbDR(i))%vd,                   &
     &                     vtype, nvd4, bs3dgrd, Aval, Vinfo, ncname,   &
     &                     SetFillVal = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
!  Define sediment porewater concentrations.
!
        DO i=1,NPWC
          IF (Lrestart.or.VarOut(idbPW(i),ng)) THEN
            Vinfo( 1)=Vname(1,idbPW(i))
	      IF (FileH.eq.ABS(AVG(ng)%pioFile%fh)) THEN
              WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idbPW(i)))
            ELSE
              Vinfo( 2)=Vname(2,idbPW(i))
            END IF
            Vinfo( 3)=Vname(3,idbPW(i))
            Vinfo(14)=Vname(4,idbPW(i))
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idbPW(i),ng),r8)
            S(ng)%pioVar(idbPW(i))%dkind=vtype
            S(ng)%pioVar(idbPW(i))%gtype=r2dvar
!
            status=def_var(ng, model, S(ng)%pioFile,                    &
     &                   S(ng)%PioVar(idbPW(i))%vd,                     &
     &                   vtype, nvd3, t2dgrd, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
!  Define sediment-water fluxes.
!
        DO i=1,NSF
          IF (Lrestart.or.VarOut(idbSF(i),ng)) THEN
            Vinfo( 1)=Vname(1,idbSF(i))
            IF (FileH.eq.ABS(AVG(ng)%pioFile%fh)) THEN
              WRITE (Vinfo( 2),'(a,1x,a)') Prefix,                        &
     &                                   TRIM(Vname(2,idbSF(i)))
            ELSE
              Vinfo( 2)=Vname(2,idbSF(i))
            END IF
            Vinfo( 3)=Vname(3,idbSF(i))
            Vinfo(14)=Vname(4,idbSF(i))
            Vinfo(16)=Vname(1,idtime)
#  if defined WRITE_WATER && defined MASKING
            Vinfo(20)='mask_rho'
#  endif
            Vinfo(22)='coordinates'
            Aval(5)=REAL(Iinfo(1,idbSF(i),ng),r8)
            S(ng)%pioVar(idbSF(i))%dkind=vtype
            S(ng)%pioVar(idbSF(i))%gtype=r2dvar
            status=def_var(ng, model, S(ng)%pioFile,                    &
     &                   S(ng)%pioVar(idbSF(i))%vd,                     &
     &                   vtype, nvd3, t2dgrd, Aval, Vinfo, ncname,      &
     &                   SetFillVal = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO

      END IF DEFINE
!
!-----------------------------------------------------------------------
!  Otherwise, check existing output file and prepare for appending
!  data.
!-----------------------------------------------------------------------
!
      QUERY : IF (.not.ldef) THEN
!
!  Initialize local logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  sediment biology variables. Get variable IDs.
!
        DO i=1,n_var
          DO itrc=1,NPOM
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idbPM(itrc)))) THEN
              got_var(idbPM(itrc))=.TRUE.
              S(ng)%pioVar(idbPM(itrc))%vd=var_desc(i)
              S(ng)%pioVar(idbPM(itrc))%dkind=PIO_FOUT
              S(ng)%pioVar(idbPM(itrc))%gtype=bs3dvar
            END IF
          END DO
          DO itrc=1,NPWC
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idbPW(itrc)))) THEN
              got_var(idbPW(itrc))=.TRUE.
              S(ng)%pioVar(idbPW(itrc))%vd=var_desc(i)
              S(ng)%pioVar(idbPW(itrc))%dkind=PIO_FOUT
              S(ng)%pioVar(idbPW(itrc))%gtype=r2dvar
            END IF
          END DO
          DO itrc=1,NDR
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idbDR(itrc)))) THEN
              got_var(idbDR(itrc))=.TRUE.
              S(ng)%pioVar(idbDR(itrc))%vd=var_desc(i)
              S(ng)%pioVar(idbDR(itrc))%dkind=PIO_FOUT
              S(ng)%pioVar(idbDR(itrc))%gtype=bs3dvar
            END IF
          END DO
          DO itrc=1,NSF
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idbSF(itrc)))) THEN
              got_var(idbSF(itrc))=.TRUE.
              S(ng)%pioVar(idbSF(itrc))%vd=var_desc(i)
              S(ng)%pioVar(idbSF(itrc))%dkind=PIO_FOUT
              S(ng)%pioVar(idbSF(itrc))%gtype=r2dvar
            END IF
          END DO        
        END DO
!
!  Check if output variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,10) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO i=1,NPOM
          IF (.not.got_var(idbPM(i)).and.                              &
     &        (Lrestart.or.VarOut(idbPM(i),ng))) THEN        
            IF (Master) WRITE (stdout,10) TRIM(Vname(1,idbPM(i))),     &
     &                                    TRIM(ncname)
                exit_flag=3
                RETURN
          END IF
        END DO
        DO i=1,NPWC
          IF (.not.got_var(idbPW(i)).and.                              &
     &        (Lrestart.or.VarOut(idbPW(i),ng))) THEN
            IF (Master) WRITE (stdout,10) TRIM(Vname(1,idbPW(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO i=1,NDR
          IF (.not.got_var(idbDR(i)).and.                              &
     &        (Lrestart.or.VarOut(idbDR(i),ng))) THEN
            IF (Master) WRITE (stdout,10) TRIM(Vname(1,idbDR(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO i=1,NSF
          IF (.not.got_var(idbSF(i)).and.                              &
     &        (Lrestart.or.VarOut(idbSF(i),ng))) THEN
            IF (Master) WRITE (stdout,10) TRIM(Vname(1,idbSF(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO

      END IF QUERY
!
  10  FORMAT (/,' SEDBIO_DEF_PIO - unable to find variable: ',      &
     &        a,2x,' in output NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE sedbio_def_pio
!
!***********************************************************************
      SUBROUTINE sedbio_wrt_pio (ng, model, tile,                   &
    &                                LBi, UBi, LBj, UBj,                &
    &                                VarOut, S)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: Lrestart
!
      integer :: FileH, gfactor, gtype, i, status
!
      real(dp) :: scale
!
      character (len=256) :: ncname

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", sedbio_wrt_pio"
!
      TYPE (IO_desc_t), pointer :: ioDesc
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out sediment biology output variables into
!  specified NetCDF file.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=S(ng)%name
      FileH=ABS(S(ng)%pioFile%fh)
!
!  Set grid type factor to write full (gfactor=1) fields or water
!  points (gfactor=-1) fields only.
!
#  if defined WRITE_WATER && defined MASKING
      gfactor=-1
#  else
      gfactor=1
#  endif
!
!  Set NetCDF variable external data type. The restart file is usually
!  in double precision.
!
      IF (FileH.eq.ABS(RST(ng)%pioFile%fh)) THEN
        Lrestart=.TRUE.
      ELSE
        Lrestart=.FALSE.
      END IF
!
!  Write out sediment particulate matter concentrations.
!
      DO i=1,NPOM
        IF (Lrestart.or.VarOut(idbPM(i),ng)) THEN
          scale=1.0_r8
          IF (S(ng)%pioVar(idbPM(i))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_bs3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_bs3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, S(ng)%pioFile, idbPM(i),        &
     &                       S(ng)%pioVar(idbPM(i)),                    &
     &                       S(ng)%Rindex,                              &
     &                       ioDesc,                                    &
     &                       LBi, UBi, LBj, UBj, 1, nspc, scale,        &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       SEDBIOL(ng) % sedPOM(:,:,:,i),             &
     &                       SetFillVal= .FALSE.)
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idbPM(i))), S(ng)%Rindex,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out sediment porewater concentrations.
!
      DO i=1,NPWC
        IF (Lrestart.or.VarOut(idbPW(i),ng)) THEN
          scale=1.0_r8
	    IF (S(ng)%pioVar(idbPW(i))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r2dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_r2dvar(ng)
          END IF
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idbPW(i),        &
     &                     S(ng)%pioVar(idbPW(i)),                      &
     &                     S(ng)%Rindex,                                &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     SEDBIOL(ng) % sedPoreWaterCon(:,:,i))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idbPW(i))), S(ng)%Rindex,  &
     &                        TRIM(ncname)
		END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out sediment particulate matter decay rates.
!
      DO i=1,NDR
        IF (Lrestart.or.VarOut(idbDR(i),ng)) THEN
          scale=1.0_r8
          IF (S(ng)%pioVar(idbDR(i))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_bs3dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_bs3dvar(ng)
          END IF
          status=nf_fwrite3d(ng, model, S(ng)%pioFile, idbDR(i),        &
     &                       S(ng)%pioVar(idbDR(i)),                    &
     &                       S(ng)%Rindex,                              &
     &                       LBi, UBi, LBj, UBj, 1, nspc, scale,        &
#  ifdef MASKING
     &                       GRID(ng) % rmask,                          &
#  endif
     &                       SEDBIOL(ng) % sedDecayRate(:,:,:,i),       &
     &                       SetFillVal= .FALSE.)
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idbDR(i))), S(ng)%Rindex,  &
     &                          TRIM(ncname)
            END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
!  Write out sediment-water fluxes.
!
      DO i=1,NSF
        IF (Lrestart.or.VarOut(idbSF(i),ng)) THEN
          scale=1.0_r8
	    IF (S(ng)%pioVar(idbSF(i))%dkind.eq.PIO_double) THEN
            ioDesc => ioDesc_dp_r2dvar(ng)
          ELSE
            ioDesc => ioDesc_sp_r2dvar(ng)
          END IF
          status=nf_fwrite2d(ng, model, S(ng)%pioFile, idbSF(i),        &
     &                     S(ng)%pioVar(idbSF(i)),                      &
     &                     S(ng)%Rindex,                                &
     &                     LBi, UBi, LBj, UBj, scale,                   &
#  ifdef MASKING
     &                     GRID(ng) % rmask,                            &
#  endif
     &                     SEDBIOL(ng) % sedFlux(:,:,i))
          IF (FoundError(status, PIO_noerr, __LINE__, MyFile)) THEN
            IF (Master) THEN
              WRITE (stdout,10) TRIM(Vname(1,idbSF(i))), S(ng)%Rindex,  &
     &                        TRIM(ncname)
		END IF
            exit_flag=3
            ioerror=status
            RETURN
          END IF
        END IF
      END DO
!
  10  FORMAT (/," SEDBIO_WRT_PIO - error while writing variable '", &
     &        a,"', time record = ",i0,/,11x,'into file: ',a)
!
      RETURN
      END SUBROUTINE sedbio_wrt_pio
# endif
!
# ifdef STATIONS
!
!***********************************************************************
      SUBROUTINE sedbio_def_station_nf90 (ng, model, ldef, VarOut, S,      &
     &                                 bbgrd, pgrd)
!***********************************************************************
!
      USE mod_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: ldef, VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: pgrd(:), bbgrd(:)
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: got_var(NV)
!
      integer, parameter :: Natt = 25

      integer :: i, j, status, itrc
!
      real(r8) :: Aval(6)
!
      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", sedbio_def_station_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Define sedbio output stations variables.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=S(ng)%name
!
      DEFINE : IF (ldef) THEN
!
!  Initialize local information variable arrays.
!
        DO i=1,Natt
          DO j=1,LEN(Vinfo(1))
            Vinfo(i)(j:j)=' '
          END DO
        END DO
        DO i=1,6
          Aval(i)=0.0_r8
        END DO
!
!  Define sediment particulate matter concentrations.
!
        DO i=1,NPOM
	    IF (VarOut(idbPM(i),ng)) THEN
            Vinfo( 1)=Vname(1,idbPM(i))
		Vinfo( 2)=Vname(2,idbPM(i))
            Vinfo( 3)=Vname(3,idbPM(i))
            Vinfo(14)=Vname(4,idbPM(i))
            Vinfo(16)=Vname(1,idtime)
            status=def_var(ng, model, S(ng)%ncid,                    &
     &                     S(ng)%Vid(idbPM(i)), NF_FOUT,              &
     &                     3, bbgrd, Aval, Vinfo, ncname,            &
     &                     SetFillVal = .FALSE.,                      &
     &                     SetParAccess = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
!  Define sediment particulate matter decay rates.
!
        DO i=1,NDR
          IF (VarOut(idbDR(i),ng)) THEN
            Vinfo( 1)=Vname(1,idbDR(i))
		Vinfo( 2)=Vname(2,idbDR(i))
            Vinfo( 3)=Vname(3,idbDR(i))
            Vinfo(14)=Vname(4,idbDR(i))
            Vinfo(16)=Vname(1,idtime)
            status=def_var(ng, model, S(ng)%ncid,                    &
     &                     S(ng)%Vid(idbDR(i)), NF_FOUT,              &
     &                     3, bbgrd, Aval, Vinfo, ncname,            &
     &                     SetFillVal = .FALSE.,                      &
     &                     SetParAccess = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
!  Define sediment porewater concentrations.
!
        DO i=1,NPWC
          IF (VarOut(idbPW(i),ng)) THEN
            Vinfo( 1)=Vname(1,idbPW(i))
		Vinfo( 2)=Vname(2,idbPW(i))
            Vinfo( 3)=Vname(3,idbPW(i))
            Vinfo(14)=Vname(4,idbPW(i))
            Vinfo(16)=Vname(1,idtime)
            status=def_var(ng, model, S(ng)%ncid,                    &
     &                     S(ng)%Vid(idbPW(i)), NF_FOUT,              &
     &                     2, pgrd, Aval, Vinfo, ncname,             &
     &                     SetFillVal = .FALSE.,                      &
     &                     SetParAccess = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
!  Define sediment-water fluxes.
!
        DO i=1,NSF
          IF (VarOut(idbSF(i),ng)) THEN
            Vinfo( 1)=Vname(1,idbSF(i))
		Vinfo( 2)=Vname(2,idbSF(i))
            Vinfo( 3)=Vname(3,idbSF(i))
            Vinfo(14)=Vname(4,idbSF(i))
            Vinfo(16)=Vname(1,idtime)
            status=def_var(ng, model, S(ng)%ncid,                    &
     &                     S(ng)%Vid(idbSF(i)), NF_FOUT,              &
     &                     2, pgrd, Aval, Vinfo, ncname,             &
     &                     SetFillVal = .FALSE.,                      &
     &                     SetParAccess = .FALSE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
      END IF DEFINE
!
!-----------------------------------------------------------------------
!  Otherwise, check existing output file and prepare for appending
!  data.
!-----------------------------------------------------------------------
!
      QUERY : IF (.not.ldef) THEN
!
!  Initialize locallogical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  Sediment Biology variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            S(ng)%Vid(idtime)=var_id(i)
          END IF
          DO itrc=1,NPOM
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idbPM(itrc)))) THEN
              got_var(idbPM(itrc))=.TRUE.
              S(ng)%Vid(idbPM(itrc))=var_id(i)
            END IF
          END DO
          DO itrc=1,NPWC
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idbPW(itrc)))) THEN
              got_var(idbPW(itrc))=.TRUE.
              S(ng)%Vid(idbPW(itrc))=var_id(i)
            END IF
          END DO
          DO itrc=1,NDR
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idbDR(itrc)))) THEN
              got_var(idbDR(itrc))=.TRUE.
              S(ng)%Vid(idbDR(itrc))=var_id(i)
            END IF
          END DO
          DO itrc=1,NSF
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idbSF(itrc)))) THEN
              got_var(idbSF(itrc))=.TRUE.
              S(ng)%Vid(idbSF(itrc))=var_id(i)
            END IF
          END DO
        END DO
!
!  Check if output variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO i=1,NPOM
          IF (.not.got_var(idbPM(i)).and.VarOut(idbPM(i),ng)) THEN        
            IF (Master) WRITE (stdout,20) TRIM(Vname(1,idbPM(i))),     &
     &                                    TRIM(ncname)
                exit_flag=3
                RETURN
          END IF
        END DO
        DO i=1,NPWC
          IF (.not.got_var(idbPW(i)).and.VarOut(idbPW(i),ng)) THEN
            IF (Master) WRITE (stdout,20) TRIM(Vname(1,idbPW(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO i=1,NDR
          IF (.not.got_var(idbDR(i)).and.VarOut(idbDR(i),ng)) THEN
            IF (Master) WRITE (stdout,20) TRIM(Vname(1,idbDR(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO i=1,NSF
          IF (.not.got_var(idbSF(i)).and.VarOut(idbSF(i),ng)) THEN
            IF (Master) WRITE (stdout,20) TRIM(Vname(1,idbSF(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO

      END IF QUERY 
!
  20  FORMAT (/,' SEDBIO_DEF_STATION_NF90 - unable to find variable:',&
     &        1x,a,2x,' in output NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE sedbio_def_station_nf90
      
!
!***********************************************************************
      SUBROUTINE sedbio_wrt_station_nf90 (ng, model, tile,            &
     &                                      LBi, UBi, LBj, UBj,         &
     &                                      VarOut, S)
!***********************************************************************
!
      USE mod_netcdf
      USE mod_sedbio
!
!  Imported variable declarations.
!
      logical, intent(in) :: VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: Cgrid
!
      integer :: NposB
      integer :: i, k, np, status
!
      real(dp) :: scale
!
      real(r8), dimension(Nstation(ng)) :: Xpos, Ypos, Zpos, psta
      real(r8), dimension(Nstation(ng)*nspc) :: XposSB, YposSB, ZposSB
      real(r8), dimension(Nstation(ng)*nspc) :: sbsta
      integer :: NposSB, indx, kk
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", sedbio_wrt_station_nf90"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out sedbio output variables into specified stations NetCDF
!  file.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set switch to extract station data at native C-grid position (TRUE)
!  or at RHO-points (FALSE).
!
#  ifdef STATIONS_CGRID
      Cgrid=.TRUE.
#  else
      Cgrid=.FALSE.
#  endif
!
!  Set positions for generic extraction routine.
!
      NposSB=Nstation(ng)*nspc
      DO i=1,Nstation(ng)
        Xpos(i)=SCALARS(ng)%SposX(i)
        Ypos(i)=SCALARS(ng)%SposY(i)
        Zpos(i)=1.0_r8
        DO k=1,nspc
          np=k+(i-1)*nspc
          XposSB(np)=SCALARS(ng)%SposX(i)
          YposSB(np)=SCALARS(ng)%SposY(i)
          ZposSB(np)=REAL(k,r8)
        END DO
      END DO
!
!  Write out sediment particulate matter concentrations.
!
      DO i=1,NPOM
        IF (VarOut(idbPM(i),ng)) THEN
          scale=1.0_dp
          CALL extract_sta3d (ng, model, Cgrid, idbPM(i), bs3dvar,      &
     &                        LBi, UBi, LBj, UBj, 1, nspc,              &
     &                        scale, SEDBIOL(ng)%sedPOM(:,:,:,i),       &
     &                        NposSB, XposSB, YposSB, ZposSB, sbsta)
          CALL netcdf_put_fvar (ng, model, S(ng)%name,                  &
     &                          TRIM(Vname(1,idbPM(i))), sbsta,         &
     &                          (/1,1,S(ng)%Rindex/),                   &
     &                          (/nspc,Nstation(ng),1/),                &
     &                          ncid = S(ng)%ncid,                      &
     &                          varid = S(ng)%Vid(idbPM(i)))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO 
!
!  Write out sediment particulate matter decay rates.
!
      DO i=1,NDR
        IF (VarOut(idbDR(i),ng)) THEN
          scale=1.0_dp
          CALL extract_sta3d (ng, model, Cgrid, idbDR(i), bs3dvar,      &
     &                        LBi, UBi, LBj, UBj, 1, nspc,              &
     &                        scale, SEDBIOL(ng)%sedDecayRate(:,:,:,i), &
     &                        NposSB, XposSB, YposSB, ZposSB, sbsta)
          CALL netcdf_put_fvar (ng, model, S(ng)%name,                  &
     &                          TRIM(Vname(1,idbDR(i))), sbsta,         &
     &                          (/1,1,S(ng)%Rindex/),                   &
     &                          (/nspc,Nstation(ng),1/),                &
     &                          ncid = S(ng)%ncid,                      &
     &                          varid = S(ng)%Vid(idbDR(i)))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
!
!  Write out sediment porewater concentrations.
!
      DO i=1,NPWC
        IF (VarOut(idbPW(i),ng)) THEN
          scale=1.0_dp
          CALL extract_sta2d (ng, model, Cgrid, idbPW(i), r2dvar,       &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        scale, SEDBIOL(ng)%sedPoreWaterCon(:,:,i),&
     &                        Nstation(ng), Xpos, Ypos, psta)
          CALL netcdf_put_fvar (ng, model, S(ng)%name,                  &
     &                          TRIM(Vname(1,idbPW(i))), psta,          &
     &                          (/1,S(ng)%Rindex/),(/Nstation(ng),1/),  &
     &                          ncid = S(ng)%ncid,                    &
     &                          varid = S(ng)%Vid(idbPW(i)))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
!
!  Write out sediment-water fluxes
!
      DO i=1,NSF
        IF (VarOut(idbSF(i),ng)) THEN
          scale=1.0_dp
          CALL extract_sta2d (ng, model, Cgrid, idbSF(i), r2dvar,       &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        scale, SEDBIOL(ng)%sedFlux(:,:,i),        &
     &                        Nstation(ng), Xpos, Ypos, psta)
          CALL netcdf_put_fvar (ng, model, S(ng)%name,                  &
     &                          TRIM(Vname(1,idbSF(i))), psta,          &
     &                          (/1,S(ng)%Rindex/),(/Nstation(ng),1/),  &
     &                          ncid = S(ng)%ncid,                    &
     &                          varid = S(ng)%Vid(idbSF(i)))
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
      
      RETURN
      END SUBROUTINE sedbio_wrt_station_nf90
      
#  if defined PIO_LIB && defined DISTRIBUTE
!
!***********************************************************************
      SUBROUTINE sedbio_def_station_pio (ng, model, ldef, VarOut, S,      &
     &                                 bbgrd, pgrd)
!***********************************************************************
!
      USE mod_pio_netcdf
!
!  Imported variable declarations.
!
      logical, intent(in) :: ldef, VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model
      integer, intent(in), optional :: pgrd(:), bbgrd(:)
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: got_var(NV)
!
      integer, parameter :: Natt = 25

      integer :: i, j, status, itrc
!
      real(r8) :: Aval(6)
!
      character (len=120) :: Vinfo(Natt)
      character (len=256) :: ncname
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", sedbio_def_station_pio"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Define sedbio output stations variables.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      ncname=S(ng)%name
!
      DEFINE : IF (ldef) THEN
!
!  Initialize local information variable arrays.
!
        DO i=1,Natt
          DO j=1,LEN(Vinfo(1))
            Vinfo(i)(j:j)=' '
          END DO
        END DO
        DO i=1,6
          Aval(i)=0.0_r8
        END DO
!
!  Define sediment particulate matter concentrations.
!
        DO i=1,NPOM
	    IF (VarOut(idbPM(i),ng)) THEN
            Vinfo( 1)=Vname(1,idbPM(i))
		Vinfo( 2)=Vname(2,idbPM(i))
            Vinfo( 3)=Vname(3,idbPM(i))
            Vinfo(14)=Vname(4,idbPM(i))
            Vinfo(16)=Vname(1,idtime)
            S(ng)%pioVar(idbPM(i))%dkind=PIO_FOUT
            S(ng)%pioVar(idbPM(i))%gtype=0
            status=def_var(ng, model, S(ng)%pioFile,                 &
     &                     S(ng)%pioVar(idbPM(i))%vd,                &
     &                     PIO_FOUT, 3, bbgrd, Aval, Vinfo, ncname,  &
     &                     SetFillVal = .TRUE.,                     &
     &                     SetParAccess = .TRUE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
!  Define sediment particulate matter decay rates.
!
        DO i=1,NDR
          IF (VarOut(idbDR(i),ng)) THEN
            Vinfo( 1)=Vname(1,idbDR(i))
		Vinfo( 2)=Vname(2,idbDR(i))
            Vinfo( 3)=Vname(3,idbDR(i))
            Vinfo(14)=Vname(4,idbDR(i))
            Vinfo(16)=Vname(1,idtime)
            S(ng)%pioVar(idbDR(i))%dkind=PIO_FOUT
            S(ng)%pioVar(idbDR(i))%gtype=0
            status=def_var(ng, model, S(ng)%pioFile,                 &
     &                     S(ng)%pioVar(idbDR(i))%vd,                &
     &                     PIO_FOUT, 3, bbgrd, Aval, Vinfo, ncname,  &
     &                     SetFillVal = .TRUE.,                     &
     &                     SetParAccess = .TRUE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
!  Define sediment porewater concentrations.
!
        DO i=1,NPWC
          IF (VarOut(idbPW(i),ng)) THEN
            Vinfo( 1)=Vname(1,idbPW(i))
		Vinfo( 2)=Vname(2,idbPW(i))
            Vinfo( 3)=Vname(3,idbPW(i))
            Vinfo(14)=Vname(4,idbPW(i))
            Vinfo(16)=Vname(1,idtime)
            S(ng)%pioVar(idbPW(i))%dkind=PIO_FOUT
            S(ng)%pioVar(idbPW(i))%gtype=0
            status=def_var(ng, model, S(ng)%pioFile,                 &
     &                     S(ng)%pioVar(idbPW(i))%vd,                &
     &                     PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,   &
     &                     SetFillVal = .TRUE.,                     &
     &                     SetParAccess = .TRUE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
!
!  Define sediment-water fluxes.
!
        DO i=1,NSF
          IF (VarOut(idbSF(i),ng)) THEN
            Vinfo( 1)=Vname(1,idbSF(i))
		Vinfo( 2)=Vname(2,idbSF(i))
            Vinfo( 3)=Vname(3,idbSF(i))
            Vinfo(14)=Vname(4,idbSF(i))
            Vinfo(16)=Vname(1,idtime)
            S(ng)%pioVar(idbSF(i))%dkind=PIO_FOUT
            S(ng)%pioVar(idbSF(i))%gtype=0
            status=def_var(ng, model, S(ng)%pioFile,                 &
     &                     S(ng)%pioVar(idbSF(i))%vd,                &
     &                     PIO_FOUT, 2, pgrd, Aval, Vinfo, ncname,   &
     &                     SetFillVal = .TRUE.,                     &
     &                     SetParAccess = .TRUE.)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
        END DO
      END IF DEFINE
!
!-----------------------------------------------------------------------
!  Otherwise, check existing output file and prepare for appending
!  data.
!-----------------------------------------------------------------------
!
      QUERY : IF (.not.ldef) THEN
!
!  Initialize locallogical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from input NetCDF and activate switches for
!  Sediment Biology variables. Get variable IDs.
!
        DO i=1,n_var
          IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idtime))) THEN
            got_var(idtime)=.TRUE.
            S(ng)%pioVar(idtime)%vd=var_desc(i)
            S(ng)%pioVar(idtime)%dkind=PIO_TOUT
            S(ng)%pioVar(idtime)%gtype=0
          END IF
          DO itrc=1,NPOM
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idbPM(itrc)))) THEN
              got_var(idbPM(itrc))=.TRUE.
              S(ng)%pioVar(idbPM(itrc))%vd=var_desc(i)
              S(ng)%pioVar(idbPM(itrc))%dkind=PIO_FOUT
              S(ng)%pioVar(idbPM(itrc))%gtype=0
            END IF
          END DO
          DO itrc=1,NPWC
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idbPW(itrc)))) THEN
              got_var(idbPW(itrc))=.TRUE.
              S(ng)%pioVar(idbPW(itrc))%vd=var_desc(i)
              S(ng)%pioVar(idbPW(itrc))%dkind=PIO_FOUT
              S(ng)%pioVar(idbPW(itrc))%gtype=0
            END IF
          END DO
          DO itrc=1,NDR
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idbDR(itrc)))) THEN
              got_var(idbDR(itrc))=.TRUE.
              S(ng)%pioVar(idbDR(itrc))%vd=var_desc(i)
              S(ng)%pioVar(idbDR(itrc))%dkind=PIO_FOUT
              S(ng)%pioVar(idbDR(itrc))%gtype=0
            END IF
          END DO
          DO itrc=1,NSF
            IF (TRIM(var_name(i)).eq.TRIM(Vname(1,idbSF(itrc)))) THEN
              got_var(idbSF(itrc))=.TRUE.
              S(ng)%pioVar(idbSF(itrc))%vd=var_desc(i)
              S(ng)%pioVar(idbSF(itrc))%dkind=PIO_FOUT
              S(ng)%pioVar(idbSF(itrc))%gtype=0
            END IF
          END DO
        END DO
!
!  Check if output variables are available in input NetCDF file.
!
        IF (.not.got_var(idtime)) THEN
          IF (Master) WRITE (stdout,20) TRIM(Vname(1,idtime)),          &
     &                                  TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
        DO i=1,NPOM
          IF (.not.got_var(idbPM(i)).and.VarOut(idbPM(i),ng))) THEN        
            IF (Master) WRITE (stdout,20) TRIM(Vname(1,idbPM(i))),     &
     &                                    TRIM(ncname)
                exit_flag=3
                RETURN
          END IF
        END DO
        DO i=1,NPWC
          IF (.not.got_var(idbPW(i)).and.VarOut(idbPW(i),ng))) THEN
            IF (Master) WRITE (stdout,20) TRIM(Vname(1,idbPW(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO i=1,NDR
          IF (.not.got_var(idbDR(i)).and.VarOut(idbDR(i),ng))) THEN
            IF (Master) WRITE (stdout,20) TRIM(Vname(1,idbDR(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO
        DO i=1,NSF
          IF (.not.got_var(idbSF(i)).and.VarOut(idbSF(i),ng))) THEN
            IF (Master) WRITE (stdout,20) TRIM(Vname(1,idbSF(i))),     &
     &                                    TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END DO

      END IF QUERY 
!
  20  FORMAT (/,' SEDBIO_DEF_STATION_PIO - unable to find variable:',&
     &        1x,a,2x,' in output NetCDF file: ',a)
!
      RETURN
      END SUBROUTINE sedbio_def_station_pio
!
!***********************************************************************
      SUBROUTINE sedbio_wrt_station_pio (ng, model, tile,               &
     &                                      LBi, UBi, LBj, UBj,         &
     &                                      VarOut, S)
!***********************************************************************
!
      USE mod_pio_netcdf
      USE mod_sedbio
!
!  Imported variable declarations.
!
      logical, intent(in) :: VarOut(NV,Ngrids)
!
      integer, intent(in) :: ng, model, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
!
      TYPE(T_IO), intent(inout) :: S(Ngrids)
!
!  Local variable declarations.
!
      logical :: Cgrid
!
      integer :: NposB
      integer :: i, k, np, status
!
      real(dp) :: scale
!
      real(r8), dimension(Nstation(ng)) :: Xpos, Ypos, Zpos, psta
      real(r8), dimension(Nstation(ng)*nspc) :: XposSB, YposSB, ZposSB
      real(r8), dimension(Nstation(ng)*nspc) :: sbsta
      integer :: NposSB, indx, kk
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", sedbio_wrt_station_pio"
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Write out sedbio output variables into specified stations NetCDF
!  file.
!-----------------------------------------------------------------------
!
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set switch to extract station data at native C-grid position (TRUE)
!  or at RHO-points (FALSE).
!
#  ifdef STATIONS_CGRID
      Cgrid=.TRUE.
#  else
      Cgrid=.FALSE.
#  endif
!
!  Set positions for generic extraction routine.
!
      NposSB=Nstation(ng)*nspc
      DO i=1,Nstation(ng)
        Xpos(i)=SCALARS(ng)%SposX(i)
        Ypos(i)=SCALARS(ng)%SposY(i)
        Zpos(i)=1.0_r8
        DO k=1,nspc
          np=k+(i-1)*nspc
          XposSB(np)=SCALARS(ng)%SposX(i)
          YposSB(np)=SCALARS(ng)%SposY(i)
          ZposSB(np)=REAL(k,r8)
        END DO
      END DO
!
!  Write out sediment particulate matter concentrations.
!
      DO i=1,NPOM
        IF (VarOut(idbPM(i),ng)) THEN
          scale=1.0_dp
          CALL extract_sta3d (ng, model, Cgrid, idbPM(i), bs3dvar,       &
     &                        LBi, UBi, LBj, UBj, 1, nspc,              &
     &                        scale, SEDBIOL(ng)%sedPOM(:,:,:,i),       &
     &                        NposSB, XposSB, YposSB, ZposSB, sbsta)
          CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,              &
     &                          TRIM(Vname(1,idbPM(i))), sbsta,         &
     &                          (/1,1,S(ng)%Rindex/),                   &
     &                          (/nspc,Nstation(ng),1/),                &
     &                          pioFile = S(ng)%pioFile,                &
     &                          pioVar = S(ng)%pioVar(idbPM(i))%vd)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO 
!
!  Write out sediment particulate matter decay rates.
!
      DO i=1,NDR
        IF (VarOut(idbDR(i),ng)) THEN
          scale=1.0_dp
          CALL extract_sta3d (ng, model, Cgrid, idbDR(i), bs3dvar,       &
     &                        LBi, UBi, LBj, UBj, 1, nspc,              &
     &                        scale, SEDBIOL(ng)%sedDecayRate(:,:,:,i), &
     &                        NposSB, XposSB, YposSB, ZposSB, sbsta)
          CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,              &
     &                          TRIM(Vname(1,idbDR(i))), sbsta,         &
     &                          (/1,1,S(ng)%Rindex/),                   &
     &                          (/nspc,Nstation(ng),1/),                &
     &                          pioFile = S(ng)%pioFile,                &
     &                          pioVar = S(ng)%pioVar(idbDR(i))%vd)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
!
!  Write out sediment porewater concentrations.
!
      DO i=1,NPWC
        IF (VarOut(idbPW(i),ng)) THEN
          scale=1.0_dp
          CALL extract_sta2d (ng, model, Cgrid, idbPW(i), r2dvar,       &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        scale, SEDBIOL(ng)%sedPoreWaterCon(:,:,i),&
     &                        Nstation(ng), Xpos, Ypos, psta)
          CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,              &
     &                          TRIM(Vname(1,idbPW(i))), psta,          &
     &                          (/1,S(ng)%Rindex/),(/Nstation(ng),1/),  &
     &                          pioFile = S(ng)%pioFile,                &
     &                          pioVar = S(ng)%pioVar(idbPW(i))%vd)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
!
!  Write out sediment-water fluxes
!
      DO i=1,NSF
        IF (VarOut(idbSF(i),ng)) THEN
          scale=1.0_dp
          CALL extract_sta2d (ng, model, Cgrid, idbSF(i), r2dvar,       &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        scale, SEDBIOL(ng)%sedFlux(:,:,i),        &
     &                        Nstation(ng), Xpos, Ypos, psta)
          CALL pio_netcdf_put_fvar (ng, model, S(ng)%name,              &
     &                          TRIM(Vname(1,idbSF(i))), psta,          &
     &                          (/1,S(ng)%Rindex/),(/Nstation(ng),1/),  &
     &                          pioFile = S(ng)%pioFile,                &
     &                          pioVar = S(ng)%pioVar(idbSF(i))%vd)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
      END DO
      
      RETURN
      END SUBROUTINE sedbio_wrt_station_pio
#  endif
# endif
#endif
      END MODULE sedbio_output_mod
