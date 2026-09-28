MODULE ionic_optim_utils
  ! The wavefunction behind sess_* stays the one INTERFACE allocated.
  ! rgsaddle and gpr_optim are loaded by name. Their source is not
  ! part of this program. The contract is include/opencpmd_ionic.h.
  USE coor,                            ONLY: fion,&
                                             tau0,&
                                             taup,&
                                             velp
  USE copot_utils,                     ONLY: copot
  USE elct,                            ONLY: crge
  USE ener,                            ONLY: ener_com
  USE error_handling,                  ONLY: stopgm
  USE forcedr_driver,                  ONLY: forcedr
  USE ions,                            ONLY: ions0,&
                                             ions1
  USE iso_c_binding,                   ONLY: c_char,&
                                             c_double,&
                                             c_f_pointer,&
                                             c_f_procpointer,&
                                             c_funloc,&
                                             c_funptr,&
                                             c_int,&
                                             c_int32_t,&
                                             c_int64_t,&
                                             c_loc,&
                                             c_null_char,&
                                             c_null_funptr,&
                                             c_null_ptr,&
                                             c_ptr,&
                                             c_associated
  USE kinds,                           ONLY: real_8
  USE machine,                         ONLY: m_walltime
  USE metr,                            ONLY: metr_com
  USE mp_interface,                    ONLY: mp_bcast,&
                                             mp_sync
  USE nlcc,                            ONLY: corel
  USE norm,                            ONLY: cnorm,&
                                             gemax
  USE parac,                           ONLY: parai,&
                                             paral
  USE phfac_utils,                     ONLY: phfac
  USE purge_utils,                     ONLY: purge
  USE rrane_utils,                     ONLY: rrane
  USE setirec_utils,                   ONLY: write_irec
  USE system,                          ONLY: cnti,&
                                             cntl,&
                                             cntr,&
                                             maxsys,&
                                             nkpt
  USE updwf_utils,                     ONLY: updwf
  USE wv30_utils,                       ONLY: zhwwf
  USE wrener_utils,                    ONLY: wrener

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: ionic_bind, ionic_rgsaddle, ionic_gpr

  COMPLEX(real_8), POINTER, SAVE :: sess_c0(:,:) => NULL()
  COMPLEX(real_8), POINTER, SAVE :: sess_c2(:,:) => NULL()
  COMPLEX(real_8), POINTER, SAVE :: sess_sc0(:,:) => NULL()
  COMPLEX(real_8), POINTER, SAVE :: sess_pme(:) => NULL()
  COMPLEX(real_8), POINTER, SAVE :: sess_gde(:) => NULL()
  REAL(real_8), POINTER, SAVE :: sess_vpp(:) => NULL()
  REAL(real_8), POINTER, SAVE :: sess_eigv(:) => NULL()
  REAL(real_8), POINTER, SAVE :: sess_rhoe(:,:) => NULL()
  COMPLEX(real_8), POINTER, SAVE :: sess_psi(:,:) => NULL()

  INTEGER, PARAMETER :: rtld_lazy = 1

  TYPE, BIND(C) :: rs_version
     INTEGER(c_int32_t) :: major
     INTEGER(c_int32_t) :: minor
  END TYPE rs_version

  TYPE, BIND(C) :: rs_request
     TYPE(rs_version) :: version
     INTEGER(c_int64_t) :: flags
     INTEGER(c_int64_t) :: n_images
     INTEGER(c_int64_t) :: n_atoms
     TYPE(c_ptr) :: positions
     TYPE(c_ptr) :: energies
     TYPE(c_ptr) :: gradients
  END TYPE rs_request

  TYPE, BIND(C) :: rs_band_config
     TYPE(rs_version) :: version
     INTEGER(c_int64_t) :: flags
     INTEGER(c_int32_t) :: tangent
     INTEGER(c_int32_t) :: spring
     INTEGER(c_int32_t) :: projection
     INTEGER(c_int32_t) :: method
     REAL(c_double) :: spring_k
     TYPE(c_ptr) :: spring_ks
     REAL(c_double) :: ci_trigger_factor
     REAL(c_double) :: ci_trigger_force
     TYPE(c_ptr) :: cell
     REAL(c_double) :: force_tol
     REAL(c_double) :: max_move
     INTEGER(c_int64_t) :: memory
  END TYPE rs_band_config

  TYPE, BIND(C) :: rs_min_config
     TYPE(rs_version) :: version
     INTEGER(c_int64_t) :: flags
     INTEGER(c_int32_t) :: kind
     INTEGER(c_int32_t) :: method
     REAL(c_double) :: dr
     REAL(c_double) :: rotation_tol
     INTEGER(c_int64_t) :: max_rotations
     INTEGER(c_int64_t) :: krylov_dim
     REAL(c_double) :: force_tol
     REAL(c_double) :: max_move
  END TYPE rs_min_config

  TYPE, BIND(C) :: rs_report
     TYPE(rs_version) :: version
     INTEGER(c_int64_t) :: flags
     INTEGER(c_int32_t) :: status
     INTEGER(c_int32_t) :: reserved
     REAL(c_double) :: max_force
     INTEGER(c_int64_t) :: ci_index
     INTEGER(c_int64_t) :: iteration
     REAL(c_double) :: curvature
     INTEGER(c_int64_t) :: rotations
  END TYPE rs_report

  TYPE, BIND(C) :: ionic_task
     INTEGER(c_int32_t) :: major
     INTEGER(c_int32_t) :: minor
     INTEGER(c_int64_t) :: n_atoms
     TYPE(c_ptr) :: positions
     TYPE(c_ptr) :: atomic_numbers
     TYPE(c_ptr) :: cell
     TYPE(c_funptr) :: eval
  END TYPE ionic_task

  ABSTRACT INTERFACE
     FUNCTION rs_band_create_i(config, n_images, n_atoms, positions) &
          BIND(C) RESULT(band)
       IMPORT :: c_ptr, c_int64_t
       TYPE(c_ptr), VALUE :: config
       INTEGER(c_int64_t), VALUE :: n_images
       INTEGER(c_int64_t), VALUE :: n_atoms
       TYPE(c_ptr), VALUE :: positions
       TYPE(c_ptr) :: band
     END FUNCTION rs_band_create_i
     FUNCTION rs_band_step_i(band, surface, user, report) BIND(C) RESULT(rc)
       IMPORT :: c_ptr, c_funptr, c_int
       TYPE(c_ptr), VALUE :: band
       TYPE(c_funptr), VALUE :: surface
       TYPE(c_ptr), VALUE :: user
       TYPE(c_ptr), VALUE :: report
       INTEGER(c_int) :: rc
     END FUNCTION rs_band_step_i
     FUNCTION rs_band_pos_i(band, out) BIND(C) RESULT(rc)
       IMPORT :: c_ptr, c_int
       TYPE(c_ptr), VALUE :: band
       TYPE(c_ptr), VALUE :: out
       INTEGER(c_int) :: rc
     END FUNCTION rs_band_pos_i
     SUBROUTINE rs_band_free_i(band) BIND(C)
       IMPORT :: c_ptr
       TYPE(c_ptr), VALUE :: band
     END SUBROUTINE rs_band_free_i
     FUNCTION rs_min_create_i(config, n_atoms, position, mode) &
          BIND(C) RESULT(session)
       IMPORT :: c_ptr, c_int64_t
       TYPE(c_ptr), VALUE :: config
       INTEGER(c_int64_t), VALUE :: n_atoms
       TYPE(c_ptr), VALUE :: position
       TYPE(c_ptr), VALUE :: mode
       TYPE(c_ptr) :: session
     END FUNCTION rs_min_create_i
     FUNCTION rs_min_step_i(session, surface, user, report) BIND(C) RESULT(rc)
       IMPORT :: c_ptr, c_funptr, c_int
       TYPE(c_ptr), VALUE :: session
       TYPE(c_funptr), VALUE :: surface
       TYPE(c_ptr), VALUE :: user
       TYPE(c_ptr), VALUE :: report
       INTEGER(c_int) :: rc
     END FUNCTION rs_min_step_i
     FUNCTION rs_min_pos_i(session, out) BIND(C) RESULT(rc)
       IMPORT :: c_ptr, c_int
       TYPE(c_ptr), VALUE :: session
       TYPE(c_ptr), VALUE :: out
       INTEGER(c_int) :: rc
     END FUNCTION rs_min_pos_i
     SUBROUTINE rs_min_free_i(session) BIND(C)
       IMPORT :: c_ptr
       TYPE(c_ptr), VALUE :: session
     END SUBROUTINE rs_min_free_i
     FUNCTION gpr_run_i(task) BIND(C) RESULT(rc)
       IMPORT :: c_ptr, c_int
       TYPE(c_ptr), VALUE :: task
       INTEGER(c_int) :: rc
     END FUNCTION gpr_run_i
     FUNCTION dlopen_i(name, flags) BIND(C, name='dlopen') RESULT(handle)
       IMPORT :: c_ptr, c_char, c_int
       CHARACTER(kind=c_char), INTENT(in) :: name(*)
       INTEGER(c_int), VALUE :: flags
       TYPE(c_ptr) :: handle
     END FUNCTION dlopen_i
     FUNCTION dlsym_i(handle, name) BIND(C, name='dlsym') RESULT(sym)
       IMPORT :: c_ptr, c_char
       TYPE(c_ptr), VALUE :: handle
       CHARACTER(kind=c_char), INTENT(in) :: name(*)
       TYPE(c_ptr) :: sym
     END FUNCTION dlsym_i
     FUNCTION dlerror_i() BIND(C, name='dlerror') RESULT(msg)
       IMPORT :: c_ptr
       TYPE(c_ptr) :: msg
     END FUNCTION dlerror_i
     FUNCTION dlclose_i(handle) BIND(C, name='dlclose') RESULT(rc)
       IMPORT :: c_ptr, c_int
       TYPE(c_ptr), VALUE :: handle
       INTEGER(c_int) :: rc
     END FUNCTION dlclose_i
  END INTERFACE

  REAL(c_double), ALLOCATABLE, TARGET, SAVE :: rs_pos(:)
  REAL(c_double), TARGET, SAVE :: deck_cell(9)
  INTEGER(c_int32_t), ALLOCATABLE, TARGET, SAVE :: deck_z(:)

CONTAINS

  SUBROUTINE ionic_bind(c0, c2, sc0, pme, gde, vpp, eigv, rhoe, psi)
    COMPLEX(real_8), TARGET, INTENT(inout) :: c0(:,:)
    COMPLEX(real_8), TARGET, INTENT(inout) :: c2(:,:)
    COMPLEX(real_8), TARGET, INTENT(inout) :: sc0(:,:)
    COMPLEX(real_8), TARGET, INTENT(inout) :: pme(:)
    COMPLEX(real_8), TARGET, INTENT(inout) :: gde(:)
    REAL(real_8), TARGET, INTENT(inout) :: vpp(:)
    REAL(real_8), TARGET, INTENT(inout) :: eigv(:)
    REAL(real_8), TARGET, INTENT(inout) :: rhoe(:,:)
    COMPLEX(real_8), TARGET, INTENT(inout) :: psi(:,:)
    sess_c0 => c0
    sess_c2 => c2
    sess_sc0 => sc0
    sess_pme => pme
    sess_gde => gde
    sess_vpp => vpp
    sess_eigv => eigv
    sess_rhoe => rhoe
    sess_psi => psi
  END SUBROUTINE ionic_bind

  INTEGER FUNCTION nat_deck()
    INTEGER :: is
    nat_deck = 0
    DO is = 1, ions1%nsp
       nat_deck = nat_deck + ions0%na(is)
    END DO
  END FUNCTION nat_deck

  SUBROUTINE flat_from_tau(pos)
    REAL(c_double), INTENT(out) :: pos(:)
    INTEGER :: is, ia, k
    k = 0
    DO is = 1, ions1%nsp
       DO ia = 1, ions0%na(is)
          pos(k+1) = tau0(1, ia, is)
          pos(k+2) = tau0(2, ia, is)
          pos(k+3) = tau0(3, ia, is)
          k = k + 3
       END DO
    END DO
  END SUBROUTINE flat_from_tau

  SUBROUTINE tau_from_flat(pos)
    REAL(c_double), INTENT(in) :: pos(:)
    INTEGER :: is, ia, k
    k = 0
    DO is = 1, ions1%nsp
       DO ia = 1, ions0%na(is)
          tau0(1, ia, is) = pos(k+1)
          tau0(2, ia, is) = pos(k+2)
          tau0(3, ia, is) = pos(k+3)
          k = k + 3
       END DO
    END DO
  END SUBROUTINE tau_from_flat

  SUBROUTINE fill_species_and_cell()
    INTEGER :: is, ia, k, nat
    nat = nat_deck()
    IF (ALLOCATED(deck_z)) DEALLOCATE(deck_z)
    ALLOCATE(deck_z(nat))
    k = 0
    DO is = 1, ions1%nsp
       DO ia = 1, ions0%na(is)
          k = k + 1
          deck_z(k) = ions0%iatyp(is)
       END DO
    END DO
    deck_cell = 0.0_c_double
    deck_cell(1) = metr_com%ht(1, 1)
    deck_cell(5) = metr_com%ht(2, 2)
    deck_cell(9) = metr_com%ht(3, 3)
  END SUBROUTINE fill_species_and_cell

  SUBROUTINE ionic_one_point()
    ! One geometry already stored in tau0. c0 is the previous point.
    CHARACTER(*), PARAMETER :: procedureN = 'ionic_one_point'
    INTEGER :: infi, irec(100)
    LOGICAL :: fnowf
    REAL(real_8) :: detot, dummy(1), etoto, tcpu, time1, time2
    CALL ropt_prepare()
    CALL phfac(tau0)
    IF (corel%tinlc) CALL copot(sess_rhoe, sess_psi, .FALSE.)
    CALL mp_sync(parai%allgrp)
    IF (cntl%trane) THEN
       CALL rrane(sess_c0, sess_c2, crge%n)
    END IF
    fnowf = .FALSE.
    etoto = 0.0_real_8
    DO infi = 1, cnti%nomore
       time1 = m_walltime()
       IF (paral%parent) THEN
          ropt_mod_engpri(infi)
       END IF
       CALL updwf(sess_c0, sess_c2, sess_sc0, tau0, fion, sess_pme, &
            sess_gde, sess_vpp, sess_eigv, sess_rhoe, sess_psi, &
            crge%n, fnowf, .TRUE.)
       IF (paral%parent) THEN
          detot = ener_com%etot + ener_com%eext - etoto
          IF (infi .EQ. 1) detot = 0.0_real_8
          time2 = m_walltime()
          tcpu = (time2 - time1) * 0.001_real_8
          IF (paral%io_parent) THEN
             CALL wrener
             WRITE(6, '(I4,F10.6,F10.6,F14.5,F14.6,F12.2)') &
                  infi, gemax, cnorm, ener_com%etot + ener_com%eext, detot, tcpu
          END IF
          etoto = ener_com%etot + ener_com%eext
       END IF
       IF (ropt_conv() .AND. fnowf) EXIT
       IF (ropt_conv()) THEN
          CALL ropt_clear_conv()
          fnowf = .TRUE.
       END IF
    END DO
    CALL forcedr(sess_c0, sess_c2, sess_sc0, sess_rhoe, sess_psi, tau0, &
         fion, sess_eigv, crge%n, 1, .TRUE., .TRUE.)
    CALL purge(tau0, fion)
    CALL write_irec(irec)
    dummy(1) = 0.0_real_8
    CALL zhwwf(2, irec, sess_c0, sess_c2, crge%n, dummy, tau0, velp, &
         taup, iteropt_nfi())
  END SUBROUTINE ionic_one_point

  ! ropt_mod and iteropt live in the same modules egointer uses.
  ! Wrappers keep this file from depending on a host association.

  SUBROUTINE ropt_prepare()
    USE ropt,                          ONLY: ropt_mod
    ropt_mod%convwf = .FALSE.
    ropt_mod%sdiis = .TRUE.
    ropt_mod%spcg = .TRUE.
    ropt_mod%modens = .FALSE.
    ropt_mod%engpri = .FALSE.
    ropt_mod%calste = .FALSE.
  END SUBROUTINE ropt_prepare

  SUBROUTINE ropt_mod_engpri(infi)
    USE ropt,                          ONLY: ropt_mod
    USE store_types,                   ONLY: cprint
    INTEGER, INTENT(in) :: infi
    ropt_mod%engpri = MOD(infi - 1, cprint%iprint_step) .EQ. 0
    IF (.NOT. paral%parent) ropt_mod%engpri = .FALSE.
  END SUBROUTINE ropt_mod_engpri

  LOGICAL FUNCTION ropt_conv()
    USE ropt,                          ONLY: ropt_mod
    ropt_conv = ropt_mod%convwf
  END FUNCTION ropt_conv

  SUBROUTINE ropt_clear_conv()
    USE ropt,                          ONLY: ropt_mod
    ropt_mod%convwf = .FALSE.
  END SUBROUTINE ropt_clear_conv

  INTEGER FUNCTION iteropt_nfi()
    USE ropt,                          ONLY: iteropt
    iteropt%nfi = iteropt%nfi + 1
    iteropt_nfi = iteropt%nfi
  END FUNCTION iteropt_nfi

  FUNCTION opencpmd_ionic_eval(n_atoms, positions, energy, forces) &
       BIND(C, name='opencpmd_ionic_eval') RESULT(rc)
    INTEGER(c_int64_t), VALUE :: n_atoms
    TYPE(c_ptr), VALUE :: positions
    TYPE(c_ptr), VALUE :: energy
    TYPE(c_ptr), VALUE :: forces
    INTEGER(c_int) :: rc
    REAL(c_double), POINTER :: pos(:)
    REAL(c_double), POINTER :: ene
    REAL(c_double), POINTER :: frc(:)
    INTEGER :: nat, is, ia, k
    rc = 0
    nat = nat_deck()
    IF (n_atoms /= nat .OR. .NOT. c_associated(positions) &
         .OR. .NOT. c_associated(energy) .OR. .NOT. c_associated(forces)) THEN
       rc = -1
       RETURN
    END IF
    CALL c_f_pointer(positions, pos, [3 * nat])
    CALL c_f_pointer(energy, ene)
    CALL c_f_pointer(forces, frc, [3 * nat])
    CALL tau_from_flat(pos)
    CALL ionic_one_point()
    ene = ener_com%etot + ener_com%eext
    k = 0
    DO is = 1, ions1%nsp
       DO ia = 1, ions0%na(is)
          frc(k+1) = fion(1, ia, is)
          frc(k+2) = fion(2, ia, is)
          frc(k+3) = fion(3, ia, is)
          k = k + 3
       END DO
    END DO
  END FUNCTION opencpmd_ionic_eval

  FUNCTION rs_surface(user, reqp) BIND(C) RESULT(rc)
    TYPE(c_ptr), VALUE :: user
    TYPE(c_ptr), VALUE :: reqp
    INTEGER(c_int) :: rc
    TYPE(rs_request), POINTER :: req
    REAL(c_double), POINTER :: pos(:)
    REAL(c_double), POINTER :: ene(:)
    REAL(c_double), POINTER :: grad(:)
    INTEGER :: img, nat, nimg, off, j
    REAL(c_double), ALLOCATABLE :: one_f(:)
    TYPE(c_ptr) :: eptr, fptr
    rc = 0
    IF (.NOT. c_associated(reqp)) THEN
       rc = -1
       RETURN
    END IF
    CALL c_f_pointer(reqp, req)
    nat = nat_deck()
    nimg = INT(req%n_images)
    IF (req%n_atoms /= nat .OR. nimg < 1) THEN
       rc = -1
       RETURN
    END IF
    CALL c_f_pointer(req%positions, pos, [3 * nat * nimg])
    CALL c_f_pointer(req%energies, ene, [nimg])
    CALL c_f_pointer(req%gradients, grad, [3 * nat * nimg])
    ALLOCATE(one_f(3 * nat))
    DO img = 1, nimg
       off = (img - 1) * 3 * nat
       eptr = c_loc(ene(img))
       fptr = c_loc(one_f(1))
       rc = opencpmd_ionic_eval(INT(nat, c_int64_t), c_loc(pos(off + 1)), &
            eptr, fptr)
       IF (rc /= 0) RETURN
       DO j = 1, 3 * nat
          grad(off + j) = -one_f(j)
       END DO
    END DO
    IF (c_associated(user)) CONTINUE
  END FUNCTION rs_surface

  SUBROUTINE lib_name(env, fallback, out)
    CHARACTER(*), INTENT(in) :: env
    CHARACTER(*), INTENT(in) :: fallback
    CHARACTER(len=512), INTENT(out) :: out
    INTEGER :: st
    CALL get_environment_variable(env, out, status=st)
    IF (st /= 0) out = fallback
  END SUBROUTINE lib_name

  SUBROUTINE c_message(ptr, text)
    TYPE(c_ptr), INTENT(in) :: ptr
    CHARACTER(len=256), INTENT(out) :: text
    CHARACTER(kind=c_char), POINTER :: raw(:)
    INTEGER :: i
    text = ''
    IF (.NOT. c_associated(ptr)) RETURN
    CALL c_f_pointer(ptr, raw, [256])
    DO i = 1, 256
       IF (raw(i) == c_null_char) EXIT
       text(i:i) = raw(i)
    END DO
  END SUBROUTINE c_message

  FUNCTION load_sym(handle, name) RESULT(sym)
    TYPE(c_ptr), INTENT(in) :: handle
    CHARACTER(*), INTENT(in) :: name
    TYPE(c_ptr) :: sym
    CHARACTER(len=256) :: err
    sym = dlsym_i(handle, TRIM(name) // c_null_char)
    IF (.NOT. c_associated(sym)) THEN
       CALL c_message(dlerror_i(), err)
       IF (paral%io_parent) WRITE(6, '(A)') TRIM(err)
       CALL stopgm('ionic_optim', 'missing optimizer symbol', &
            __LINE__, __FILE__)
    END IF
  END FUNCTION load_sym

  SUBROUTINE ionic_rgsaddle()
    CHARACTER(*), PARAMETER :: procedureN = 'ionic_rgsaddle'
    CHARACTER(len=512) :: lib
    CHARACTER(len=256) :: err
    TYPE(c_ptr) :: handle
    TYPE(c_funptr) :: fp
    PROCEDURE(rs_band_create_i), POINTER :: band_create
    PROCEDURE(rs_band_step_i), POINTER :: band_step
    PROCEDURE(rs_band_pos_i), POINTER :: band_pos
    PROCEDURE(rs_band_free_i), POINTER :: band_free
    PROCEDURE(rs_min_create_i), POINTER :: min_create
    PROCEDURE(rs_min_step_i), POINTER :: min_step
    PROCEDURE(rs_min_pos_i), POINTER :: min_pos
    PROCEDURE(rs_min_free_i), POINTER :: min_free
    TYPE(rs_band_config), TARGET :: bcfg
    TYPE(rs_min_config), TARGET :: mcfg
    TYPE(rs_report), TARGET :: report
    TYPE(c_ptr) :: session
    INTEGER :: nat, nimg, ismode
    INTEGER(c_int) :: rc
    LOGICAL :: exists
    CALL lib_name('OPENCPMD_RGSADDLE_LIB', 'librgsaddle.so', lib)
    handle = dlopen_i(TRIM(lib) // c_null_char, rtld_lazy)
    IF (.NOT. c_associated(handle)) THEN
       CALL c_message(dlerror_i(), err)
       IF (paral%io_parent) WRITE(6, '(A)') TRIM(err)
       CALL stopgm(procedureN, 'rgsaddle library did not load', &
            __LINE__, __FILE__)
    END IF
    nat = nat_deck()
    CALL fill_species_and_cell()
    INQUIRE(file='rgsaddle_mode', exist=exists)
    ismode = 0
    IF (exists) ismode = 1
    CALL mp_bcast(ismode, parai%io_source, parai%cp_grp)
    IF (ismode == 1) THEN
       CALL read_mode(nat, mcfg)
       fp = load_as_fun(handle, 'rgsaddle_minmode_create')
       CALL c_f_procpointer(fp, min_create)
       fp = load_as_fun(handle, 'rgsaddle_minmode_step')
       CALL c_f_procpointer(fp, min_step)
       fp = load_as_fun(handle, 'rgsaddle_minmode_position')
       CALL c_f_procpointer(fp, min_pos)
       fp = load_as_fun(handle, 'rgsaddle_minmode_free')
       CALL c_f_procpointer(fp, min_free)
       session = min_create(c_loc(mcfg), INT(nat, c_int64_t), &
            c_loc(rs_pos(1)), c_loc(rs_pos(3 * nat + 1)))
       IF (.NOT. c_associated(session)) &
            CALL stopgm(procedureN, 'rgsaddle rejected the mode', &
            __LINE__, __FILE__)
       DO
          report%version%major = 1
          report%version%minor = 0
          rc = min_step(session, c_funloc(rs_surface), c_null_ptr, c_loc(report))
          IF (rc /= 0) CALL stopgm(procedureN, 'rgsaddle step failed', &
               __LINE__, __FILE__)
          IF (report%status == 1) EXIT
       END DO
       rc = min_pos(session, c_loc(rs_pos(1)))
       CALL min_free(session)
    ELSE
       INQUIRE(file='rgsaddle_task', exist=exists)
       IF (.NOT. exists) CALL stopgm(procedureN, &
            'rgsaddle_task or rgsaddle_mode is required', &
            __LINE__, __FILE__)
       CALL read_band(nat, nimg, bcfg)
       fp = load_as_fun(handle, 'rgsaddle_band_create')
       CALL c_f_procpointer(fp, band_create)
       fp = load_as_fun(handle, 'rgsaddle_band_step')
       CALL c_f_procpointer(fp, band_step)
       fp = load_as_fun(handle, 'rgsaddle_band_positions')
       CALL c_f_procpointer(fp, band_pos)
       fp = load_as_fun(handle, 'rgsaddle_band_free')
       CALL c_f_procpointer(fp, band_free)
       session = band_create(c_loc(bcfg), INT(nimg, c_int64_t), &
            INT(nat, c_int64_t), c_loc(rs_pos(1)))
       IF (.NOT. c_associated(session)) &
            CALL stopgm(procedureN, 'rgsaddle rejected the band', &
            __LINE__, __FILE__)
       DO
          report%version%major = 1
          report%version%minor = 0
          rc = band_step(session, c_funloc(rs_surface), c_null_ptr, c_loc(report))
          IF (rc /= 0) CALL stopgm(procedureN, 'rgsaddle step failed', &
               __LINE__, __FILE__)
          IF (report%status == 1) EXIT
       END DO
       rc = band_pos(session, c_loc(rs_pos(1)))
       CALL band_free(session)
    END IF
    IF (ismode == 1) CALL tau_from_flat(rs_pos(1:3 * nat))
    rc = dlclose_i(handle)
  END SUBROUTINE ionic_rgsaddle

  FUNCTION load_as_fun(handle, name) RESULT(fp)
    TYPE(c_ptr), INTENT(in) :: handle
    CHARACTER(*), INTENT(in) :: name
    TYPE(c_funptr) :: fp
    fp = TRANSFER(load_sym(handle, name), c_null_funptr)
  END FUNCTION load_as_fun

  SUBROUTINE read_band(nat, nimg, cfg)
    INTEGER, INTENT(in) :: nat
    INTEGER, INTENT(out) :: nimg
    TYPE(rs_band_config), INTENT(out) :: cfg
    INTEGER :: ios, i, nread
    REAL(c_double) :: x, y, z
    nimg = 0
    IF (paral%io_parent) THEN
       OPEN(unit=79, file='rgsaddle_task', status='OLD', iostat=ios)
       IF (ios /= 0) nimg = -1
       IF (nimg == 0) THEN
          READ(79, *, iostat=ios) nimg
          IF (ios /= 0 .OR. nimg < 2) nimg = -1
       END IF
       IF (nimg > 0) THEN
          cfg%version%major = 1
          cfg%version%minor = 0
          cfg%flags = 0
          READ(79, *, iostat=ios) cfg%tangent, cfg%spring, cfg%projection, cfg%method
          IF (ios /= 0) nimg = -1
       END IF
       IF (nimg > 0) THEN
          READ(79, *, iostat=ios) cfg%spring_k, cfg%force_tol, cfg%max_move, cfg%memory
          IF (ios /= 0) nimg = -1
       END IF
       IF (nimg > 0) THEN
          IF (ALLOCATED(rs_pos)) DEALLOCATE(rs_pos)
          ALLOCATE(rs_pos(3 * nat * nimg))
          nread = 0
          DO i = 1, nat * nimg
             READ(79, *, iostat=ios) x, y, z
             IF (ios /= 0) THEN
                nimg = -1
                EXIT
             END IF
             rs_pos(nread + 1) = x
             rs_pos(nread + 2) = y
             rs_pos(nread + 3) = z
             nread = nread + 3
          END DO
       END IF
       IF (ios == 0) CLOSE(79)
    END IF
    CALL mp_bcast(nimg, parai%io_source, parai%cp_grp)
    IF (nimg < 2) CALL stopgm('read_band', 'rgsaddle_task is incomplete', &
         __LINE__, __FILE__)
    CALL mp_bcast(cfg%tangent, parai%io_source, parai%cp_grp)
    CALL mp_bcast(cfg%spring, parai%io_source, parai%cp_grp)
    CALL mp_bcast(cfg%projection, parai%io_source, parai%cp_grp)
    CALL mp_bcast(cfg%method, parai%io_source, parai%cp_grp)
    CALL mp_bcast(cfg%spring_k, parai%io_source, parai%cp_grp)
    CALL mp_bcast(cfg%force_tol, parai%io_source, parai%cp_grp)
    CALL mp_bcast(cfg%max_move, parai%io_source, parai%cp_grp)
    CALL mp_bcast(cfg%memory, parai%io_source, parai%cp_grp)
    IF (.NOT. paral%io_parent) THEN
       IF (ALLOCATED(rs_pos)) DEALLOCATE(rs_pos)
       ALLOCATE(rs_pos(3 * nat * nimg))
    END IF
    CALL mp_bcast(rs_pos, 3 * nat * nimg, parai%io_source, parai%cp_grp)
    cfg%version%major = 1
    cfg%version%minor = 0
    cfg%flags = 0
    cfg%spring_ks = c_null_ptr
    cfg%ci_trigger_factor = 0.0_c_double
    cfg%ci_trigger_force = 0.0_c_double
    cfg%cell = c_loc(deck_cell)
  END SUBROUTINE read_band

  SUBROUTINE read_mode(nat, cfg)
    INTEGER, INTENT(in) :: nat
    TYPE(rs_min_config), INTENT(out) :: cfg
    INTEGER :: ios, i, ok
    REAL(c_double) :: x, y, z
    ok = 1
    IF (paral%io_parent) THEN
       OPEN(unit=79, file='rgsaddle_mode', status='OLD', iostat=ios)
       IF (ios /= 0) ok = 0
       IF (ok == 1) THEN
          READ(79, *, iostat=ios) cfg%kind, cfg%method
          IF (ios /= 0) ok = 0
       END IF
       IF (ok == 1) THEN
          READ(79, *, iostat=ios) cfg%dr, cfg%rotation_tol, cfg%max_rotations, &
               cfg%krylov_dim, cfg%force_tol, cfg%max_move
          IF (ios /= 0) ok = 0
       END IF
       IF (ok == 1) THEN
          IF (ALLOCATED(rs_pos)) DEALLOCATE(rs_pos)
          ALLOCATE(rs_pos(6 * nat))
          DO i = 1, 2 * nat
             READ(79, *, iostat=ios) x, y, z
             IF (ios /= 0) THEN
                ok = 0
                EXIT
             END IF
             rs_pos(3 * (i - 1) + 1) = x
             rs_pos(3 * (i - 1) + 2) = y
             rs_pos(3 * (i - 1) + 3) = z
          END DO
       END IF
       IF (ok == 1) CLOSE(79)
    END IF
    CALL mp_bcast(ok, parai%io_source, parai%cp_grp)
    IF (ok /= 1) CALL stopgm('read_mode', 'rgsaddle_mode is incomplete', &
         __LINE__, __FILE__)
    CALL mp_bcast(cfg%kind, parai%io_source, parai%cp_grp)
    CALL mp_bcast(cfg%method, parai%io_source, parai%cp_grp)
    CALL mp_bcast(cfg%dr, parai%io_source, parai%cp_grp)
    CALL mp_bcast(cfg%rotation_tol, parai%io_source, parai%cp_grp)
    CALL mp_bcast(cfg%max_rotations, parai%io_source, parai%cp_grp)
    CALL mp_bcast(cfg%krylov_dim, parai%io_source, parai%cp_grp)
    CALL mp_bcast(cfg%force_tol, parai%io_source, parai%cp_grp)
    CALL mp_bcast(cfg%max_move, parai%io_source, parai%cp_grp)
    IF (.NOT. paral%io_parent) THEN
       IF (ALLOCATED(rs_pos)) DEALLOCATE(rs_pos)
       ALLOCATE(rs_pos(6 * nat))
    END IF
    CALL mp_bcast(rs_pos, 6 * nat, parai%io_source, parai%cp_grp)
    cfg%version%major = 1
    cfg%version%minor = 0
    cfg%flags = 0
  END SUBROUTINE read_mode

  SUBROUTINE ionic_gpr()
    CHARACTER(*), PARAMETER :: procedureN = 'ionic_gpr'
    CHARACTER(len=512) :: lib
    CHARACTER(len=256) :: err
    TYPE(c_ptr) :: handle
    TYPE(c_funptr) :: fp
    PROCEDURE(gpr_run_i), POINTER :: gpr_run
    TYPE(ionic_task), TARGET :: task
    INTEGER :: nat
    INTEGER(c_int) :: rc
    REAL(c_double), ALLOCATABLE, TARGET :: pos(:)
    CALL lib_name('OPENCPMD_GPR_LIB', 'libgpr_optim.so', lib)
    handle = dlopen_i(TRIM(lib) // c_null_char, rtld_lazy)
    IF (.NOT. c_associated(handle)) THEN
       CALL c_message(dlerror_i(), err)
       IF (paral%io_parent) WRITE(6, '(A)') TRIM(err)
       CALL stopgm(procedureN, 'gpr_optim library did not load', &
            __LINE__, __FILE__)
    END IF
    fp = load_as_fun(handle, 'opencpmd_gpr_run')
    CALL c_f_procpointer(fp, gpr_run)
    nat = nat_deck()
    CALL fill_species_and_cell()
    ALLOCATE(pos(3 * nat))
    CALL flat_from_tau(pos)
    task%major = 1
    task%minor = 0
    task%n_atoms = nat
    task%positions = c_loc(pos(1))
    task%atomic_numbers = c_loc(deck_z(1))
    task%cell = c_loc(deck_cell)
    task%eval = c_funloc(opencpmd_ionic_eval)
    rc = gpr_run(c_loc(task))
    IF (dlclose_i(handle) /= 0) CONTINUE
    IF (rc /= 0) CALL stopgm(procedureN, 'gpr_optim returned an error', &
         __LINE__, __FILE__)
  END SUBROUTINE ionic_gpr

END MODULE ionic_optim_utils
