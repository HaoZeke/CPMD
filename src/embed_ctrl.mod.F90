MODULE embed_ctrl

  IMPLICIT NONE

  PRIVATE

  ! cpmd.x leaves this true. rwfopt skips zhwwf, and rwfopt and initrun
  ! skip geofile, when a host sets it false. wrgeof still prints
  ! coordinates on the output.
  LOGICAL, PUBLIC, SAVE :: embed_write_files = .TRUE.

END MODULE embed_ctrl
