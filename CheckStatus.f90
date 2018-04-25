subroutine check(status)

    use netcdf

    Implicit none

    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
        print *, trim(nf90_strerror(status))
        stop 2
    end if
end subroutine check
