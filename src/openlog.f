      subroutine openlog(lunit,fnm,lfnm, stamp, lstamp)
      integer lunit, lfnm, lstamp, I
      character*(*) fnm, stamp

      open(lunit, file=fnm(1:lfnm))
      call xsetf(1)
      call xsetun(lunit)
      write(lunit, '(1x,80a1)') (stamp(I:I), I = 1, lstamp)
      return
      end
      
      subroutine closelog(lunit)
      integer lunit

      write(lunit,'(1x,8HAll Done)')
      close(lunit)
      call xsetf(0)
      return
      end
