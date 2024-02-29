! file input/out functions and subroutines, for genotype files and allele frequencies.
!
! Jisca Huisman, jisca.huisman@gmail.com
!
! This code is available under GNU General Public License v3
!
!===============================================================================
module sqa_fileIO
  use sqa_general
  implicit none
   
  character(len=3), dimension(4) :: valid_formats = (/'SEQ', 'PED', 'RAW', 'LMT'/)
  character(len=1), allocatable :: alleles(:,:)
  
contains
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! check if file exists
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine checkFile(FileName)
    character(len=*), intent(IN) :: FileName
    logical :: file_exists
    
    inquire(file=trim(FileName), exist = file_exists)
    if (.not. file_exists)  stop 'File '//trim(FileName)//' not found!'
  
  end subroutine checkFile
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! handle IOSTAT errors
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine IOstat_handler(ios, x, FileName)
    integer, intent(IN) :: ios, x
    character(len=*), intent(IN) :: FileName
    
    if (ios > 0) then
      print *, "ERROR: Wrong input in file "//trim(FileName)//" on line ", x
      stop
    else if (ios < 0) then
      print *, "ERROR: file "//trim(FileName)//" has too few rows"
      stop
    endif
  
  end subroutine IOstat_handler


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! determine the number of columns in a file
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  integer function FileNumCol(FileName, sep, has_header)
      character(len=*), intent(IN) :: FileName
      character(len=*), intent(IN), optional :: sep
      logical, intent(IN), optional :: has_header   ! skip header (if any), which has more characters due to long SNP names
      integer :: j, strLen, numcol, ios
      character(len=:), allocatable :: line
      character(len=1) :: s(2)
      logical :: skip_header
      
      if (present(sep)) then   ! separator
        s = (/ sep, '' /)
      else
        s = (/ ' ', achar(9) /)  ! achar(9) = \t (tab)
      endif
      
      if (present(has_header)) then
        skip_header = has_header
      else
        skip_header = .FALSE.
      endif
      
      allocate(character(len=500000) :: line)

      open(unit=102, file=trim(FileName), status="old")
      if (skip_header) then
        read(102,*,IOSTAT=ios)   
        if (ios < 0) then
          FileNumCol = 0
          return
        endif
      endif
      read(102, '(a)',IOSTAT=ios) line
      if (ios < 0) then
        FileNumCol = 0
        return
      endif
      close(102) 

      strLen = len_trim(line)
      if (strLen >= 500000)  print *, 'WARNING: '//trim(FileName)//' EXCEEDING MAXIMUM NUMBER OF COLUMNS!'
      
      if (strLen  == 0) then
        FileNumCol = 0
        return
      endif
      
      if (all(s=='')) then
        FileNumCol = strLen
        return
      endif

      numcol = 0   ! first column (no space 'after')  achar(9) = \t
      do j=1, strLen-1
        if (j==1 .and. .not. any(s == line(j:j))) then
          numcol = numcol +1
        endif
        if (any(s == line(j:j))) then
          if (.not. any(s == line((j+1):(j+1)))) then
            numcol = numcol +1    ! new column starts at j+1
          endif
        endif
      enddo
      FileNumCol = numcol

    end function FileNumCol


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! determine the number of rows in a file
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  integer function FileNumRow(FileName)
    character(len=*), intent(IN) :: FileName
    integer :: nrow, i, maxRow, ios
    character(len=42) :: dumC

    maxRow = 1e8  ! fail safe
    nrow = 0
    open(unit=102, file=trim(FileName), status="old", action='read')
      do i=1, maxRow
        read(102,*,IOSTAT=ios) dumC
        if (ios < 0) then
          exit  ! end of file
        else
          nrow = nrow +1  
        end if
      enddo
    close(102)
    if (nrow >= maxRow) print *, 'WARNING: '//trim(FileName)//' EXCEEDING MAXIMUM NUMBER OF ROWS!'
    FileNumRow = nrow

  end function FileNumRow  
  
 
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! add file extension to genotype file
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  pure function add_extension(FileName, FileFormat)
    character(len=*), intent(IN) :: FileName
    character(len=3), intent(IN) :: FileFormat
    character(len=len(FileName)+5) :: add_extension
    
    select case (FileFormat)
      case ('SEQ')
        add_extension = trim(FileName)//'.txt'
      case('PED')
        add_extension = trim(FileName)//'.ped'
      case('RAW')
        add_extension = trim(FileName)//'.raw'
      case('LMT')
        add_extension = trim(FileName)//'.geno'
    end select   
  
  end function add_extension
  
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! read genotype file
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine read_geno(Geno,ID, SNP_names, FileName, FileFormat, transp)  
    integer(kind=ishort), allocatable, intent(OUT) :: Geno(:,:)
    character(len=nchar_ID), allocatable, intent(OUT) :: ID(:), SNP_names(:)
    character(len=*), intent(IN) :: FileName
    character(len=3), intent(IN) :: FileFormat
    logical, intent(IN), optional :: transp   ! transpose matrix (default=yes)
    integer :: nSnp, nInd, i, l, ios,m
    integer(kind=ishort), allocatable :: G_int(:), Geno_tmp(:,:)
    character(len=LEN(FileName)+5) :: GenoFile
    character(len=2), allocatable :: G_char(:)
    character(len=1), allocatable :: G_duos(:,:)
    
    character(len=nchar_ID) :: dumC, dumV(4)
    logical :: do_transpose, mapfile_found
    
    ! check input
    if (.not. any(valid_formats == FileFormat)) then
      print *, 'ERROR: Genotype FileFormat must be one of: ', valid_formats
      stop
    endif    
    GenoFile = add_extension(FileName, FileFormat)
    call CheckFile(GenoFile)
          
    ! determine size of genotype matrix
    nInd = FileNumRow(GenoFile) 
    if (FileFormat == 'RAW') then
      nInd = nInd -1   ! header row
    endif
    allocate(ID(1:nInd))
    ID = "NA"
    
    select case (FileFormat)
      case ('SEQ')
        nSnp = FileNumCol(trim(GenoFile)) -1  ! column 1 = IDs
      case('PED')
        nSnp = (FileNumCol(trim(GenoFile)) -6)/2  ! columns FID-IID-sire-dam-sex-pheno
        allocate(G_duos(1:nInd, nSnp*2))   ! 2 columns per SNP
      case('RAW')
        nSnp = FileNumCol(trim(GenoFile), has_header=.TRUE.) -6
        allocate(G_char(nSnp))   ! uses NA for missing values
      case('LMT')
        nSnp = FileNumCol(trim(GenoFile), sep='')
    end select  
  
    allocate(G_int(nSnp))
    allocate(Geno(0:nInd, nSnp))   
    Geno = -1
    
    allocate(SNP_names(nSnp))
    SNP_names = ''
    
    ! read in genotype data
    open (unit=101,file=trim(GenoFile),status="old", action='read')
      if (FileFormat=='RAW') then
        read (101,*) dumC, dumC, dumV, SNP_names   ! header
      endif
      ! TODO: select once, not every iteration (speed effect?)
      do i=1,nInd
        select case (FileFormat)
          case ('SEQ')
            read (101,*,IOSTAT=ios)  Id(i), G_int
          case('PED')
            read (101,*,IOSTAT=ios) dumC, Id(i), dumV, G_duos(i,:)
          case('RAW')
            read (101,*,IOSTAT=ios) dumC, Id(i), dumV, G_char
            WHERE (G_char == 'NA')  G_char = '-9'
            do l=1,nSnp
              read(G_char(l), '(i2)',IOSTAT=ios) G_int(l)
            enddo
          case('LMT')
            read (101,'(500000i1)',IOSTAT=ios)  G_int
        end select
        if (ios /= 0)  call IOstat_handler(ios, i, trim(GenoFile))

        WHERE (G_int >= 0 .and. G_int <= 2)  Geno(i,:) = G_int     
      enddo
    close (101)
    
    if (FileFormat=='PED') then
      call Two2One()
    endif
    
    deallocate(G_int)
    if(allocated(G_duos))  deallocate(G_duos)
    if(allocated(G_char))  deallocate(G_char)
       
    ! read in auxiliary files with SNP names (PED) or IDs (LMT)
    mapfile_found = .FALSE.
    if (FileFormat=='PED')  inquire(file=trim(FileName)//'.map', exist = mapfile_found)
    if (FileFormat=='PED' .and. mapfile_found) then 
      open (unit=303, file=trim(FileName)//'.map', status='old', action='read')
        do l=1,nSnp
          read (303,*,IOSTAT=ios)  dumC, SNP_names(l) ! SNP names in 2nd column
          if (ios /= 0)  call IOstat_handler(ios, l, trim(FileName)//'.map')
        enddo     
      close(303)
      
    else if (FileFormat /= 'RAW') then  ! RAW: SNP names are header in genotype file
      ! create fake SNP names
      do l=1,nSnp
        write(SNP_names(l), '("SNP",i6.6)') l 
      enddo   
    endif
    
    if (FileFormat == 'LMT') then
      call CheckFile(trim(FileName)//'.id')
      open (unit=303, file=trim(FileName)//'.id', status='old', action='read')
        do i=1,nInd
          read (303,*,IOSTAT=ios)  ID(i) 
          if (ios /= 0)  call IOstat_handler(ios, i, trim(FileName)//'.id')
        enddo
      close(303) 
    endif
   
    ! transpose
    if (present(transp)) then
      do_transpose = transp
    else
      do_transpose = .TRUE.
    endif
    if (do_transpose) then
      Geno_tmp = transpose(Geno)
      call move_alloc(Geno_tmp, Geno)   ! from, to
    endif
    
    
  contains
    !~~~~~~~~~~~~~~
    ! From 2-columns-per-SNP to 1-column-per-SNP
    subroutine Two2One()
      integer(kind=ishort), allocatable :: Gint(:)
      character(len=1) :: valid_alleles(6), Gi(2)
      character(len=1), allocatable :: Gl(:,:)
      integer :: x,a,i
      
      valid_alleles = (/'A','C','T','G','1','2'/)
      allocate(alleles(2,nSnp))
      alleles = '0'
      allocate(Gl(2,1:nInd))
      allocate(Gint(1:nInd))
      
      do l=1, nSnp
        m = 2*(l-1)+1
        Gl = transpose(G_duos(:,m:(m+1)))
             
        a = 0
        do x=1,6
          if (any(Gl == valid_alleles(x))) then
            a = a+1
            alleles(a,l) = valid_alleles(x)
            if (a==2)  exit
          endif
        enddo
      
        ! TODO: try various alternatives to find fastest way
        do i=1,nInd
          Gi = Gl(:,i)
          if (all(Gi == '0')) then
            Gint(i) = -1
          else if (all(Gi==alleles(1,l))) then
            Gint(i) = 0
          else if (all(Gi==alleles(2,l))) then
            Gint(i) = 2
          else if (any(Gi==alleles(1,l)) .and. any(Gi==alleles(2,l))) then
            Gint(i) = 1
          else  ! shouldn't happen
            Gint(i) = -1
          endif
        enddo 
        Geno(1:nInd,l) = Gint
      enddo      
 !     print *, 'genotype counts:', COUNT(Gint==0), COUNT(Gint==1), COUNT(Gint==2), COUNT(Gint==-1)
      
    end subroutine Two2One

  end subroutine read_geno
  
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! write genotype file
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine write_geno(Geno, nInd, nSnp, ID, SNP_names, FileName, FileFormat, make_map)
    integer, intent(IN) :: nInd, nSnp
    integer(kind=ishort), intent(IN) :: Geno(1:nSnp, 1:nInd)
    character(len=nchar_ID), intent(IN) :: ID(nInd)
    character(len=nchar_ID), intent(IN), optional :: SNP_names(nSnp)
    logical, intent(IN), optional :: make_map
    character(len=*), intent(IN) :: FileName
    character(len=3), intent(IN) :: FileFormat
    character(len=LEN(FileName)+5) :: GenoFile
    character(len=nchar_ID) :: SNP_names_out(nSnp)
    character(len=2), allocatable :: G_char(:)
    character(len=1), allocatable :: G_duos(:,:)
    integer :: i, l
    integer(kind=ishort), allocatable :: G_int(:)
    logical :: did_warn, do_make_map
  
    if (.not. any(valid_formats == FileFormat)) then
      print *, 'ERROR: Genotype FileFormat must be one of: ', valid_formats
      stop
    endif    
    GenoFile = add_extension(FileName, FileFormat)
    
    if (FileFormat == 'RAW')  allocate(G_char(nSnp))   ! uses NA for missing values
    
    if (present(make_map)) then
      do_make_map = make_map
    else if (FileFormat=='PED') then
      do_make_map = .TRUE.
    else
      do_make_map = .FALSE.
    endif
    
    if (present(SNP_names)) then
      SNP_names_out = SNP_names
    else if (FileFormat == 'RAW' .or. do_make_map) then
      ! make fake SNP names
      do l=1,nSnp
        write(SNP_names_out(l), '("SNP",i6.6)') l 
      enddo    
    endif
    
    if (FileFormat=='PED') then
      allocate(G_duos(1:(nSnp*2), 1:nInd))     
      call One2Two()
    endif
    
    did_warn = .FALSE.

    open (unit=202, file=trim(GenoFile), status='unknown', action='write')
      if (FileFormat=='RAW') then   ! header
        write (202,'(a4,2x,a3,37x,3a5, a6, 2x, 200000a25)') 'FID', 'IID', 'PAT', &
         'MAT', 'SEX', 'PHENO', SNP_names_out
      endif
      do i=1,nInd
        G_int = Geno(:,i)
        select case (FileFormat)
          case ('SEQ')
            write(202, '(a40, 100000i3)') Id(i), G_int
          case('PED')
            write(202, '(i3,2x,a40,4i3,2x, 200000a2)') 0, Id(i), 0,0,0,-9, G_duos(:,i) 
          case('RAW')
            do l=1,nSnp
              write(G_char(l), '(i2)')  G_int(l)
            enddo
            WHERE (G_char == '-9')  G_char = 'NA'
            write(202, '(i4,2x,a40,4i5,2x, 200000a3)') 0, Id(i), 0,0,0,0, G_char
          case('LMT')
            if (any(G_int < 0) .and. .not. did_warn) then  ! warning once is enough
              print *, 'WARNING: LMT does not support missing values! Coded as 9 in output'
              did_warn = .TRUE.         
            endif
            WHERE (G_int < 0)  G_int = 9
            write(202,'(500000i1)')  G_int
        end select 
      enddo
    close(202)
    
    deallocate(G_int)
    if(allocated(G_char))  deallocate(G_char)
        
    ! write auxiliary files
    if (do_make_map) then  ! typically in combination with FileFormat=='PED'
      open(505, file=trim(FileName)//'.map', status='unknown', action='write')
        ! no header
        do l=1,nSnp
          write(505,'(i2,1x, a25, 2i2)') 0, SNP_names_out(l), 0, 0
        enddo
      close(505)
    endif
    
    if (FileFormat=='LMT') then
      open(505, file=trim(FileName)//'.id', status='unknown', action='write')
        do i=1,nInd
          write(505,'(a40)') ID(i)
        enddo
      close(505)
    endif
    
    
  contains
    !~~~~~~~~~~~~~~
    ! From 1-column-per-SNP to 2-columns-per-SNP
    subroutine One2Two() 
      character(len=1), allocatable :: Gl(:,:)
      integer :: l,m,i
      
      allocate(Gl(2,1:nInd))
      
      do l=1, nSnp
        m = 2*(l-1)+1
        
        if (allocated(alleles)) then
          do i=1,nInd
            select case (Geno(l,i))
              case (0)
                Gl(:,i) = (/alleles(1,l), alleles(1,l)/)
              case (1)
                Gl(:,i) = (/alleles(1,l), alleles(2,l)/) 
              case (2)
                Gl(:,i) = (/alleles(2,l), alleles(2,l)/)
              case default 
                Gl(:,i) = (/'0','0'/)     ! missing
            end select
          enddo
        else
          do i=1,nInd
            select case (Geno(l,i))
              case (0)
                Gl(:,i) = (/'1','1'/)
              case (1)
                Gl(:,i) = (/'1','2'/) 
              case (2)
                Gl(:,i) = (/'2','2'/)
              case default 
                Gl(:,i) = (/'0','0'/)     ! missing
            end select   
          enddo
        endif
        G_duos(m:(m+1),:) = Gl
      enddo      
    end subroutine One2Two
  
  end subroutine write_geno
  

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! read allele frequency file
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function readAF(FileName)
    character(len=*), intent(IN), optional :: FileName
    double precision, allocatable :: readAF(:)
    integer :: l, nCol, nRow, AFcol, k, ios
    character(len=50), allocatable :: header(:), tmpC(:)
    
    call CheckFile(FileName)

    nCol = FileNumCol(FileName)
    nRow = FileNumRow(FileName)
    if (nCol==1) then
      allocate(readAF(nRow))  ! no header
    else
      allocate(readAF(nRow -1))  ! with header
    endif
    allocate(header(nCol))
    header = 'NA'
    AFcol = 0
    
    readAF = 0D0
    open(unit=3, file=trim(FileName), status="old", action='read')     
      if (nCol == 1) then
        AFcol = 1
      else
        read(103,*)  header
        do k=1, nCol
          if (header(k) == 'MAF' .or. header(k)=='AF' .or. header(k)=='Frequency') then
            AFcol = k
            exit
          endif
        enddo
      endif
      if (AFcol > 1)  allocate(tmpC(AFcol -1))
      
      do l=1, nRow
        if (AFcol == 1) then
          read(3, *,IOSTAT=ios)  readAF(l)
        else
          read(3, *,IOSTAT=ios)  tmpC, readAF(l)
        endif
        if (ios > 0) then
          print *, "Wrong input in file "//trim(FileName)//" on line ", l
          stop
        end if
      enddo
    close(3)
    
    deallocate(header)
    if (allocated(tmpC))  deallocate(tmpC)

  end function readAF
  
  
end module sqa_fileIO

!===============================================================================
! end.