module soc_mo
  
  use prec_mod
  use data_input
  use init_prep

  implicit none
    
  real(dpr), allocatable    :: soc_moint(:,:,:)

  
  contains

  subroutine cal_soc_moint(aoint_temp1d, mo_coeff_temp0)
      integer                :: i, k, j
      real(dpr)              :: diag_ele
      real(dpr), intent(in)  :: aoint_temp1d(input%num_bov(1)*input%num_bov(1)*3)
      real(dpr), intent(in)  :: mo_coeff_temp0(input%num_bov(1), ndim1)
      real(dpr)              :: aoint_temp(input%num_bov(1), input%num_bov(1), 3)
      real(dpr)              :: mo_coeff_temp(input%num_bov(1), ndim1)
      real(dpr), allocatable :: tp_mo_coeff_temp(:,:)
     
      print *, "input%num_bov(1), ndim1" , input%num_bov(1), ndim1
      print *, "mo_coeff1", mo_coeff_temp0(input%num_bov(1),ndim1)
      print *, "mo_coeff2", mo_coeff_temp0(input%num_bov(1),1)
      print *, "mo_coeff3", mo_coeff_temp0(1,1)
      print *, "mo_coeff4", mo_coeff_temp0(2,1)
  
      allocate(soc_moint(ndim1, ndim1, 3))
      allocate(tp_mo_coeff_temp(ndim1, input%num_bov(1)))
      
      aoint_temp = reshape(aoint_temp1d, [input%num_bov(1), input%num_bov(1), 3])
      tp_mo_coeff_temp = transpose(mo_coeff_temp0)
      
      k = input%num_bov(1)
      do i=1, 3
        mo_coeff_temp = matmul(aoint_temp(:,:,i), mo_coeff_temp0)
        soc_moint(:,:,i) = matmul(tp_mo_coeff_temp, mo_coeff_temp)
        print *, "aoint_temp_diag", aoint_temp(k,k,i)
        print *, "aoint_temp_symm", aoint_temp(k,k-1,i), aoint_temp(k-1,k,i)
        print *, "aoint_temp_size", size(aoint_temp(:,1,i)), size(aoint_temp(:,:,i))
        print *, "soc_moint_diag",  soc_moint(ndim1,ndim1,i), soc_moint(ndim1/2,ndim1/2,i)
        print *, "soc_moint_size",  size(soc_moint(:,1,i)), size(soc_moint(:,:,i))
        print *, "soc_moint_symm",  soc_moint(ndim1,ndim1-1,i), soc_moint(ndim1-1,ndim1,i)
      enddo
      soc_moint = 0.5*fine_stru**2*soc_moint
      do k=1, 3
        do j=1, ndim1
          do i=j, ndim1
            if(abs(soc_moint(i,j,k)) > 1.0e-4) then
              ! print *, ...
            endif
          enddo
          diag_ele = soc_moint(j,j,k)
          if( abs(diag_ele) > 1.0e-5) then
            print *, "soc mo integral diagional: ", diag_ele
            stop
          endif
        enddo
      enddo
      deallocate(tp_mo_coeff_temp)
  end subroutine cal_soc_moint

end module soc_mo
