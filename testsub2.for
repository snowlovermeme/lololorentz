      Subroutine dload(f, kstep, kinc, time, noel, npt, layer, kspt, 
     1coords, jltyp, sname)

c    
      include 'ABA_PARAM.INC'
c      
c    形参定义
      DIMENSION TIME(2), COORDS (3)
      Character *80 sname
      Real :: f
      Integer :: noel, kinc, kspt, jltyp, kstep, layer, npt
      
c     内参定义
c    形参说明：
c	X：双精度实型一维数组，存放矩形区域X方向的N个坐标值
c	Y：双精度实型一维数组，存放矩形区域Y方向的M个坐标值
c	Z: 双精度实型二维数组，体积为NXM,输入参数，存放矩形区域上的N*M个节点上的函数值
c	N：整型变量，输入参数，矩形区域X方向的坐标个数
c	M: 整型变量，输入参数，矩形区域X方向的坐标个数
c	U，V：双精度实型变量，输入参数，存放指定插值点的X、Y坐标
c	W：双精度实型变量，输出参数，存放指定（U，Y）处的函数近似值。
      Real :: z(86, 431),x(431),y(76) 
      Real :: uu,fresult
      Integer :: rc,i,j
      Character :: blankLine
      Dimension b(10)
      Real :: u, v, w, b, hh
c     Implicit None
       jltyp=2
c     读取Z文件信息
      Open (101,FILE="D:\abaqusTemp\Fz.txt",STATUS='OLD',ACTION='READ')
 
      Read (101, *) x
C        Write (*, *) x
      Read (101, *) y
C        Write (*, *) y
      Do j = 1, 75
        Read (101, *,IOSTAT=ios) z(:, j)
C          Write (*, *) z(:, j)
        if (ios.NE.0) then
            write(*,*) "******some thing wrong happened****"
            f=10000
            exit
          endif
      EndDo
       close(101)
       write(*,*) "***************read load***********************"
      x = x/1000
      y = y/1000
C       Close (101)
C c     读取R文件信息
C c
C c
C c
C c
C c
C c
C c     开始计算力的大小,插值计算
      u = coords(1)
      v = coords(2)
      


      If (u<=x(1)) Then
        ip = 1
        ipp = 4
      Else If (u>=x(n)) Then
      ip = n - 3
      ipp = n
      Else
        i = 1
        j = n
10      If (iabs(i-j)/=1) Then
          l = (i+j)/2
          If (u<x(l)) Then
            j = l
          Else
            i = l
          End If
          Goto 10
        End If
        ip = i - 3
        ipp = i + 4
      End If
      If (ip<1) ip = 1
      If (ipp>n) ipp = n
      If (v<=y(1)) Then
        iq = 1
        iqq = 4
      Else If (v>=y(m)) Then
        iq = m - 3
        iqq = m
      Else
        i = 1
        j = m

20    If (iabs(j-i)/=1) Then
          l = (i+j)/2
          If (v<y(l)) Then
            j = l
          Else
            i = l
          End If
          Goto 20
        End If
        iq = i - 3
        iqq = i + 4
      End If
      If (iq<1) iq = 1
      If (iqq>m) iqq = m
      Do i = ip, ipp
        b(i-ip+1) = 0.0
        Do j = iq, iqq
          hh = z(i, j)
          Do k = iq, iqq
            If (k/=j) Then
              hh = hh*(v-y(k))/(y(j)-y(k))
            End If
          End Do
          b(i-ip+1) = b(i-ip+1) + hh
        End Do
      End Do
      w = 0.0
      Do i = ip, ipp
        hh = b(i-ip+1)
        Do j = ip, ipp
          If (j/=i) Then
            hh = hh*(u-x(j))/(x(i)-x(j))
          End If
        End Do
        w = w + hh
      End Do
      f = w
       write(*,*) "*************get load ******************"
      
      Return
      End 




