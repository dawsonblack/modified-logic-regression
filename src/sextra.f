      ! *****************************************************************
      ! utility for debugging - printing model
      SUBROUTINE dprintmodel(nkn,ntr,nsep,conc,negs,pick,term,
     #    betas,score,k2,m2,kk)
      IMPLICIT NONE
      INTEGER nkn,ntr,k2,i,j,k,m2,nsep,j2,kk
      INTEGER conc(nkn,ntr,3)
      INTEGER negs(nkn,ntr,3)
      INTEGER pick(nkn,ntr,3)
      INTEGER term(nkn,ntr,3)
      REAL score,betas(0:(nsep+ntr))
      CHARACTER (LEN=125) aa
      aa(1:11)="###     ###"
      CALL makeistring(5,7,aa,kk,3)
      CALL stringprint(aa,11)
      DO j2=1,ntr
         IF(pick(1,j2,k2).gt.0)THEN
      aa(1:11)="print tree "
      CALL makeistring(12,13,aa,j2,2)
      CALL makeistring(14,15,aa,k2,2)
      aa(16:22)=" out of "
      CALL makeistring(23,24,aa,ntr,2)
      IF(m2.gt.0)THEN
         aa(25:27)=" b "
         CALL makerstring(28,41,aa,betas(nsep+j2),9,4)
         aa(42:44)=" s "
         CALL makerstring(45,58,aa,score,9,4)
         CALL stringprint(aa,58)
      ELSE
         CALL stringprint(aa,24)
      END IF
      END IF
      DO i=1,nkn
         IF(pick(i,j2,k2).gt.0)THEN
         CALL makeistring(1,2,aa,i,2)
         CALL makeistring(3,8,aa,conc(i,j2,k2),6)
         CALL makeistring(9,13,aa,negs(i,j2,k2),5)
         CALL makeistring(14,18,aa,pick(i,j2,k2),5)
         CALL makeistring(19,23,aa,term(i,j2,k2),5)
         CALL stringprint(aa,23)
         END IF
      END DO
      END DO
      END

      SUBROUTINE iprintmodel(i,aa,l)
      IMPLICIT NONE
      INTEGER i,l
      CHARACTER (LEN=125) aa
      CALL makeistring(l+1,l+7,aa,i,7)
      CALL stringprint(aa,l+7)
      END

      SUBROUTINE rprintmodel(r,aa,l)
      IMPLICIT NONE
      INTEGER l
      REAL r
      CHARACTER (LEN=125) aa
      CALL makerstring(l+1,l+14,aa,r,9,4)
      CALL stringprint(aa,l+14)
      END
       
      SUBROUTINE betaprint(sc,b,aax,l)
      IMPLICIT NONE
      REAL sc,b(0:4)
      INTEGER l
      CHARACTER (LEN=l) aax
      CHARACTER (LEN=125) aa
      aa(1:l)=aax
      CALL makerstring(l+1,l+14,aa,sc,9,4)
      CALL makerstring(l+15,l+28,aa,b(0),9,4)
      CALL makerstring(l+29,l+42,aa,b(1),9,4)
      CALL makerstring(l+43,l+56,aa,b(2),9,4)
      CALL makerstring(l+57,l+70,aa,b(3),9,4)
      CALL makerstring(l+71,l+84,aa,b(4),9,4)
      CALL stringprint(aa,l+84)
      END
       
      SUBROUTINE prtrprint(prtr,i,aax,l,n1,ntr)
      IMPLICIT NONE
      INTEGER i,l,n1,ntr,j1,j2
      INTEGER prtr(n1,ntr)
      CHARACTER (LEN=l) aax
      CHARACTER (LEN=125) aa
      aa(1:l)=aax
      DO j1=1,40
         CALL makeistring(j1*2+l-1,j1*2+l,aa,prtr(j1,i),2)
      END DO
      CALL stringprint(aa,l+80)
      END

      SUBROUTINE storprint(stor,aax,l,n1)
      IMPLICIT NONE
      INTEGER l,n1,ntr,j1,j2
      INTEGER stor(n1)
      CHARACTER (LEN=l) aax
      CHARACTER (LEN=125) aa
      aa(1:l)=aax
      DO j1=1,40
         CALL makeistring(j1*2+l-1,j1*2+l,aa,stor(j1),2)
      END DO
      CALL stringprint(aa,l+80)
      END

      SUBROUTINE emprint(aax,l,n1,n2,n3,n4,n5,n6,n7,n8)
      IMPLICIT NONE
      INTEGER l,n1,n2,n3,n4,n5,n6,n7,n8,i
      CHARACTER (LEN=l) aax
      CHARACTER (LEN=125) aa
      aa(1:l)=aax
      CALL stringprint(aa,l)
      IF(n1.ne.999)THEN
          CALL makeistring(l+1,l+3,aa,n1,3)
          l=l+3
      CALL stringprint(aa,l)
      END IF
      IF(n2.ne.999)THEN
          CALL makeistring(l+1,l+3,aa,n2,3)
          l=l+3
      CALL stringprint(aa,l)
      END IF
      IF(n3.ne.999)THEN
          CALL makeistring(l+1,l+3,aa,n3,3)
          l=l+3
      CALL stringprint(aa,l)
      END IF
      IF(n4.ne.999)THEN
          CALL makeistring(l+1,l+3,aa,n4,3)
          l=l+3
      CALL stringprint(aa,l)
      END IF
      IF(n5.ne.999)THEN
          CALL makeistring(l+1,l+3,aa,n5,3)
          l=l+3
      CALL stringprint(aa,l)
      END IF
      IF(n6.ne.999)THEN
          CALL makeistring(l+1,l+3,aa,n6,3)
          l=l+3
      CALL stringprint(aa,l)
      END IF
      IF(n7.ne.999)THEN
          CALL makeistring(l+1,l+3,aa,n7,3)
          l=l+3
      CALL stringprint(aa,l)
      END IF
      IF(n8.ne.999)THEN
          CALL makeistring(l+1,l+3,aa,n8,3)
          l=l+3
      CALL stringprint(aa,l)
      END IF
      CALL stringprint(aa,l)
      END
