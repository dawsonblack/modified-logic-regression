! This file contains Fortran code for fitting a trio logic regression model
! written by Ingo Ruczinski and Qing Li (and slightly modified by
! Holger Schwender).

! =====================================================================
! Copyright (C) 2009-2011  Ingo Ruczinski and Qing Li

! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! The text of the GNU General Public License, version 2, is available
! as http://www.gnu.org/copyleft or by writing to the Free Software
! Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

! The main reference for Trio Logic Regression is
! Li Q, Fallin MD, Louis TA, Lasseter VK, McGrath JA, Avramopoulos D,
! Wolyniec PS, Valle D, Liang KY, Pulver AE, Ruczinski I (2010).
! Detection of SNP-SNP Interactions in Trios of Parents with
! Schizophrenic Children. Genetic Epdemiology, 34, 396-406.

! Contact: ingo@jhu.edu
! =======================================================================

      SUBROUTINE myphxxz(delta,idx,covs,np,n1,nsep,ntr,logl,beta,strata,
     #                   reject,wud1,wud3,bmax)
      IMPLICIT none
        ! i/o
          INTEGER n1,nsep,ntr,zolala,np,reject,bmax
          INTEGER delta(n1),idx(n1),strata(n1)
          DOUBLE PRECISION beta(nsep+ntr),logl,covs(n1*(nsep+ntr))
        ! workspaces
          DOUBLE PRECISION wud1(n1*5),wud3(n1*(bmax+bmax*bmax))
        ! local
          DOUBLE PRECISION ologl,nlogl,alpha,pp,alphap1,alphap2,prec
          DOUBLE PRECISION nbeta(bmax),grad(bmax)
          DOUBLE PRECISION hess(bmax,bmax)
          INTEGER i,iter

          DO i=1,np
             beta(i)=0.
          END DO
          alphap1=0.001
          alphap2=0.00001
          alpha=1
          prec=0.00001
          pp=prec+10
          iter=0
          DO WHILE(iter.LT.10 .AND. pp.GT.prec .AND. alpha.GT.alphap2)
            iter=iter+1
            CALL mygradphz(grad,hess,beta,delta,idx,covs,np,n1,ologl,
     #           strata,bmax,
     #                 wud1(1:n1),
     #                 wud1((n1+1):(2*n1)),
     #                 wud1((2*n1+1):(3*n1)),
     #                 wud1((3*n1+1):(4*n1)),
     #                 wud1((4*n1+1):(5*n1)),
     #                 wud3(1:(n1*bmax)),
     #                 wud3((n1*bmax+1):(n1*(bmax+bmax*bmax))))
            DO i=1,np
              IF(hess(i,i).LT.1.0e-10 .AND. hess(i,i).GT.-1.0e-10) THEN
                CALL mypllxxz(logl,beta,delta,idx,covs,np,n1,strata,
     #                 wud1(1:n1),
     #                 wud1((n1+1):(2*n1)),
     #                 wud1((2*n1+1):(3*n1)),
     #                 wud1((3*n1+1):(4*n1)))
                GOTO 1234
              END IF
            END DO
            CALL lusolveph(hess,grad,np,reject,bmax)
            IF(reject.eq.1)RETURN
            alpha=1
            zolala = 0
            DO WHILE ((alpha.GT.alphap2 .AND. nlogl.LT.ologl)
     #                                                 .OR.zolala.EQ.0)
              zolala = 1
              DO i=1,np
                nbeta(i)=beta(i)+alpha*grad(i)
              END DO
              CALL mypllxxz(nlogl,nbeta,delta,idx,covs,np,n1,strata,
     #                 wud1(1:n1),
     #                 wud1((n1+1):(2*n1)),
     #                 wud1((2*n1+1):(3*n1)),
     #                 wud1((3*n1+1):(4*n1)))
              IF(nlogl.LT.ologl) alpha=alpha/2.
            END DO
            IF(alpha.GT.alphap1) THEN
              pp=0.
              DO i=1,np
                pp=pp+((nbeta(i)-beta(i))*(nbeta(i)-beta(i)))
                beta(i)=nbeta(i)
              END DO
              pp=SQRT(pp)
              IF(iter.LT.3) pp=prec+10.
            END IF
          END DO
          CALL mygradphz(grad,hess,beta,delta,idx,covs,np,n1,logl,
     #           strata,bmax,
     #                 wud1(1:n1),
     #                 wud1((n1+1):(2*n1)),
     #                 wud1((2*n1+1):(3*n1)),
     #                 wud1((3*n1+1):(4*n1)),
     #                 wud1((4*n1+1):(5*n1)),
     #                 wud3(1:(n1*bmax)),
     #                 wud3((n1*bmax+1):(n1*(bmax+bmax*bmax))))
1234      CONTINUE
      END
      ! **********************************************************
      ! **********************************************************
      SUBROUTINE mygradphz(grad,hess,beta,delta,idx,covs,np,n1,logl,
     #           strata,bmax,ff,s1s,gg,ff2,s0,s1,s2)
      IMPLICIT none
        ! i/o
          INTEGER n1,np,bmax
          INTEGER delta(n1),idx(n1),strata(n1)
          DOUBLE PRECISION beta(np),grad(np),covs(n1*np),hess(bmax,np)
          DOUBLE PRECISION logl
        ! local
          DOUBLE PRECISION ff(n1),s1s(n1),gg(n1)
          DOUBLE PRECISION ff2(n1),s0(n1),myexp,mylog
          DOUBLE PRECISION s1(bmax*n1),s1r,u,z
          DOUBLE PRECISION s2(bmax*bmax*n1)
          INTEGER i,i2,j,k,r,s,it,sx
          u=0
          DO i=1,n1
            s0(i)=0
            ff(i)=0.
            DO k=1,np
              s1((k-1)*n1+i)=0
              DO j=1,np
                s2(((k-1)*np+j-1)*n1+i)=0
              END DO
              ff(i)=ff(i)+beta(k)*covs(i+n1*(k-1))
            END DO
          END DO
          DO i=1,np
            grad(i)=0.
            DO k=1,np
              hess(i,k)=0.
            END DO
          END DO
          DO i=1,n1
             gg(i)=ff(idx(i))
             ff2(i)=myexp(gg(i))
          END DO
          DO i2=1,n1
            i=n1+1-i2
            j=idx(i)
            sx=strata(j)
            IF(sx.GT.0)THEN
            s0(sx)=s0(sx)+ff2(i)
            DO r=1,np
              s1r=ff2(i)*covs(j+n1*(r-1))
              s1((r-1)*n1+sx)=s1((r-1)*n1+sx)+s1r
              it=(r-1)*np
              DO s=r,np
                s2((it+s-1)*n1+sx)=s2((it+s-1)*n1+sx)
     #                               +s1r*covs(j+n1*(s-1))
              END DO
            END DO
            IF(delta(idx(i)).EQ.1) THEN
              DO r=1,np
                s1s(r)=s1((r-1)*n1+sx)/s0(sx)
              END DO
              DO r=1,np
                it=(r-1)*np
                grad(r)=grad(r)+covs(idx(i)+n1*(r-1))-s1s(r)
                DO s=r,np
                   hess(r,s)=hess(r,s)-s1s(r)*s1s(s)
     #                       +s2((it+s-1)*n1+sx)/s0(sx)
                END DO
              END DO
              z=ff2(i)/s0(sx)
              z=mylog(z)
              u=u+z
            END IF
            END IF
          END DO
          DO r=1,np
            DO s=1,r
              hess(r,s)=hess(s,r)
            END DO
          END DO
          logl=u
      END
      ! ***********************************************************
      ! ***********************************************************
      SUBROUTINE mypllxxz(logl,beta,delta,idx,covs,np,n1,strata,ff,
     #          ff2,gg,s0)
      IMPLICIT none
        ! i/o
          INTEGER n1,np
          INTEGER delta(n1),idx(n1),strata(n1)
          DOUBLE PRECISION beta(np),covs(n1*np),logl
        ! local
          INTEGER i,k,sx
          DOUBLE PRECISION z,ff(n1),ff2(n1),gg(n1)
          DOUBLE PRECISION s0(n1),myexp,mylog
          logl=0.
          DO i=1,n1
            ff(i)=0.
            DO k=1,np
              ff(i)=ff(i)+beta(k)*covs(i+n1*(k-1))
            END DO
          END DO
          DO i=1,n1
             s0(i)=0.
             gg(i)=ff(idx(i))
             ff2(i)=myexp(gg(i))
          END DO
          DO k=1,n1
            i=n1+1-k
            sx=strata(idx(i))
            IF(sx.GT.0)THEN
            s0(sx)=s0(sx)+ff2(i)
            IF(delta(idx(i)).EQ.1) THEN
              z=ff2(i)/s0(sx)
              z=mylog(z)
              logl=logl+z
            END IF
            END IF
          END DO
      END
      ! ********************************************************
      ! ********************************************************

      SUBROUTINE triofitting(prtr,rsp,dchp,ordrs,weight,n1,ntr,
     #               nop,wh,nsep,seps,score,betas,reject,
     #               wud1,covsf,wud3,strata,idx,delta,bmax)
      IMPLICIT NONE

        ! arguments in
          INTEGER n1,nop,nsep,ntr,wh,reject,bmax
          INTEGER dchp(n1),ordrs(n1)
          REAL rsp(n1),weight(n1),seps(nsep,n1)
          INTEGER prtr(n1,ntr)
        ! workspaces
          DOUBLE PRECISION wud1(5*n1),wud3(n1*(bmax+bmax*bmax))
        ! local
          INTEGER i,j,k,l,m,nnf(2)
          INTEGER myausfahrt,myausfahrt0,myausfahrt1,myicheck
          INTEGER strata(n1),idx(n1),delta(n1)
          DOUBLE PRECISION loglf,betaf(bmax)
          DOUBLE PRECISION covsf(n1*bmax)
        ! arguments out
          REAL score(3),betas(0:(nsep+ntr)),r
        DO i=1,n1
           IF((dchp(i).NE.0).and.(dchp(i).NE.1)) THEN
              CALL REXIT("Response not correctly specified.")
           END IF
        END DO
        ordrs(1)=ordrs(1)
        weight(1)=weight(1)
        wh=wh

       ! reorganize the data
         j=0
         k=0
         DO i=1,n1
            idx(i)=i
            IF(j.EQ.0)THEN
               l=INT(rsp(i))
               IF(l.GT.0)THEN
                  k=k+1
                  j=l
                  delta(i)=1
                  strata(i)=k
               ELSE
                  delta(i)=0
                  strata(i)=-1
               END IF
            ELSE
               j=j-1
               delta(i)=0
               strata(i)=k
            END IF
         END DO

       ! convergence check

         myicheck=0
         myausfahrt=0
         myausfahrt0=0
         myausfahrt1=0

         DO i=1,n1
            myicheck=myicheck+prtr(i,1)
         END DO

         IF (myicheck.GT.0) THEN
            DO i=1,n1
               IF (rsp(i).GT.0) THEN
                  m=INT(rsp(i))
                  DO j=1,m
                     IF (prtr(i,1).NE.prtr(i+j,1)) THEN
                        IF (prtr(i,1).EQ.0) THEN
                           myausfahrt0=1
                        ELSE
                           myausfahrt1=1
                        END IF
                     END IF
                  END DO
               END IF
            END DO
            IF ((myausfahrt0.EQ.0).OR.(myausfahrt1.EQ.0)) THEN
               myausfahrt=1
            END IF
         END IF

         nnf(1)=nop+nsep
         nnf(2)=n1
         DO i=1,(n1*(nsep+ntr))
            covsf(i)=0.D0
         END DO
         IF (nnf(1).GT.0) THEN
            IF (nsep.GT.0) THEN
               DO k=1,nsep
                  DO j=1,n1
                     covsf((k-1)*n1+j)=seps(k,j)
                  END DO
               END DO
            END IF
            DO k=(nsep+1),nnf(1)
               DO j=1,n1
                  covsf((k-1)*n1+j)=REAL(prtr(j,k-nsep))
               END DO
            END DO
         END IF
         reject=0


      ! calculate partial likelihood
         CALL myphxxz(delta,idx,covsf,nnf(1),n1,nsep,ntr,
     #           loglf,betaf,strata,reject,wud1,wud3,bmax)

         r=n1
         score(1)=REAL(-loglf/r)
         DO i=1,(nsep+nop)
            betas(i)=REAL(betaf(i))
         END DO

         IF (myausfahrt.EQ.1) THEN
            reject=1
         END IF

      END