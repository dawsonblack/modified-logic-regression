! =====================================================================
! Copyright (C) 1999-2005  Ingo Ruczinski and Charles Kooperberg

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

! The main reference for Logic Regression is
! Ruczinski I, Kooperberg C, LeBlanc M (2003). Logic Regression,
! Journal of Computational and Graphical Statistics, 12, 475-511
! Other references can be found on our homepages
! http://www.biostat.jhsph.edu/~iruczins/
! http://bear.fhcrc.org/~clk
! You can contact us at ingo@jhu.edu and clk@fhcrc.org
! =======================================================================


      ! this subroutine performs an annealing step
      ! last modification 06/05/03

      SUBROUTINE annealing(n1,n2,mdl,nkn,nsp,ntr,conc,negs,pick,term,
     #                     storage,slprbc,datri,weight,tstr,tend,tint,
     #                     ehm,msz,nsep,seps,cnc,score,betas,
     #                     ssize,dcph,ordrs,nfcnt,penalty,resp,mtm,
     #                     mcmc,hyperpars,rd1,rd2,rd3,rd4,bout)
      IMPLICIT NONE

        ! parameters
          INTEGER LGCnknMAX,LGCntrMAX,LGCn1MAX,LGCn2MAX,LGCbetaMAX
          PARAMETER (LGCn1MAX   = 10000)
          PARAMETER (LGCn2MAX   =  1000)
          PARAMETER (LGCnknMAX  =   128)
          PARAMETER (LGCntrMAX  =     5)
          PARAMETER (LGCbetaMAX =    55)

        ! arguments in
          INTEGER ehm,n1,n2,mdl,msz,nkn,nsp,ntr,nsep,nfcnt,mtm
          INTEGER cnc(3),dcph(n1),ordrs(n1),storage(2*ntr*nkn*n1) 
          REAL tstr,tend,tint,slprbc(25),weight(n1),rd2(1),rd3(1)
          REAL resp(n1),seps(nsep,n1)
          REAL penalty,hyperpars(10)
          INTEGER datri(n2,n1),mcmc,rd1(1),rd4(1),bout
        ! local
          INTEGER i,j,accept,fcnt,hm,nac,nop,npertemp,tcnt,nrj,nde
          INTEGER nac2,ntot,npckmv(6,LGCntrMAX)
          INTEGER pickmv(6,LGCnknMAX,LGCntrMAX),nsame,nsame2
          REAL acr,scav,sco,temp,rsp(LGCn1MAX),ttemp
          INTEGER iearly(0:6,LGCntrMAX,0:LGCnknMAX,LGCn2MAX,0:1,2)
          INTEGER earlyup
          REAL rearly(0:6,LGCntrMAX,0:LGCnknMAX,LGCn2MAX,0:1,2)
          INTEGER prtr(LGCn1MAX,LGCntrMAX)
          REAL cbetas(0:LGCbetaMAX),xtxsep(LGCbetaMAX+1,LGCbetaMAX+1)
          INTEGER visit(2+ntr*nkn),new
          REAL lvisit(2)

        ! arguments out
          INTEGER conc(nkn,ntr,3)
          INTEGER negs(nkn,ntr,3)
          INTEGER pick(nkn,ntr,3)
          INTEGER term(nkn,ntr,3)
          INTEGER ssize,xstop
          REAL score(3),betas(3,0:(nsep+ntr))
          CHARACTER *80 astring

        ! npckmv: how many available for a move of a particular type 
        ! pickmv: which available for a move of a particular type
        ! prtr: logic trees predictions
        xstop=0
        CALL stopper(LGCn1MAX,n1,"LGCn1MAX","annealing()",8,11,xstop,0)
        CALL stopper(LGCn2MAX,n2,"LGCn2MAX","annealing()",8,11,xstop,0)
        CALL stopper(LGCnknMAX,nkn,"LGCnknMAX","annealing()",9,11,
     #               xstop,0)
        CALL stopper(LGCntrMAX,ntr,"LGCntrMAX","annealing()",9,11,
     #               xstop,0)
        CALL stopper(LGCbetaMAX,nsep+ntr,"LGCbetaMAX","annealing()",10,
     #          11,xstop,1)

      ! calculate response for scoring function
        DO i=1,LGCn1MAX
           DO j=1,LGCntrMAX
              prtr(i,j)=0
           END DO
        END DO
        nac2=0
        DO j=0,(nsep+ntr)
           betas(1,j)=0
           betas(2,j)=0
           betas(3,j)=0
        END DO
        DO j=1,n1
          rsp(j)=resp(j) 
        END DO
      ! carry out an initial simulated annealing step
        CALL annealing_init(n1,n2,mdl,nkn,ntr,conc,negs,pick,term,
     #                      storage,datri,rsp,weight,npckmv,
     #                      pickmv,score,betas,ssize,nsep,seps,
     #                      nop,dcph,ordrs,penalty,prtr,cbetas,xtxsep,
     #                      mtm,mcmc,hyperpars)
        IF (msz.EQ.0) THEN
          IF(ehm.GE.0)THEN
             CALL stringprint("  ",2)
             astring(1:31)='The model of size 0 has score  '
             CALL makerstring(32,43,astring,score(1),7,4)
             CALL stringprint(astring,43)
          END IF
          GOTO 2424
        END IF
        earlyup=0
        CALL clearly(iearly,ntr,nkn,n2)

      ! carry out an annealing chain
        IF (tstr.EQ.0.0.AND.tend.EQ.0.0) THEN
      !   carry out an annealing chain with defaults
      !   find starting temperature
          tstr=3
          acr=1.0 
          tcnt=0
          hm=100
          temp=(10.0**tstr)
          DO WHILE (acr.GT.0.75)
            temp=(10.0**tstr)*(10.0**(-tcnt/5.0))
            nac=0
            nsame=0
            nrj=0
            nde=0
            DO j=1,hm        
              sco=score(1)
              new=nac
              CALL annealing_step(n1,n2,mdl,nkn,nsp,ntr,conc,negs,
     #                            pick,term,storage,npckmv,
     #                            pickmv,temp,slprbc,datri,rsp,
     #                            weight,msz,nsep,seps,cnc,score,
     #                            betas,ssize,nop,dcph,ordrs,accept,
     #                            penalty,iearly,rearly,earlyup,prtr,
     #                            cbetas,xtxsep,mtm,mcmc,hyperpars)
              IF(accept.LT.0)score(1)=sco
              IF(accept.LT.0)nrj=nrj+1
              IF(accept.EQ.0)nde=nde+1
              IF(accept.GT.0)THEN
                earlyup=earlyup+2
                  nac=nac+1
                  IF(sco.EQ.score(1))nsame=nsame+1
              END IF
            END DO
            acr=REAL(nac)/REAL(hm-nrj)
            tcnt=tcnt+1
          END DO

      !   start annealing using the scheme just found
          tstr=REAL(INT(LOG10(temp)))+1+1
          temp=(10.0**tstr)
          CALL writeinfo(ehm,mdl,0,nsep,ntr,temp,score,score(1),betas,
     #         nac,nrj,nde,nsame)
          npertemp=MIN(5000,50*n2)
          if(ehm.GT.0) ehm=npertemp
          tcnt=1
          fcnt=0
          CALL writeinfo(ehm,mdl,0,nsep,ntr,temp,score,score(1),betas,
     #                   nac,nrj,nde,nsame)
      !   if for 5 series in a row there are <10 acceptances stop
          DO WHILE (fcnt.LT.5)
            nac=0
            nsame=0
            nrj=0
            nde=0
            temp=(10.0**tstr)*(10.0**(-tcnt/20.0))
            scav=0.0
            DO j=1,npertemp        
              sco=score(1)
              new=nac
              CALL annealing_step(n1,n2,mdl,nkn,nsp,ntr,conc,negs,
     #                            pick,term,storage,npckmv,
     #                            pickmv,temp,slprbc,datri,rsp,
     #                            weight,msz,nsep,seps,cnc,score,
     #                            betas,ssize,nop,dcph,ordrs,accept,
     #                            penalty,iearly,rearly,earlyup,prtr,
     #                            cbetas,xtxsep,mtm,mcmc,hyperpars)
              IF(accept.LT.0)score(1)=sco
              IF(accept.LT.0)nrj=nrj+1
              IF(accept.GT.0)THEN
                earlyup=earlyup+2
                nac=nac+1
                IF(sco.EQ.score(1))nsame=nsame+1
              END IF
              IF (sco.NE.score(1))THEN
                earlyup=earlyup+2
              END IF
              nde=j-nrj-nac
              scav=(REAL(j-1)*scav+score(1))/REAL(j)
              nac2=nac
              nsame2=nsame
              CALL writeinfo(ehm,mdl,j,nsep,ntr,temp,score,scav,betas,
     #                       nac,nrj,nde,nsame)
              IF (nac.GT.500.AND.j.LT.npertemp) THEN
                CALL writeinfo(ehm,mdl,ehm,nsep,ntr,temp,score,scav,
     #                         betas,nac,nrj,nde,nsame)
                GOTO 1212
              END IF
            END DO
 1212       CONTINUE
            IF ((nac2-nsame2).LE.10) fcnt=fcnt+1
            IF ((nac2-nsame2).GT.10) fcnt=0
            tcnt=tcnt+1
          END DO
          tstr=0.0
          tend=0.0
        ELSE
      !   carry out an annealing chain with pre-specified values
          fcnt=0
          nac=0
          nrj=0
          nde=0
          nac2=0
          nsame=0
          ntot=0
          DO i=0,INT(tint)
            temp=(10.0**tstr)*(10.0**((tend-tstr)/tint))**REAL(i)
            if(mcmc.GT.0)temp= -(REAL(i)/10000)-1
            CALL writeinfo(ehm,mdl,i,nsep,ntr,temp,score,score(1),betas,
     #                     nac,nrj,nde,nsame)
            temp=(10.0**tstr)*(10.0**((tend-tstr)/tint))**REAL(i)
            sco=score(1)
            new=nac
            CALL annealing_step(n1,n2,mdl,nkn,nsp,ntr,conc,negs,pick,
     #                          term,storage,npckmv,pickmv,
     #                          temp,slprbc,datri,rsp,weight,msz,
     #                          nsep,seps,cnc,score,betas,ssize,nop,
     #                          dcph,ordrs,accept,penalty,iearly,
     #                          rearly,earlyup,prtr,cbetas,xtxsep,mtm,
     #                          mcmc,hyperpars)
             IF(accept.LT.0)score(1)=sco
             IF(accept.LT.0)nrj=nrj+1
             IF(accept.GT.0)THEN
               earlyup=earlyup+2
               nac=nac+1
               IF(sco.EQ.score(1))nsame=nsame+1
             END IF
             IF(sco.NE.score(1))THEN
               earlyup=earlyup+2
             END IF
             IF(accept.EQ.0)nde=nde+1
             IF(sco.NE.score(1).and.accept.GT.0) nac2=nac2+1
             ntot=ntot+1
             IF(ntot.EQ.nfcnt.AND.mcmc.EQ.0) THEN
               IF(nac2.LE.10)fcnt=fcnt+1
               IF(nac2.GT.10)fcnt=0
               nac2=0
               ntot=0
               IF(fcnt.EQ.5)goto 1216
             END IF
             IF(mcmc.GT.0.AND.i.GT.nfcnt)THEN
                IF(i.EQ.(nfcnt+1))THEN
                   CALL storeone(mcmc,new,hyperpars,lvisit,visit,
     #                           ntr,nkn,conc,negs,pick,term,-17,
     #                           rd1,rd2,rd3,rd4,bout,n2)
                   IF(bout.GE.0)CALL ciwrite()
                ELSE
                   CALL storeone(mcmc,new,hyperpars,lvisit,visit,
     #                           ntr,nkn,conc,negs,pick,term,nac,
     #                           rd1,rd2,rd3,rd4,bout,n2)
                END IF
             END IF
          END DO
 1216     CONTINUE
        END IF
        IF(nac+nrj+nde.GT.1)THEN
          CALL writeinfo(ehm,mdl,ehm,nsep,ntr,temp,score,score(1),
     #                   betas,nac,nrj,nde,nsame)
        END IF
 2424   CONTINUE
      END 

      ! *****************************************************************
      ! *****************************************************************

      ! this subroutine initializes the first model
      ! last modification 06/05/03

      SUBROUTINE annealing_init(n1,n2,mdl,nkn,ntr,conc,negs,pick,term,
     #                          storage,datri,rsp,weight,npckmv,
     #                          pickmv,score,betas,ssize,nsep,
     #                          seps,nop,dcph,ordrs,penalty,prtr,cbetas,
     #                          xtxsep,mtm,mcmc,hyperpars)
      IMPLICIT NONE

      ! arguments in
        INTEGER n1,n2,mdl,nkn,nsep,ntr,dcph(n1),ordrs(n1),mtm
        INTEGER conc(nkn,ntr,3),mcmc
        INTEGER negs(nkn,ntr,3)
        INTEGER pick(nkn,ntr,3)
        INTEGER term(nkn,ntr,3)
        INTEGER storage(2*ntr*nkn*n1)
        INTEGER prtr(n1,ntr),datri(n2,n1)
        REAL rsp(n1),weight(n1),hyperpars(10)
        REAL seps(nsep,n1),cbetas(0:(nsep+ntr)),xtxsep(0:nsep,0:nsep)
      ! local
        INTEGER j,i,k,wh
        REAL penalty,pp
      ! arguments out
        INTEGER nop,reject,ssize,npckmv(6,ntr)
        INTEGER pickmv(6,nkn,ntr)
        REAL score(3),betas(3,0:(nsep+ntr))
          real tmp
          character *100 astring

      ! precompute separate predictors and intercept for linear regress
        IF(mdl.EQ.2)THEN
          DO i=0,nsep
            DO j=0,nsep
               xtxsep(i,j)=0
            END DO
          END DO
          DO i=1,n1
            xtxsep(0,0)=xtxsep(0,0)+weight(i)
            DO j=1,nsep
              pp=weight(i)*seps(j,i)
              xtxsep(0,j)=xtxsep(0,j)+pp
              DO k=j,nsep
                xtxsep(k,j)=xtxsep(k,j)+pp*seps(k,i)
              END DO
            END DO
          END DO
          DO j=1,nsep
            xtxsep(j,0)=xtxsep(0,j)
            DO k=1,j
               xtxsep(k,j)=xtxsep(j,k)
            END DO
          END DO
        END IF
            

      ! initialize and score the first model
        CALL initialize(n1,ntr,nkn,conc,term,negs,pick,storage,score)
        CALL storing(nkn,ntr,conc,pick,npckmv,pickmv,ssize,nop)
        DO wh=1,ntr
          CALL evaluate_first(wh,n1,n2,nkn,ntr,conc,term,negs,pick,
     #                        datri,prtr)
        END DO
        CALL scoring(prtr,rsp,dcph,ordrs,weight,n1,ntr,mdl,nop,
     #               wh,nsep,seps,score,cbetas,reject,xtxsep,mtm,0)
        DO j=0,(nsep+ntr)
          betas(1,j)=cbetas(j)
        END DO
        IF (reject.EQ.1) THEN
          CALL stringprint("Initial model could not be fitted!",34)
          STOP
        END IF
        IF(mdl.EQ.2)THEN
           score(1)=score(1)+penalty/(REAL(n1))*ssize
        ELSE
           score(1)=score(1)+penalty*ssize
        END IF
        DO j=1,3
          IF(mcmc.EQ.0)score(j)=score(1)
          DO k=0,(nsep+ntr)
            betas(j,k)=betas(1,k)
          END DO
        END DO
        IF(mcmc.GT.0)THEN
           CALL smackonprior(score,ssize,ntr,nkn,conc,term,negs,pick,
     #           hyperpars,n2,-1,weight,(REAL(1.)),0)
        END IF
      END

      ! *****************************************************************
      ! *****************************************************************

      ! this subroutine performs an annealing step
      ! last modification 06/05/03

      SUBROUTINE annealing_step(n1,n2,mdl,nkn,nsp,ntr,conc,negs,pick,
     #                          term,storage,npckmv,pickmv,
     #                          temp,slprbc,datri,rsp,weight,msz,
     #                          nsep,seps,cnc,score,betas,ssize,nop,
     #                          dcph,ordrs,accept,penalty,iearly,rearly,
     #                          earlyup,prtr,cbetas,xtxsep,mtm,mcmc,
     #                          hyperpars)
      IMPLICIT NONE
        ! parameters
          INTEGER LGCnknMAX
          PARAMETER (LGCnknMAX  =   128)
        ! arguments in
          INTEGER n1,n2,msz,mdl,nkn,nsep,nsp,ntr,mtm
          INTEGER cnc(3),mcmc
          INTEGER dcph(n1),ordrs(n1)
          INTEGER npckmv(6,ntr)
          INTEGER pickmv(6,nkn,ntr)
          INTEGER conc(nkn,ntr,3)
          INTEGER negs(nkn,ntr,3)
          INTEGER pick(nkn,ntr,3)
          INTEGER term(nkn,ntr,3)
          INTEGER storage(2*ntr*nkn*n1)
          INTEGER prtr(n1,ntr),datri(n2,n1)
          REAL temp,penalty,hyperpars(10)
          REAL score(3),cbetas(0:(nsep+ntr))
          REAL rsp(n1),weight(n1)
          REAL slprbc(25)
          REAL xtxsep(0:nsep,0:nsep)
          REAL seps(nsep,n1)
          INTEGER iearly(0:6,ntr,0:nkn,n2,0:1,2),earlyup
          REAL rearly(0:6,ntr,0:nkn,n2,0:1,2)
        ! local
          INTEGER j,l1,l2,letter,neg,opper
          INTEGER knt,mtp,nwkv,reject,wh
          INTEGER wkv(LGCnknMAX),xrej,xstop,nopold
          REAL xsc(3)
        ! arguments out
          INTEGER accept,nop,ssize,mctry
          REAL betas(3,0:(nsep+ntr)),rr,rr2,rr3
          

      ! carry out one annealing step
        xstop=0
        hyperpars(8)=0.
        if(mcmc.GT.0)temp=1.
        CALL stopper(LGCnknMAX,nkn,"LGCnknMAX","annealing_step()",9,16,
     #          xstop,1)

        mctry=0
147     CALL moving(n2,nkn,nsp,ntr,conc,negs,pick,term,slprbc,cnc,mcmc,
     #             npckmv,pickmv,msz,ssize,nop,wh,knt,mtp,mctry)
        letter=1
        neg=1
        opper=1
        IF(mtp.EQ.1)THEN
           letter=term(knt,wh,1)
           neg=negs(knt,wh,1)
        END IF
        IF(mtp.EQ.4.or.mtp.EQ.5)THEN
           letter=term(2*knt+1,wh,1)
           neg=negs(2*knt+1,wh,1)
           opper=conc(knt,wh,1)
        END IF

        IF(mtp.GE.3)THEN
           rr=npckmv(mtp,wh)
           rr2=slprbc(1)
           IF(npckmv(2,wh).GT.0)rr2=rr2+slprbc(2)-slprbc(1)
           IF(npckmv(3,wh).GT.0)rr2=rr2+slprbc(3)-slprbc(2)
           IF(npckmv(4,wh).GT.0)rr2=rr2+slprbc(4)-slprbc(3)
           IF(npckmv(5,wh).GT.0)rr2=rr2+slprbc(5)-slprbc(4)
           IF(npckmv(6,wh).GT.0)rr2=rr2+slprbc(6)-slprbc(5)
           rr2=(slprbc(mtp)-slprbc(mtp-1))/rr2
           rr=rr/rr2
         END IF
        IF(iearly(mtp,wh,knt,letter,neg,opper).LE.
     #               earlyup.OR.(mcmc.GT.0).OR.(mtp.eq.0))THEN
        nopold=nop
        CALL storing(nkn,ntr,conc,pick,npckmv,pickmv,ssize,nop)
        CALL evaluating(wh,knt,mtp,n1,n2,nkn,ntr,conc,term,negs,
     #                  datri,prtr,storage,nwkv,wkv)
        CALL scoring(prtr,rsp,dcph,ordrs,weight,n1,ntr,mdl,nop,wh,
     #               nsep,seps,score,cbetas,reject,xtxsep,mtm,nopold)

      ! take care that the there is real convergence
          l1=0 
          l2=0
          IF(betas(1,1).GT. -10000.)l1=1
          IF(betas(1,1).LT. 10000.)l2=1
          IF(l1+l2.EQ.0)reject=1
          DO j=0,(nsep+ntr)
            betas(1,j)=cbetas(j)
          END DO
          IF (reject.EQ.1) THEN
            accept=-1
            iearly(mtp,wh,knt,letter,neg,opper)=earlyup+2
          ELSE 
            IF(mdl.EQ.2)THEN
               score(1)=score(1)+penalty/(REAL(n1))*ssize
            ELSE 
               score(1)=score(1)+penalty*ssize
            END IF
            iearly(mtp,wh,knt,letter,neg,opper)=earlyup+1
            rearly(mtp,wh,knt,letter,neg,opper)=score(1)
            IF(mcmc.GT.0) THEN
               IF(mtp.GE.3)THEN
                  IF(mtp.EQ.3)rr=rr/(REAL(npckmv(4,wh)))
                  IF(mtp.EQ.4)rr=rr/(REAL(npckmv(3,wh)))
                  IF(mtp.EQ.5)rr=rr/(REAL(npckmv(6,wh)))
                  IF(mtp.EQ.6)rr=rr/(REAL(npckmv(5,wh)))
                  rr3=slprbc(1)
                  IF(npckmv(2,wh).GT.0)rr3=rr3+slprbc(2)-slprbc(1)
                  IF(npckmv(3,wh).GT.0)rr3=rr3+slprbc(3)-slprbc(2)
                  IF(npckmv(4,wh).GT.0)rr3=rr3+slprbc(4)-slprbc(3)
                  IF(npckmv(5,wh).GT.0)rr3=rr3+slprbc(5)-slprbc(4)
                  IF(npckmv(6,wh).GT.0)rr3=rr3+slprbc(6)-slprbc(5)
                  IF(mtp.EQ.3)rr3=(slprbc(4)-slprbc(3))/rr3
                  IF(mtp.EQ.4)rr3=(slprbc(3)-slprbc(2))/rr3
                  IF(mtp.EQ.5)rr3=(slprbc(6)-slprbc(5))/rr3
                  IF(mtp.EQ.6)rr3=(slprbc(5)-slprbc(4))/rr3
                  rr=rr*rr3
               ELSE
                  rr=1.
               END IF
               CALL smackonprior(score,ssize,ntr,nkn,conc,term,negs,
     #           pick,hyperpars,n2,mtp,slprbc,rr,nopold-nop)
            END IF
            CALL deciding(score,temp,accept,hyperpars(8),mcmc)
          END IF
          CALL recording(accept,wh,nkn,ntr,nsep,score,betas,
     #                   conc,negs,pick,term,mcmc)
          CALL restoring(accept,wh,n1,nkn,ntr,storage,nwkv,wkv)
          CALL storing(nkn,ntr,conc,pick,npckmv,pickmv,ssize,nop)
          IF(reject.EQ.1.AND.mcmc.GE.0.AND.mctry.LT.10)THEN
             mctry=mctry+1
             GOTO 147
          END IF
        ELSE
          IF(iearly(mtp,wh,knt,letter,neg,opper).EQ.(earlyup+2))THEN
            reject=1
            accept=-1
          ELSE
            reject=0
            score(1)=rearly(mtp,wh,knt,letter,neg,opper)
            CALL deciding(score,temp,accept,hyperpars(8),mcmc)
            IF(accept.EQ.1)THEN
              CALL storing(nkn,ntr,conc,pick,npckmv,pickmv,ssize,nop)
              CALL evaluating(wh,knt,mtp,n1,n2,nkn,ntr,conc,term,negs,
     #                        datri,prtr,storage,nwkv,wkv)
              CALL scoring(prtr,rsp,dcph,ordrs,weight,n1,ntr,mdl,nop,wh,
     #                     nsep,seps,xsc,cbetas,xrej,xtxsep,mtm,1000)
              CALL recording(accept,wh,nkn,ntr,nsep,score,betas,
     #                       conc,negs,pick,term,mcmc)
              CALL restoring(accept,wh,n1,nkn,ntr,storage,nwkv,wkv)
              CALL storing(nkn,ntr,conc,pick,npckmv,pickmv,ssize,nop)
            ELSE
              score(1)=score(2)
            END IF
          END IF
        END IF
           
      END 

      ! *****************************************************************
      ! *****************************************************************

      ! this subroutine creates the output during the annealing run
      ! last modification 10/16/02

      SUBROUTINE writeinfo(ehm,mdl,i,nsep,ntr,temp,score,scav,betas,
     #                   nac,nrj,nde,nsame)
      IMPLICIT NONE

        ! arguments in
          INTEGER ehm,i,mdl,nsep,ntr,nac,nrj,nde,nsame
          REAL scav,temp,score(3),betas(3,0:(ntr+nsep)),t2
          DOUBLE PRECISION t1,mylog
        ! local
          INTEGER j
          CHARACTER *120 astring
        ! arguments out
        ! none

      ! print info
        t1=temp
        IF(t1.GT.0.0)  THEN
            t1=mylog(t1)/2.302585092994046
        ELSE 
            t1= -temp-1 
        END IF
        t2=t1
        IF (ehm.GT.0) THEN
          IF (mdl.EQ.1) THEN
            IF (i.EQ.0) THEN
              CALL stringprint('  ',2)
              astring(1:26)="log-temp current score    "
              astring(27:59)="best score        acc / rej /sing"
              CALL stringprint(astring,59)
            END IF
            IF (MOD(i,ehm).EQ.0) THEN
              CALL makerstring(1,8,astring,t2,4,3)
              CALL makerstring(9,22,astring,score(1),9,4)
              CALL makerstring(23,36,astring,score(3),9,4)
              CALL makeistring(37,42,astring,nac-nsame,6)
              astring(43:43)="("
              CALL makeistring(44,46,astring,nsame,3)
              astring(47:47)=")"
              CALL makeistring(48,53,astring,nde,6)
              CALL makeistring(54,59,astring,nrj,6)
              CALL stringprint(astring,59)
            END IF
          ELSE 
            IF (i.EQ.0) THEN
              CALL stringprint('  ',2)
              astring(1:26)="log-temp current score    "
              IF(temp.LT.0)astring(1:26)="iter(10k)  current scr    "
              astring(27:59)="best score        acc / rej /sing"
              astring(60:81)="    current parameters"
              CALL stringprint(astring,81)
            END IF
            IF (MOD(i,ehm).EQ.0) THEN
              CALL makerstring(1,8,astring,t2,4,3)
              CALL makerstring(9,22,astring,score(1),9,4)
              CALL makerstring(23,36,astring,score(3),9,4)
              IF(i.eq.0.AND.temp.LT.0)THEN
              CALL makerstring(23,36,astring,score(1),9,4)
              END IF
              CALL makeistring(37,42,astring,nac-nsame,6)
              astring(43:43)="("
              CALL makeistring(44,46,astring,nsame,3)
              astring(47:47)=")"
              CALL makeistring(48,53,astring,nde,6)
              CALL makeistring(54,59,astring,nrj,6)
              IF((nsep+ntr).LT.4)THEN
                DO j=0,(nsep+ntr)
                  astring((60+j*9):(60+j*9))=" "
                  CALL makerstring(61+j*9,68+j*9,astring,
     #                               betas(1,j),4,3)
                END DO
                CALL stringprint(astring,68+9*(nsep+ntr))
              ELSE
                DO j=0,4
                  astring((60+j*9):(60+j*9))=" "
                  CALL makerstring(61+j*9,68+j*9,astring,
     #                               betas(1,j),4,3)
                END DO
                CALL stringprint(astring,68+9*4)
                astring(1:12)="===========>"
                DO j=5,(nsep+ntr)
                  IF(j.LT.13)THEN
                    astring((j*9-26):(j+9-26))=" "
                    CALL makerstring(j*9-13,j*9-6,astring,
     #                               betas(1,j),4,3)
                  END IF
                  IF(j.EQ.13) astring(85:89)="...."
                END DO
                IF(nsep+ntr.GE.13)THEN
                   CALL stringprint(astring,89)
                ELSE
                   IF(nsep+ntr.GT.4)
     #               CALL stringprint(astring,(nsep+ntr-2)*9-5)
                END IF
              END IF
            END IF
          END IF
          IF (MOD(i,ehm).EQ.0) THEN
            nac=0
            nde=0
            nrj=0
            nsame=0
          END IF
        END IF

      END 

      ! *****************************************************************
      ! *****************************************************************

      ! this subroutine does the crossvalidation 
      ! last modification 10/16/02

      SUBROUTINE crossvalx(kfold,n1,n2,mdl,
     #                    nkn,nsp,ntr,conc,negs,pick,term,storage,
     #                    slprbc,datri,weight,tstr,tend,tint,ehm,msz,
     #                    nsep,seps,cnc,score,betas,ssize,dcph,ordrs,
     #                    nfcnt,seed,resp,rnumsrr,mtm,ltree,ioscores)
      IMPLICIT NONE

        ! parameters
          INTEGER LGCnsepMAX,LGCn1MAX,LGCn2MAX
          PARAMETER (LGCn1MAX   = 10000)
          PARAMETER (LGCn2MAX  =   1000)
          PARAMETER (LGCnsepMAX =    50)
        ! arguments in
          INTEGER ehm,kfold,n1,n2,mdl,msz,nkn,nsp,ntr,nsep,nfcnt,cnc(3)
          INTEGER dcph(n1),ordrs(n1),seed
          REAL tstr,tend,tint,slprbc(25),weight(n1)
          INTEGER datri(n2,n1) ,mtm
          REAL seps(nsep,n1),resp(n1),rnumsrr(n1)
        ! arguments out
          INTEGER conc(nkn,ntr,3),ltree
          INTEGER negs(nkn,ntr,3)
          INTEGER pick(nkn,ntr,3)
          INTEGER term(nkn,ntr,3)
          INTEGER ssize,storage(2*ntr*nkn*n1)
          REAL score(3),betas(3,0:(nsep+ntr)),ioscores(1)
        ! local
          REAL sepstt(LGCnsepMAX,LGCn1MAX)
          INTEGER datritt(LGCn2MAX,LGCn1MAX)
          INTEGER xstop
        xstop=0
        CALL stopper(LGCn1MAX,n1,"LGCn1MAX","crossvalx()",8,11,xstop,0)
        CALL stopper(LGCn2MAX,n2,"LGCn2MAX","crossvalx()",8,11,xstop,0)
        CALL stopper(LGCnsepMAX,nsep,"LGCnsepMAX","crossvalx()",10,11,
     #          xstop,1)

        CALL crossval(kfold,n1,n2,mdl,
     #                    nkn,nsp,ntr,conc,negs,pick,term,storage,
     #                    slprbc,datri,weight,tstr,tend,tint,ehm,msz,
     #                    nsep,seps,cnc,score,betas,ssize,dcph,ordrs,
     #                    nfcnt,seed,sepstt,datritt,resp,rnumsrr,mtm,
     #                    ltree,ioscores)
        END


      ! *****************************************************************
      ! *****************************************************************

      ! this subroutine does the crossvalidation 
      ! last modification 06/05/03

      SUBROUTINE crossval(kfold,n1,n2,mdl,
     #                    nkn,nsp,ntr,conc,negs,pick,term,storage,
     #                    slprbc,datri,weight,tstr,tend,tint,ehm,msz,
     #                    nsep,seps,cnc,score,betas,ssize,dcph,ordrs,
     #                    nfcnt,seed,sepstt,datritt,resp,rnumsrr,mtm,
     #                    ltree,ioscores)
      IMPLICIT NONE

        ! parameters
          INTEGER LGCn1MAX
          PARAMETER (LGCn1MAX   = 10000)
        ! arguments in
          INTEGER ehm,kfold,n1,n2,mdl,msz,nkn,nsp,ntr,nsep,nfcnt,cnc(3)
          INTEGER dcph(n1),ordrs(n1),seed,mtm
          REAL tstr,tend,tint,slprbc(25),weight(n1)
          INTEGER datri(n2,n1) 
          REAL seps(nsep,n1),resp(n1),rnumsrr(n1)
        ! arguments out
          INTEGER conc(nkn,ntr,3),ltree
          INTEGER negs(nkn,ntr,3)
          INTEGER pick(nkn,ntr,3)
          INTEGER term(nkn,ntr,3)
          INTEGER ssize,storage(2*ntr*nkn*n1)
          REAL score(3),betas(3,0:(nsep+ntr)),ioscores(1)
        ! local
          INTEGER hmpg,hmlo,ngrp,hmpgx,ngrpx
          INTEGER i,j,k,l,cnt,cnt2,dummy(LGCn1MAX),rnumsi(LGCn1MAX)
          INTEGER ngrphere,hmpghere
          REAL cvscore(4),rnumsr(LGCn1MAX)
          REAL weighttt(LGCn1MAX)
          INTEGER datritt(n2,n1)
          INTEGER dcphtt(LGCn1MAX),mcmc
          INTEGER ordrstt(LGCn1MAX),xstop,rd1(2),rd4(2),bout
          REAL datrtokentt(LGCn1MAX),resptt(LGCn1MAX)
          REAL sepstt(nsep,n1),penalty,rdummy(10),rd2(2),rd3(2)
          CHARACTER *125 astring
          xstop=0
          CALL stopper(LGCn1MAX,n1,"LGCn1MAX","crossval()",8,10,xstop,1)

        ! cvscore: cross-validation scores (tr,travg,test,testavg)

          mcmc=0
      ! randomize all other cases and allocate training and test data
          hmpg=INT(REAL(n1)/REAL(kfold))
          hmlo=n1-kfold*hmpg
          ngrp=n1-hmpg
      ! nondividor
          IF(hmlo.GT.0)THEN
            hmpgx=hmpg+1
            ngrpx=ngrp-1
          ELSE
            hmpgx=hmpg
            ngrpx=ngrp
          END IF
          DO j=1,n1
            dummy(j)=j
            rnumsr(j)=rnumsrr(j)
          END DO
          l=0
          DO j=1,kfold
            DO k=1,hmpg
              l=l+1
              rnumsi(l)=j
            END DO
            IF(j.LE.hmlo)THEN
              l=l+1
              rnumsi(l)=j
            END IF
          END DO
          IF(seed.GE.0)CALL psort2(rnumsr,n1,dummy,n1,rnumsi)
      ! cross validation steps
          DO k=1,4
            cvscore(k)=0.0
          END DO
          penalty=0
          IF(ehm.EQ.0)THEN
             DO k=1,48
                astring(k:k)=" "
             END DO
             astring(49:61)=" training-now"
             astring(62:74)=" training-ave"
             astring(75:87)="     test-now"
             astring(88:100)="     test-ave"
             CALL stringprint(astring,100)
          END IF
          
          DO k=1,kfold
            cnt=0
            cnt2=0
            DO j=1,n1
              IF(rnumsi(j).NE.k)THEN
                cnt=cnt+1
                weighttt(cnt)=weight(j)
                dcphtt(cnt)=dcph(j)
                ordrstt(cnt)=cnt
                IF (mdl.EQ.4.or.mdl.EQ.5) datrtokentt(cnt)=resp(j)
                DO i=1,n2
                  datritt(i,cnt)=datri(i,j)
                END DO
                resptt(cnt)=resp(j)
                IF (nsep.GT.0) THEN
                  DO l=1,nsep
                    sepstt(l,cnt)=seps(l,j)
                  END DO
                END IF
              END IF
            END DO
            IF(k.GT.hmlo)THEN
              ngrphere=ngrp
              hmpghere=hmpg
            ELSE
              ngrphere=ngrpx
              hmpghere=hmpgx
            END IF

            IF (mdl.EQ.4.or.mdl.EQ.5) THEN
              CALL psort2(datrtokentt,ngrphere,dummy,ngrphere,ordrstt)
            END IF
            CALL annealing(ngrphere,n2,mdl,nkn,nsp,ntr,conc,negs,pick,
     #                     term,storage,slprbc,datritt,weighttt,tstr,
     #                     tend,tint,ehm,msz,nsep,sepstt,cnc,score,
     #                     betas,ssize,dcphtt,ordrstt,nfcnt,penalty,
     #                     resptt,mtm,mcmc,rdummy,rd1,rd2,rd3,rd4,bout)
            cvscore(1)=score(3)
            cvscore(2)=(REAL(k-1)*cvscore(2)+cvscore(1))/REAL(k)
            DO j=1,n1
              IF(rnumsi(j).EQ.k)THEN
                cnt2=cnt2+1
                weighttt(cnt2)=weight(j)
                dcphtt(cnt2)=dcph(j)
                ordrstt(cnt2)=cnt2
                IF (mdl.EQ.4.or.mdl.EQ.5) datrtokentt(cnt2)=resp(j)
                DO i=1,n2
                  datritt(i,cnt2)=datri(i,j)
                END DO
                resptt(cnt2)=resp(j)
                IF (nsep.GT.0) THEN
                  DO l=1,nsep
                    sepstt(l,cnt2)=seps(l,j)
                  END DO
                END IF
              END IF
            END DO
            IF (mdl.EQ.4.or.mdl.EQ.5) THEN
              CALL psort2(datrtokentt,hmpghere,dummy,hmpghere,ordrstt)
            END IF
            CALL testsetx(hmpghere,n2,mdl,nkn,ntr,conc,negs,pick,term,
     #                   betas,datritt,weighttt,dcphtt,ordrstt,
     #                   nsep,sepstt,score,resptt)
            cvscore(3)=score(1)
            cvscore(4)=(REAL(k-1)*cvscore(4)+cvscore(3))/REAL(k)

            IF(ehm.GE.0)THEN
               IF(ehm.GT.0)THEN
                  DO l=1,48
                     astring(l:l)=" "
                  END DO
                  astring(49:61)=" training-now"
                  astring(62:74)=" training-ave"
                  astring(75:87)="     test-now"
                  astring(88:100)="     test-ave"
                  CALL stringprint(astring,100)
               END IF
               astring(1:4)="Step"
               CALL makeistring(5,7,astring,k,3)
               astring(8:10)=" of"
               CALL makeistring(11,13,astring,kfold,3)
               astring(14:15)=" ["
               CALL makeistring(16,18,astring,ntr,3)
               astring(19:27)=" trees;"
               CALL makeistring(27,29,astring,msz,3)
               astring(30:49)=" leaves] CV score: "
               CALL makerstring(50,62,astring,cvscore(1),8,4)
               CALL makerstring(63,75,astring,cvscore(2),8,4)
               CALL makerstring(76,88,astring,cvscore(3),8,4)
               CALL makerstring(89,101,astring,cvscore(4),8,4)
               CALL stringprint(astring,101)
            END IF
            ioscores(ltree*8+1)=ntr
            ioscores(ltree*8+2)=msz
            ioscores(ltree*8+3)=k
            ioscores(ltree*8+4)=kfold
            ioscores(ltree*8+5)=cvscore(1)
            ioscores(ltree*8+6)=cvscore(2)
            ioscores(ltree*8+7)=cvscore(3)
            ioscores(ltree*8+8)=cvscore(4)
            ltree=ltree+1
          END DO

      END 

      ! *****************************************************************
      ! *****************************************************************


      ! this subroutine makes the decision whether or not to accept the move
      ! last modification 06/05/03

      SUBROUTINE deciding(score,temp,accept,postrat,mcmc)
      IMPLICIT NONE

        ! arguments in
          REAL temp,score(3),postrat
        ! local
          REAL rnum,myrand
          DOUBLE PRECISION crit,z,myexp
        ! arguments out
          INTEGER accept,mcmc

      ! make decision using standard acceptance function
        rnum=myrand(0)
        accept=0
        IF(mcmc.EQ.0)THEN
           z=-(score(1)-score(2))/temp
        ELSE
      ! Bayesian
           z= -score(1)+score(2)+postrat
        END IF
        IF(z.GT.0)THEN
          accept=1
        ELSE
          crit=myexp(z)
          IF (rnum.LT.crit) accept=1
        END IF

      END 

      ! *****************************************************************
      ! *****************************************************************

      ! this subroutine evaluates the initial (empty) logic trees.
      ! to start with non-empty trees, use commented out statement
      ! in the last line for some more programming.
      ! it is rather inefficient, but only used in the very first step
      ! of the algorithm, so it doesn't matter.
      ! last modification 10/16/02

      SUBROUTINE evaluate_first(wh,n1,n2,nkn,ntr,conc,term,negs,pick,
     #                          datri,prtr)
      IMPLICIT NONE

        ! parameters
          INTEGER LGCnknMAX,LGCn1MAX
          PARAMETER (LGCn1MAX   = 10000)
          PARAMETER (LGCnknMAX  =   128)
        ! arguments in
          INTEGER n1,n2,nkn,ntr,wh
          INTEGER conc(nkn,ntr,3)
          INTEGER negs(nkn,ntr,3)
          INTEGER pick(nkn,ntr,3)
          INTEGER term(nkn,ntr,3)
          INTEGER datri(n2,n1)
        ! local
          INTEGER dim,dp,i,j,k,pc,shuffle,rank(LGCnknMAX),iterm
          INTEGER mat(LGCn1MAX+2,LGCnknMAX),xstop
        ! arguments out
          INTEGER prtr(n1,ntr)
          xstop=0
          CALL stopper(LGCn1MAX,n1,"LGCn1MAX","evaluate_first()",8,16,
     #            xstop,0)
          CALL stopper(LGCnknMAX,nkn,"LGCnknMAX","evaluate_first()",9,
     #            16,xstop,1)


        ! initialize mat and rank
          DO k=1,nkn
            rank(k)=0
            DO j=1,(n1+2)
              mat(j,k)=0
            END DO
          END DO

      ! write leaf information into mat
        pc=1
        DO i=1,nkn
          IF (pick(i,wh,1).EQ.1) THEN
            IF (conc(i,wh,1).EQ.3) THEN
             iterm=term(i,wh,1)
              IF (negs(i,wh,1).EQ.0) THEN
                DO dp=1,n1
                  mat(dp,pc)=datri(iterm,dp)           
                END DO
              ELSE
                DO dp=1,n1
                  mat(dp,pc)=(1-datri(iterm,dp))**2  
                END DO
              END IF
              mat(n1+1,pc)=i
              IF (i.GT.1) mat(n1+2,pc)=conc(INT(i/2),wh,1)
              pc=pc+1
            END IF
          END IF
        END DO
        dim=pc-1

      ! evaluate mat sequentially in the second argument as good as possible
 1000   CONTINUE
        IF (dim.GT.1) THEN
          DO pc=dim,2,-1
            IF (((MOD(mat(n1+1,pc),2).EQ.1).AND.
     #           (mat(n1+1,pc).EQ.mat(n1+1,pc-1)+1)).OR.
     #         ((MOD(mat(n1+1,pc),2).EQ.0).AND.
     #           (mat(n1+1,pc).EQ.mat(n1+1,pc-1)-1))) THEN
              IF (mat(n1+2,pc).EQ.1) THEN
                DO j=1,n1
                  mat(j,pc-1)=mat(j,pc-1)*mat(j,pc)
                END DO
              ELSE
                DO j=1,n1
                  mat(j,pc-1)= MAX(mat(j,pc-1),mat(j,pc))
                END DO
              END IF
              mat(n1+1,pc-1)=INT(mat(n1+1,pc)/2)
              IF (INT(mat(n1+1,pc)/2).GT.1) THEN
                 mat(n1+2,pc-1)=conc(INT(INT(mat(n1+1,pc)/2)/2),wh,1)
              ELSE
                mat(n1+2,pc-1)=0
              END IF
              DO k=pc,dim
                DO j=1,n1+2
                  mat(j,k)=mat(j,k+1)
                END DO
              END DO
              dim=dim-1
              GOTO 1000
            END IF
          END DO
        END IF

      ! shuffle mat in second argument if necessary
        IF (dim.GT.1) THEN
          shuffle=0
 2000     CONTINUE
          DO WHILE (shuffle.LT.dim)
            shuffle=shuffle+1
            DO i=1,dim
              rank(i)=1
              DO j=1,dim
                IF (mat(n1+1,i).GT.mat(n1+1,j)) rank(i)=rank(i)+1
              END DO
            END DO 
            DO i=1,dim
              DO j=1,(dim-1)
                IF ((rank(j).EQ.i).AND.(i.NE.j)) THEN
                  DO k=1,n1+2
                    mat(k,dim+1)=mat(k,i)
                    mat(k,i)=mat(k,j)
                    mat(k,j)=mat(k,dim+1)                    
                    mat(k,dim+1)=0
                  END DO
                  GOTO 2000
                END IF 
              END DO
            END DO 
          END DO
          GOTO 1000
        END IF         

      ! write out tree prediction
        DO j=1,n1
          prtr(j,wh)=mat(j,1)
        END DO

      END 

      ! *****************************************************************
      ! *****************************************************************


      ! this subroutine chooses the appropriate evaluation routine
      ! last modification 06/05/03

      SUBROUTINE evaluating(wh,knt,mtp,n1,n2,nkn,ntr,conc,term,negs,
     #                      datri,prtr,storage,nwkv,wkv)
      IMPLICIT NONE

        ! arguments in
          INTEGER knt,mtp,n1,n2,nkn,ntr,wh
          INTEGER conc(nkn,ntr,3)
          INTEGER negs(nkn,ntr,3)
          INTEGER term(nkn,ntr,3)
          INTEGER datri(n2,n1)
        ! local
          INTEGER j,k
        ! arguments out
      ! choose routine
          INTEGER nwkv,wkv(nkn),prtr(n1,ntr)
          INTEGER storage(2*ntr*nkn*n1)

        IF (mtp.EQ.0) THEN
          CALL evaluate_firstknot(wh,n1,n2,nkn,ntr,term,negs,
     #                            datri,storage,nwkv,wkv)
        ELSE IF (mtp.EQ.1) THEN
          CALL evaluate_altlf(wh,knt,n1,n2,nkn,ntr,conc,term,negs,
     #                        datri,storage,nwkv,wkv)
        ELSE IF (mtp.EQ.2) THEN
          CALL evaluate_altop(wh,knt,n1,nkn,ntr,conc,term,negs,
     #                        storage,nwkv,wkv)
        ELSE IF (mtp.EQ.3) THEN
          CALL evaluate_delete(wh,knt,n1,nkn,ntr,conc,term,negs,
     #                         storage,nwkv,wkv)
        ELSE IF (mtp.EQ.4) THEN
          CALL evaluate_split(wh,knt,n1,n2,nkn,ntr,conc,term,negs,
     #                        datri,storage,nwkv,wkv)
        ELSE IF (mtp.EQ.5) THEN
          CALL evaluate_branch(wh,knt,n1,n2,nkn,ntr,conc,term,negs,
     #                         datri,storage,nwkv,wkv)
        ELSE IF (mtp.EQ.6) THEN
          CALL evaluate_prune(wh,knt,n1,n2,nkn,ntr,conc,term,negs,
     #                        datri,storage,nwkv,wkv)
        END IF

        DO j=1,ntr
          DO k=1,n1
            prtr(k,j)=storage(k+n1*(j-1)*nkn)
          END DO
        END DO

      END 

      ! *****************************************************************
      ! *****************************************************************

      ! this subroutine evaluates the tree that has been started
      ! last modification 10/16/02

      SUBROUTINE evaluate_firstknot(wh,n1,n2,nkn,ntr,term,negs,
     #                              datri,storage,nwkv,wkv)
      IMPLICIT NONE

        ! arguments in
          INTEGER n1,n2,nkn,ntr,wh,negs(nkn,ntr,3)
          INTEGER term(nkn,ntr,3)
          INTEGER datri(n2,n1)
        ! local
          INTEGER j,iterm,pp
        ! arguments out
          INTEGER nwkv,wkv(nkn)
          INTEGER storage(2*ntr*nkn*n1)

      ! wkv stands for "Which Knots Visited" (used in restoring routine)
        nwkv=0
        DO j=1,nkn
          wkv(j)=0
        END DO

      ! store predictor in first knot
        iterm=term(1,wh,1)
        pp=(wh-1)*nkn*n1
        IF (negs(1,wh,1).EQ.0) THEN
          DO j=1,n1
            storage(j+pp)=datri(iterm,j)
          END DO
        ELSE 
          DO j=1,n1
            storage(j+pp)=1-datri(iterm,j)
          END DO
        END IF
        nwkv=nwkv+1
        wkv(nwkv)=1

      END

      ! *****************************************************************
      ! *****************************************************************



      ! this subroutine evaluates the tree after a leaf was changed
      ! last modification 10/16/02

      SUBROUTINE evaluate_altlf(wh,knt,n1,n2,nkn,ntr,conc,term,negs,
     #                          datri,storage,nwkv,wkv)
      IMPLICIT NONE

        ! arguments in
          INTEGER knt,n1,n2,nkn,ntr,wh
          INTEGER conc(nkn,ntr,3)
          INTEGER negs(nkn,ntr,3)
          INTEGER term(nkn,ntr,3)
          INTEGER datri(n2,n1)
        ! local
          INTEGER j,knac,iterm,knac2,pp,pp2,pp2b
        ! arguments out
          INTEGER nwkv,wkv(nkn)
          INTEGER storage(2*ntr*nkn*n1)

      ! wkv stands for "Which Knots Visited" (used in restoring routine)
        nwkv=0
        DO j=1,nkn
          wkv(j)=0
        END DO

      ! store predictor in knt
        knac=knt
        iterm=term(knac,wh,1)
        pp=(wh-1)*nkn*n1+(knac-1)*n1
        IF (negs(knac,wh,1).EQ.0) THEN
          DO j=1,n1
            storage(pp+j)=datri(iterm,j)
          END DO
        ELSE 
          DO j=1,n1
            storage(pp+j)=1-datri(iterm,j)
          END DO
        END IF
        nwkv=nwkv+1
        wkv(nwkv)=knac
        knac=INT(REAL(knac)/2.0)

      ! calculate predictions in storage
        DO WHILE (knac.GT.0)
          knac2=2*knac
          pp=(wh-1)*nkn*n1+(knac-1)*n1
          pp2=(wh-1)*nkn*n1+(knac2-1)*n1
          pp2b=(wh-1)*nkn*n1+(knac2)*n1
          IF (conc(knac,wh,1).EQ.1) THEN
            DO j=1,n1
              storage(pp+j)=storage(pp2+j)*storage(pp2b+j)
            END DO
          ELSE  
            DO j=1,n1
              storage(pp+j)=1-(1-storage(pp2+j))*(1-storage(pp2b+j))
            END DO
          END IF
          nwkv=nwkv+1
          wkv(nwkv)=knac
          knac=INT(REAL(knac)/2.0)
        ENDDO

      END

      ! *****************************************************************
      ! *****************************************************************

      ! this subroutine evaluates the tree after an operator was changed
      ! last modification 10/16/02

      SUBROUTINE evaluate_altop(wh,knt,n1,nkn,ntr,conc,term,negs,
     #                          storage,nwkv,wkv)
      IMPLICIT NONE

        ! arguments in
          INTEGER knt,n1,nkn,ntr,wh
          INTEGER conc(nkn,ntr,3)
          INTEGER negs(nkn,ntr,3)
          INTEGER term(nkn,ntr,3)
        ! local
          INTEGER j,knac,knac2,pp,pp2,pp2b
        ! arguments out
          INTEGER nwkv,wkv(nkn)
          INTEGER storage(2*ntr*nkn*n1)

      ! wkv stands for "Which Knots Visited" (used in restoring routine)
        nwkv=0
        DO j=1,nkn
          wkv(j)=0
        END DO

      ! calculate predictions in storage
        knac=knt
        DO WHILE (knac.GT.0)
          knac2=2*knac
          pp=(wh-1)*nkn*n1+(knac-1)*n1
          pp2=(wh-1)*nkn*n1+(knac2-1)*n1
          pp2b=(wh-1)*nkn*n1+(knac2)*n1
          IF (conc(knac,wh,1).EQ.1) THEN
            DO j=1,n1
              storage(pp+j)=storage(pp2+j)*storage(pp2b+j)
            END DO
          ELSE  
            DO j=1,n1
              storage(pp+j)=1-(1-storage(pp2+j))*(1-storage(pp2b+j))
            END DO
          END IF
          nwkv=nwkv+1
          wkv(nwkv)=knac
          knac=INT(REAL(knac)/2.0)
        ENDDO

      END

      ! *****************************************************************
      ! *****************************************************************

      ! this subroutine evaluates the tree after a leaf was deleted
      ! last modification 10/16/02

      SUBROUTINE evaluate_delete(wh,knt,n1,nkn,ntr,conc,term,negs,
     #                           storage,nwkv,wkv)
      IMPLICIT NONE

        ! arguments in
          INTEGER knt,n1,nkn,ntr,wh
          INTEGER conc(nkn,ntr,3)
          INTEGER negs(nkn,ntr,3)
          INTEGER term(nkn,ntr,3)
        ! local
          INTEGER j,knac,sibling,knac2,pp,pp2,pp2b
        ! arguments out
          INTEGER nwkv,wkv(nkn)
          INTEGER storage(2*ntr*nkn*n1)

      ! wkv stands for "Which Knots Visited" (used in restoring routine)
        nwkv=0
        DO j=1,nkn
          wkv(j)=0
        END DO

      ! calculate predictions in storage
        IF (knt.EQ.1) THEN
          pp=(wh-1)*nkn*n1
          DO j=1,n1
            storage(pp+j)=0
          END DO
          nwkv=1
          wkv(nwkv)=1
        ELSE 
          knac=INT(REAL(knt)/2.0)
          IF (MOD(knt,2).EQ.0) THEN
            sibling=knt+1
          ELSE
            sibling=knt-1
          END IF
          pp=(wh-1)*nkn*n1+(knac-1)*n1
          pp2=(wh-1)*nkn*n1+(sibling-1)*n1
          DO j=1,n1
            storage(pp+j)=storage(pp2+j)
          END DO
          nwkv=nwkv+1
          wkv(nwkv)=knac
          IF (knac.GT.1) THEN
            knac=INT(REAL(knac)/2.0)
            DO WHILE (knac.GT.0)
              knac2=2*knac
              pp=(wh-1)*nkn*n1+(knac-1)*n1
              pp2=(wh-1)*nkn*n1+(knac2-1)*n1
              pp2b=(wh-1)*nkn*n1+(knac2)*n1
              IF (conc(knac,wh,1).EQ.1) THEN
                DO j=1,n1
                  storage(pp+j)=storage(pp2+j)*storage(pp2b+j)
                END DO
              ELSE  
                DO j=1,n1
                  storage(pp+j)=1-(1-storage(pp2+j))*(1-storage(pp2b+j))
                END DO
              END IF
              nwkv=nwkv+1
              wkv(nwkv)=knac
              knac=INT(REAL(knac)/2.0)
            ENDDO
          ENDIF 
        ENDIF

      END 

      ! *****************************************************************
      ! *****************************************************************

      ! this subroutine evaluates the tree after a leaf was split
      ! last modification 10/16/02

      SUBROUTINE evaluate_split(wh,knt,n1,n2,nkn,ntr,conc,term,negs,
     #                          datri,storage,nwkv,wkv)
      IMPLICIT NONE

        ! arguments in
          INTEGER knt,n1,n2,nkn,ntr,wh
          INTEGER conc(nkn,ntr,3)
          INTEGER negs(nkn,ntr,3)
          INTEGER term(nkn,ntr,3)
          INTEGER datri(n2,n1)
        ! local
          INTEGER j,knac,iterm,knac2,pp,pp2,pp2b
        ! arguments out
          INTEGER nwkv,wkv(nkn)
          INTEGER storage(2*ntr*nkn*n1)

      ! wkv stands for "Which Knots Visited" (used in restoring routine)
        nwkv=0
        DO j=1,nkn
          wkv(j)=0
        END DO

      ! store predictor in 2*knt and 2*knt+1
        knac=2*knt
        iterm=term(knac,wh,1)
        pp=(wh-1)*nkn*n1+(knac-1)*n1
        IF (negs(knac,wh,1).EQ.0) THEN
          DO j=1,n1
            storage(pp+j)=datri(iterm,j)
          END DO
        ELSE 
          DO j=1,n1
            storage(pp+j)=1-datri(iterm,j)
          END DO
        END IF
        nwkv=nwkv+1
        wkv(nwkv)=knac
        iterm=term(knac+1,wh,1)
        pp=(wh-1)*nkn*n1+knac*n1
        IF (negs(knac+1,wh,1).EQ.0) THEN
          DO j=1,n1
            storage(pp+j)=datri(iterm,j)
          END DO
        ELSE 
          DO j=1,n1
            storage(pp+j)=1-datri(iterm,j)
          END DO
        END IF
        nwkv=nwkv+1
        wkv(nwkv)=knac+1
        knac=INT(REAL(knac)/2.0)

      ! calculate predictions in storage
        DO WHILE (knac.GT.0)
          knac2=2*knac
          pp=(wh-1)*nkn*n1+(knac-1)*n1
          pp2=(wh-1)*nkn*n1+(knac2-1)*n1
          pp2b=(wh-1)*nkn*n1+(knac2)*n1
          IF (conc(knac,wh,1).EQ.1) THEN
            DO j=1,n1
              storage(pp+j)=storage(pp2+j)*storage(pp2b+j)
            END DO
          ELSE  
            DO j=1,n1
              storage(pp+j)=1-(1-storage(pp2+j))*(1-storage(pp2b+j))
            END DO
          END IF
          nwkv=nwkv+1
          wkv(nwkv)=knac
          knac=INT(REAL(knac)/2.0)
        ENDDO

      END

      ! *****************************************************************
      ! *****************************************************************

      ! this subroutine evaluates the tree after a leaf was split
      ! last modification 10/16/02

      SUBROUTINE evaluate_branch(wh,knt,n1,n2,nkn,ntr,conc,term,negs,
     #                           datri,storage,nwkv,wkv)
      IMPLICIT NONE

        ! arguments in
          INTEGER knt,n1,n2,nkn,ntr,wh
          INTEGER conc(nkn,ntr,3)
          INTEGER negs(nkn,ntr,3)
          INTEGER term(nkn,ntr,3)
          INTEGER datri(n2,n1)
        ! local
          INTEGER j,knac,iterm,knac2,pp,pp2,pp2b
        ! arguments out
          INTEGER nwkv,wkv(nkn)
          INTEGER storage(2*ntr*nkn*n1)

      ! wkv stands for "Which Knots Visited" (used in restoring routine)
        nwkv=0
        DO j=1,nkn
          wkv(j)=0
        END DO

      ! store predictor in 4*knt and 4*knt+1 and 2*knt
        knac=2*knt+1
        iterm=term(knac,wh,1)
        pp=(wh-1)*nkn*n1+(knac-1)*n1
        IF (negs(knac,wh,1).EQ.0) THEN
          DO j=1,n1
            storage(pp+j)=datri(iterm,j)
          END DO
        ELSE 
          DO j=1,n1
            storage(pp+j)=1-datri(iterm,j)
          END DO
        END IF
        nwkv=nwkv+1
        wkv(nwkv)=knac
        knac=4*knt+1
        iterm=term(knac,wh,1)
        pp=(wh-1)*nkn*n1+(knac-1)*n1
        IF (negs(knac,wh,1).EQ.0) THEN
          DO j=1,n1
            storage(pp+j)=datri(iterm,j)
          END DO
        ELSE 
          DO j=1,n1
            storage(pp+j)=1-datri(iterm,j)
          END DO
        END IF
        nwkv=nwkv+1
        wkv(nwkv)=knac
        knac=4*knt
        iterm=term(knac,wh,1)
        pp=(wh-1)*nkn*n1+(knac-1)*n1
        IF (negs(knac,wh,1).EQ.0) THEN
          DO j=1,n1
            storage(pp+j)=datri(iterm,j)
          END DO
        ELSE 
          DO j=1,n1
            storage(pp+j)=1-datri(iterm,j)
          END DO
        END IF
        nwkv=nwkv+1
        wkv(nwkv)=knac
        knac=INT(REAL(knac)/2.0)

      ! calculate predictions in storage
        DO WHILE (knac.GT.0)
          knac2=2*knac
          pp=(wh-1)*nkn*n1+(knac-1)*n1
          pp2=(wh-1)*nkn*n1+(knac2-1)*n1
          pp2b=(wh-1)*nkn*n1+(knac2)*n1
          IF (conc(knac,wh,1).EQ.1) THEN
            DO j=1,n1
              storage(pp+j)=storage(pp2+j)*storage(pp2b+j)
            END DO
          ELSE  
            DO j=1,n1
              storage(pp+j)=1-(1-storage(pp2+j))*(1-storage(pp2b+j))
            END DO
          END IF
          nwkv=nwkv+1
          wkv(nwkv)=knac
          knac=INT(REAL(knac)/2.0)
        ENDDO

      END 

      ! *****************************************************************
      ! *****************************************************************



      ! this subroutine evaluates the tree after a leaf was split
      ! last modification 10/16/02

      SUBROUTINE evaluate_prune(wh,knt,n1,n2,nkn,ntr,conc,term,negs,
     #                          datri,storage,nwkv,wkv)
      IMPLICIT NONE

        ! arguments in
          INTEGER knt,n1,n2,nkn,ntr,wh
          INTEGER conc(nkn,ntr,3)
          INTEGER negs(nkn,ntr,3)
          INTEGER term(nkn,ntr,3)
          INTEGER datri(n2,n1)
        ! local
          INTEGER j,knac,iterm,knac2,pp,pp2,pp2b
        ! arguments out
          INTEGER nwkv,wkv(nkn)
          INTEGER storage(2*ntr*nkn*n1)

      ! wkv stands for "Which Knots Visited" (used in restoring routine)
        nwkv=0
        DO j=1,nkn
          wkv(j)=0
        END DO

      ! store predictor in 4*knt and 4*knt+1 and 2*knt
        knac=2*knt
        iterm=term(knac,wh,1)
        pp=(wh-1)*nkn*n1+(knac-1)*n1
        IF (negs(knac,wh,1).EQ.0) THEN
          DO j=1,n1
            storage(pp+j)=datri(iterm,j)
          END DO
        ELSE 
          DO j=1,n1
            storage(pp+j)=1-datri(iterm,j)
          END DO
        END IF
        nwkv=nwkv+1
        wkv(nwkv)=knac
        knac=2*knt+1
        iterm=term(knac,wh,1)
        pp=(wh-1)*nkn*n1+(knac-1)*n1
        IF (negs(knac,wh,1).EQ.0) THEN
          DO j=1,n1
            storage(pp+j)=datri(iterm,j)
          END DO
        ELSE 
          DO j=1,n1
            storage(pp+j)=1-datri(iterm,j)
          END DO
        END IF
        nwkv=nwkv+1
        wkv(nwkv)=knac
        knac=INT(REAL(knac)/2.0)

      ! calculate predictions in storage
        DO WHILE (knac.GT.0)
          knac2=2*knac
          pp=(wh-1)*nkn*n1+(knac-1)*n1
          pp2=(wh-1)*nkn*n1+(knac2-1)*n1
          pp2b=(wh-1)*nkn*n1+(knac2)*n1
          IF (conc(knac,wh,1).EQ.1) THEN
            DO j=1,n1
              storage(pp+j)=storage(pp2+j)*storage(pp2b+j)
            END DO
          ELSE  
            DO j=1,n1
              storage(pp+j)=1-(1-storage(pp2+j))*(1-storage(pp2b+j))
            END DO
          END IF
          nwkv=nwkv+1
          wkv(nwkv)=knac
          knac=INT(REAL(knac)/2.0)
        ENDDO

      END

      ! *****************************************************************
      ! *****************************************************************


      ! this subroutine initializes the variables needed in logic
      ! regression
      ! last modification 10/16/02

      SUBROUTINE initialize(n1,ntr,nkn,conc,term,negs,pick,storage,
     #                      score)
      IMPLICIT NONE

        ! arguments in
        INTEGER n1,ntr,nkn
        ! local
        INTEGER i,j,k,l,i2
        ! arguments out
        INTEGER conc(nkn,ntr,3)
        INTEGER negs(nkn,ntr,3)
        INTEGER pick(nkn,ntr,3)
        INTEGER term(nkn,ntr,3)
        INTEGER storage(2*ntr*nkn*n1)
        REAL score(3)

        ! tree specifications
        DO i=1,nkn
          DO j=1,ntr
            DO k=1,3
              conc(i,j,k)=0
              term(i,j,k)=0
              negs(i,j,k)=0
              pick(i,j,k)=0
            END DO
          END DO
        END DO

        ! stored values to compute the predictions of the trees
          DO j=1,ntr
            DO k=1,nkn
              i=n1*(k-1)+n1*nkn*(j-1)
              i2=n1*(k-1)+n1*nkn*(j-1)+ntr*n1*nkn
              DO l=1,n1
                storage(l+i)=0
                storage(l+i2)=0
             END DO
            END DO
          END DO

        ! others
        DO i=1,3
          score(i)=100000000.0
        ENDDO

      END

      ! *****************************************************************
      ! *****************************************************************

      ! this subroutine reads the selection probabilities/relations
      ! last modification 10/16/02

      SUBROUTINE selprob(nsp,cnc,slprbc)
      IMPLICIT NONE

        ! arguments in
          INTEGER nsp,cnc(3)
        ! local
          INTEGER j,k
          REAL denom,slprb(25)  
        ! arguments out
          REAL slprbc(25)

        slprb(1)=10
        slprb(2)=1
        slprb(3)=3
        slprb(4)=3
        slprb(5)=3
        slprb(6)=3

      ! correct for one operator only
        IF (cnc(2).EQ.1.OR.cnc(1).EQ.2) THEN
          slprb(2)=0.0
        END IF

      ! compute cumulative selection probabilities
        denom=0.0        
        DO j=1,nsp
          denom=denom+slprb(j)
        END DO
        DO j=1,nsp
          slprbc(j)=0
        END DO
        DO j=1,nsp
          DO k=1,j
            slprbc(j)=slprbc(j)+slprb(k)
          END DO
        END DO
        DO j=1,nsp
          slprb(j)=slprb(j)/denom
          slprbc(j)=slprbc(j)/denom
        END DO
        GOTO 2000
 2000   CONTINUE

      END 

      ! *****************************************************************
      ! *****************************************************************

      ! LOGIC REGRESSION MASTER VERSION (2)
      ! LAST MODIFICATION 10/16/02
      SUBROUTINE slogreg(n1,n2,nsep,intpars,rpars,seps,dcph,orders,resp,
     #                   weight,datri,iotrees,iocoef,ioscores,rd4)

      INTEGER LGCnknMAX,LGCntrMAX
      PARAMETER (LGCnknMAX  =   128)
      PARAMETER (LGCntrMAX  =     5)
      INTEGER n1,n2,nsep,cnc(3),mdl,msz,nkn,ehm,mszlo,mszup
      INTEGER ntrlo,ntrup,seed,kfold,nrep,choice,nfcnt,mtm
      INTEGER intpars(17),dcph(n1),orders(n1),datri(n2,n1)
      INTEGER iotrees(1),bout
      REAL seps(nsep,n1),resp(n1),weight(n1),rpars(14),tstr,tend,tint
      REAL penalty,iocoef(1),ioscores(1),hyperpars(10)
      INTEGER conc(LGCnknMAX,LGCntrMAX,3),xstop,rd4(1)
      INTEGER negs(LGCnknMAX,LGCntrMAX,3),i,j,k
      INTEGER pick(LGCnknMAX,LGCntrMAX,3)
      INTEGER term(LGCnknMAX,LGCntrMAX,3)

      DO i=1,LGCnknMAX
         DO j=1,LGCntrMAX
            DO k=1,3
               conc(i,j,k)=0
               negs(i,j,k)=0
               pick(i,j,k)=0
               term(i,j,k)=0
            END DO
         END DO
      END DO
      mdl=intpars(1)
      msz=intpars(2)
      mszlo=intpars(2)
      mszup=intpars(3)
      ntr=intpars(6)
      nkn=intpars(4)
      if(msz.LT.0) msz=ntr*nkn
      nkn=2*nkn-1
      ntrlo=intpars(5)
      ntrup=intpars(6)
      IF (intpars(7).EQ.2) THEN
          cnc(1)=1
          cnc(2)=1
        ELSE IF (intpars(7).EQ.3) THEN
          cnc(1)=2
          cnc(2)=2
        ELSE
          cnc(1)=1
          cnc(2)=2
        END IF

      tstr=rpars(1)
      tend=rpars(2)
      tint=rpars(3)
      ehm=intpars(8)
      seed=intpars(9)
      kfold=intpars(10)
      nrep=intpars(10)
      choice=intpars(11)
      nfcnt=intpars(12)
      penalty=rpars(4)
      mtm=intpars(13)
      xstop=0
      bout=0
      DO i=1,10
         hyperpars(i)=rpars(i+4)
      END DO
      IF(choice.EQ.7)THEN
         tstr=1.
         tend=1.
         tint=intpars(14)+intpars(15)
         nfcnt=intpars(14)
         bout=intpars(16)
      END IF
      CALL stopper(LGCnknMAX,nkn,"LGCnknMAX","slogreg()",9,9,xstop,0)
      CALL stopper(LGCntrMAX,ntr,"LGCntrMAX","slogreg()",9,9,xstop,1)
      CALL logreg(mdl,msz,n1,n2,nkn,ntr,cnc,nsep,tstr,tend,tint,ehm,
     #                  mszlo,mszup,ntrlo,ntrup,seed,kfold,nrep,
     #                  choice,nfcnt,penalty,mtm,seps,dcph,orders,
     #                  resp,weight,datri,iotrees,ltree,iocoef,
     #                  ioscores,conc,negs,pick,term,hyperpars,rd4,bout)
      intpars(1)=ltree
      END
      

      SUBROUTINE logreg(mdl,msz,n1,n2,nkn,ntr,cnc,nsep,tstr,tend,tint,
     #                  ehm,mszlo,mszup,ntrlo,ntrup,seed,kfold,
     #                  nrep,choice,nfcnt,penalty,mtm,seps,dcph,
     #                  ordrs,resp,weight,datri,iotree,iosclast,iocoef,
     #                  ioscores,conc,negs,pick,term,hyperpars,rd4,
     #                  bout)
      IMPLICIT NONE
        
        ! parameters to be determined by the user
        ! INTEGER n1    number of subjects
        ! INTEGER n2    number of predictors
        ! INTEGER mdl   model type under consideration
        ! INTEGER ntr   number of trees
        ! INTEGER nkn   number of knots
        ! INTEGER msz   maximal size of model
        ! INTEGER nsep  how many separate predictors

        ! how many selection probabilities/relations:
  
        ! parameters
          INTEGER LGCn1MAX,LGCnknMAX,LGCntrMAX,ltree,iotree(1)
          INTEGER LGCbetaMAX
          PARAMETER (LGCn1MAX   = 10000)
          PARAMETER (LGCnknMAX  =   128)
          PARAMETER (LGCntrMAX  =     5)
          PARAMETER (LGCbetaMAX =    55)

          INTEGER i,j,ssize,ehm,seed,nfcnt,error,nsp,k,rd4(1),bout
          INTEGER choice,kfold,mdl,msz,n1,n2,nkn,nsep,ntr
          INTEGER ntrnew,mszlo,mszup,nrep,ntrlo,ntrup
          INTEGER cnc(3),mtm
          INTEGER dcph(n1),ordrs(n1)
          INTEGER conc(nkn,ntrup,3)
          INTEGER negs(nkn,ntrup,3)
          INTEGER pick(nkn,ntrup,3)
          INTEGER term(nkn,ntrup,3)
          INTEGER storage(2*LGCntrMAX*LGCnknMAX*LGCn1MAX)
          INTEGER*4 sseed(3)
          REAL tstr,tend,tint,score(3),slprbc(25),penalty,hyperpars(10)
          REAL weight(n1),betas(3,0:LGCbetaMAX),resp(n1),iocoef(1)
          INTEGER datri(n2,n1),mcmc
          REAL seps(nsep,n1),rnumsrr(LGCn1MAX),ioscores(n1),myrand
          INTEGER mszmax,xstop,iolast,iosclast
          INTEGER prtr(LGCn1MAX,LGCntrMAX)
          REAL cbetas(0:LGCbetaMAX),xtxsep(LGCbetaMAX+1,LGCbetaMAX+1)
          INTEGER npckmv(6,LGCntrMAX)
          INTEGER pickmv(6,LGCnknMAX,LGCntrMAX)
          CHARACTER *100 astring
          
        ! betas: the parameter estimates in the model - current/stored/best
        ! choice: which type of logic regression
        ! cnc: specifies which operators will be considered
        ! conc: operators in the logic trees
        ! datr: data stored as reals
        ! dcph: indicator for censoring (mdl=4 only)
        ! ehm: every how many iterations you get a score update
        ! kfold: k in k-fold cross validation
        ! negs: indicator (complement 0/1) in the logic trees
        ! nrep: number of permutations in the randomization test
        ! ordrs: order of responses (mdl=4 only)
        ! pick: indicator (taken 0/1) in the logic trees
        ! score: scores of current/stored/best model,resp(LGCn1MAX)
        ! seps: separate variables to condition on
        ! ssize: size of the current model
        ! slprbc: cumulative selection probabilities
        ! storage: array that stores the values in the knots of the 
        !   logic trees
        ! term: predictors in the logic trees
        ! weight: case weights
        ! nfcnt: number of iterations after which to check whether there
        !        are 10 changes
        ! penalty: penalty parameter for model size
        iolast=0
        iosclast=0
        ltree=0
        xstop=0
        mcmc=0
        CALL stopper(LGCnknMAX,nkn,"LGCnknMAX","logreg()",9,9,xstop,0)
        CALL stopper(LGCntrMAX,ntr,"LGCntrMAX","logreg()",9,9,xstop,0)
        CALL stopper(LGCn1MAX,n1,"LGCn1MAX","logreg()",8,9,xstop,0)
        CALL stopper(LGCbetaMAX,nsep+ntr,"LGCbetaMAX","logreg()",10,9,
     #          xstop,1)
        DO i=1,3
          DO j=0,LGCbetaMAX
             betas(i,j)=0.
          END DO
        END DO
      ! read in data, weights, selection probabilities, etc
        DO i=1,LGCn1MAX
           DO j=1,LGCntrMAX
              prtr(i,j)=0
           END DO
        END DO
        error=0
        nsp=6
        CALL selprob(nsp,cnc,slprbc)

      ! set random seed
        
        IF(seed.NE.0)THEN
          sseed(1)=seed
          IF(sseed(1).LT.0)THEN
             sseed(1)= -seed
          END IF
          i=myrand(sseed(1))
        ELSE
          i=myrand(1773)
        END IF

      ! find the single best model -----------------------------------
        IF (choice.EQ.1) THEN
          ltree=ntr
          CALL annealing(n1,n2,mdl,nkn,nsp,ntr,conc,negs,pick,term,
     #                   storage,slprbc,datri,weight,tstr,tend,tint,ehm,
     #                   msz,nsep,seps,cnc,score,betas,ssize,dcph,ordrs,
     #                   nfcnt,penalty,resp,mtm,mcmc,hyperpars,iotree,
     #                   iocoef,ioscores,rd4,bout)
          DO j=0,(nsep+ntr)
             iolast=iolast+1
             iocoef(iolast)=betas(3,j)
          END DO
          DO j=1,ntr
            iotree((j-1)*(4*nkn+3)+1)=msz
            iotree((j-1)*(4*nkn+3)+2)=ntr
            iotree((j-1)*(4*nkn+3)+3)=j
            DO k=1,nkn
              iotree((j-1)*(4*nkn+3)+(k-1)*4+4)=conc(k,j,3)
              iotree((j-1)*(4*nkn+3)+(k-1)*4+5)=term(k,j,3)
              iotree((j-1)*(4*nkn+3)+(k-1)*4+6)=negs(k,j,3)
              iotree((j-1)*(4*nkn+3)+(k-1)*4+7)=pick(k,j,3)
            END DO
          END DO
          ioscores(1)=score(3)

      ! MCMC -----------------------------------
        ELSE IF (choice.EQ.7) THEN
          mcmc=1
          ltree=ntr
          CALL annealing(n1,n2,mdl,nkn,nsp,ntr,conc,negs,pick,term,
     #                   storage,slprbc,datri,weight,tstr,tend,tint,ehm,
     #                   msz,nsep,seps,cnc,score,betas,ssize,dcph,ordrs,
     #                   nfcnt,penalty,resp,mtm,mcmc,hyperpars,iotree,
     #                   iocoef,ioscores,rd4,bout)

      ! find multiple models -----------------------------------------
        ELSE IF (choice.EQ.6) THEN
          iolast=0
          iosclast=0
          ntr=ntrup
          DO msz=0,mszup
            IF(msz.eq.0)THEN
              CALL annealing_init(n1,n2,mdl,nkn,ntr,conc,negs,pick,term,
     #        storage,datri,resp,weight,npckmv,pickmv,
     #        score,betas,ssize,nsep,seps,ntrnew,dcph,ordrs,penalty,
     #        prtr,cbetas,xtxsep,mtm,mcmc,hyperpars)
            ELSE
              CALL greedyonestep(nkn,ntr,conc,negs,pick,term,n1,n2,
     #        storage,prtr,datri,nsep,seps,cbetas,xtxsep,mdl,
     #        ntrnew,ordrs,dcph,resp,weight,score(1),mtm,cnc,
     #        penalty)
              IF(ntr.LT.0)THEN
                IF(ehm.GE.0)THEN
                  astring(1:32)="No further improvement possible"
                  CALL stringprint(astring,32)
                END IF
                GOTO 6677
              END IF
            END IF
            IF (ehm.GE.0) THEN
              astring(1:6)="Model "
              CALL makeistring(7,9,astring,msz,3)
              astring(10:24)=" has a score of"
              CALL makerstring(25,37,astring,score(1),8,4)
              CALL stringprint(astring,37)
            END IF
            DO j=0,(nsep+ntrnew)
              iolast=iolast+1
              iocoef(iolast)=cbetas(j)
            END DO
            DO j=(ntrnew+1),ntr
              iolast=iolast+1
              iocoef(iolast)=0.
            END DO
            DO j=1,ntrnew
              ltree=ltree+1
              iotree((ltree-1)*(4*nkn+3)+1)=msz
              iotree((ltree-1)*(4*nkn+3)+2)=ntr
              iotree((ltree-1)*(4*nkn+3)+3)=j
              DO k=1,nkn
                CALL reorder(iotree,conc,term,negs,pick,ltree,k,j,
     #                         nkn,ntrlo)
              END DO
            END DO
            DO j=(ntrnew+1),ntr
              ltree=ltree+1
              iotree((ltree-1)*(4*nkn+3)+1)=msz
              iotree((ltree-1)*(4*nkn+3)+2)=ntr
              iotree((ltree-1)*(4*nkn+3)+3)=j
              DO k=1,nkn
                iotree((ltree-1)*(4*nkn+3)+(k-1)*4+4)=0
                iotree((ltree-1)*(4*nkn+3)+(k-1)*4+5)=0
                iotree((ltree-1)*(4*nkn+3)+(k-1)*4+6)=0
                iotree((ltree-1)*(4*nkn+3)+(k-1)*4+7)=0
              END DO
            END DO
            iosclast=iosclast+1
            ioscores(iosclast)=score(1)
          END DO
6677      CONTINUE

      ! find multiple models -----------------------------------------
        ELSE IF (choice.EQ.2) THEN
          ltree=0
          DO ntr=ntrlo,ntrup
            IF (ehm.GE.0) THEN
              CALL stringprint('  ',2)
              astring(1:39)="The number of trees in these models is "
              CALL makeistring(40,42,astring,ntr,3)
              CALL stringprint(astring,42)
            ENDIF
            DO msz=mszlo,mszup
              IF (msz.GE.0.AND.msz.LT.ntr.AND.ntr.GT.ntrlo) THEN
                IF(ehm.GE.0)THEN
                  CALL stringprint('  ',2)
                 astring(1:40)='The size for this model is smaller than'
                  astring(41:76)=' the number of trees you requested.'
                  CALL stringprint(astring,76)
                  astring(1:1)="("
                  CALL makeistring(2,5,astring,msz,2)
                  astring(4:10)=" versus"
                  CALL makeistring(11,13,astring,ntr,3)
                  astring(14:14)=")"
                  CALL stringprint(astring,14)
                  CALL stringprint(
     #               'To save CPU time, we will skip this run.',40)
                  CALL stringprint("On to the next model...",23)
                END IF
              ELSE
                IF (ehm.GE.0.AND.msz.GT.0) THEN
                  CALL stringprint('  ',2)
                  astring(1:18)="The model size is "
                  CALL makeistring(19,21,astring,msz,3)
                  CALL stringprint(astring,21)
                ENDIF
                ntrnew=MIN(ntr,msz)
                CALL annealing(n1,n2,mdl,nkn,nsp,ntrnew,
     #                       conc,negs,pick,term,
     #                       storage,slprbc,datri,weight,tstr,tend,tint,
     #                       ehm,msz,nsep,seps,cnc,score,betas,ssize,
     #                       dcph,ordrs,nfcnt,penalty,resp,mtm,
     #                       mcmc,hyperpars,iotree,
     #                       iocoef,ioscores,rd4,bout)

                IF (ehm.EQ.0) THEN
                  astring(1:20)="The best model with "
                  CALL makeistring(21,23,astring,ntr,3)
                  astring(24:37)=" trees of size "
                  CALL makeistring(38,40,astring,msz,3)
                  astring(41:55)=" has a score of"
                  CALL makerstring(56,68,astring,score(3),8,4)
                  CALL stringprint(astring,68)
                END IF
                DO j=0,(nsep+ntrnew)
                   iolast=iolast+1
                   iocoef(iolast)=betas(3,j)
                END DO
                DO j=(ntrnew+1),ntr
                   iolast=iolast+1
                   iocoef(iolast)=0.
                END DO
                DO j=1,ntrnew
                  ltree=ltree+1
                  iotree((ltree-1)*(4*nkn+3)+1)=msz
                  iotree((ltree-1)*(4*nkn+3)+2)=ntr
                  iotree((ltree-1)*(4*nkn+3)+3)=j
                  DO k=1,nkn
                    CALL reorder(iotree,conc,term,negs,pick,ltree,k,j,
     #                           nkn,ntrnew)
                  END DO
                END DO
                DO j=(ntrnew+1),ntr
                  ltree=ltree+1
                  iotree((ltree-1)*(4*nkn+3)+1)=msz
                  iotree((ltree-1)*(4*nkn+3)+2)=ntr
                  iotree((ltree-1)*(4*nkn+3)+3)=j
                  DO k=1,nkn
                    iotree((ltree-1)*(4*nkn+3)+(k-1)*4+4)=0
                    iotree((ltree-1)*(4*nkn+3)+(k-1)*4+5)=0
                    iotree((ltree-1)*(4*nkn+3)+(k-1)*4+6)=0
                    iotree((ltree-1)*(4*nkn+3)+(k-1)*4+7)=0
                  END DO
                END DO
                iosclast=iosclast+1
                ioscores(iosclast)=score(3)
              END IF
            END DO
          END DO

      ! use cross validation for model selection ---------------------
        ELSE IF (choice.EQ.3) THEN
          DO i=1,n1
             rnumsrr(i)=myrand(0)
          END DO
          ltree=0
          DO ntr=ntrlo,ntrup
            IF (ehm.GE.0) THEN
              CALL stringprint('  ',2)
              astring(1:39)="The number of trees in these models is "
              CALL makeistring(40,42,astring,ntr,3)
              CALL stringprint(astring,42)
            ENDIF
            DO msz=mszlo,mszup
              IF (msz.GE.0.AND.msz.LT.ntr.AND.ntr.GT.ntrlo) THEN
                IF (ehm.GE.0) THEN
                  CALL stringprint('  ',2)
                 astring(1:40)='The size for this model is smaller than'
                  astring(41:76)=' the number of trees you requested.'
                  CALL stringprint(astring,76)
                  astring(1:1)="("
                  CALL makeistring(2,5,astring,msz,2)
                  astring(4:10)=" versus"
                  CALL makeistring(11,13,astring,ntr,3)
                  astring(14:14)=")"
                  CALL stringprint(astring,14)
                  CALL stringprint(
     #               'To save CPU time, we will skip this run.',40)
                  CALL stringprint("On to the next model...",23)
                END IF
              ELSE
              IF (ehm.GE.0) THEN
                CALL stringprint('  ',2)
                astring(1:18)="The model size is "
                CALL makeistring(19,21,astring,msz,3)
                CALL stringprint(astring,21)
              ENDIF
              CALL crossvalx(kfold,n1,n2,mdl,nkn,nsp,ntr,conc,negs,
     #                     pick,term,storage,slprbc,datri,weight,
     #                     tstr,tend,tint,ehm,msz,nsep,seps,cnc,score,
     #                     betas,ssize,dcph,ordrs,nfcnt,seed,resp,
     #                     rnumsrr,mtm,ltree,ioscores)
              ENDIF
            END DO
          END DO

      ! run the null model test for signal in the data ---------------
        ELSE IF (choice.EQ.4) THEN
          CALL annealing(n1,n2,mdl,nkn,nsp,ntr,conc,negs,pick,term,
     #                   storage,slprbc,datri,weight,tstr,tend,tint,ehm,
     #                   0,nsep,seps,cnc,score,betas,ssize,dcph,ordrs,
     #                   nfcnt,penalty,resp,mtm,mcmc,hyperpars,iotree,
     #                   iocoef,ioscores,rd4,bout)
          ioscores(1)=score(1)
          CALL annealing(n1,n2,mdl,nkn,nsp,ntr,conc,negs,pick,term,
     #                   storage,slprbc,datri,weight,tstr,tend,tint,ehm,
     #                   msz,nsep,seps,cnc,score,betas,ssize,dcph,ordrs,
     #                   nfcnt,penalty,resp,mtm,mcmc,hyperpars,iotree,
     #                   iocoef,ioscores,rd4,bout)
          ioscores(2)=score(3)
          IF (ehm.EQ.0) THEN
            CALL stringprint('  ',2)
            astring(1:25)="The best model has score "
            CALL makerstring(26,40,astring,score(3),10,4)
            CALL stringprint(astring,40)
          END IF
          DO j=1,nrep
            CALL nullmodelx(n1,n2,mdl,nkn,nsp,ntr,conc,negs,pick,term,
     #                     storage,slprbc,datri,weight,
     #                     tstr,tend,tint,ehm,msz,nsep,seps,cnc,
     #                     score,betas,ssize,dcph,ordrs,nfcnt,resp,mtm)
            ioscores(j+2)=score(3)
            IF (ehm.GE.0) THEN  
              astring(1:18)="Permutation number"
              CALL makeistring(19,23,astring,j,5)
              astring(24:31)="  out of"
              CALL makeistring(32,36,astring,nrep,5)
              astring(37:48)="  has score"
              CALL makerstring(49,62,astring,score(3),9,4)
              CALL stringprint(astring,62)
            END IF
          END DO

      ! run the null model test for signal in the data ---------------
        ELSE IF (choice.EQ.5) THEN
          msz=nkn*ntrup
          CALL annealing(n1,n2,mdl,nkn,nsp,ntrup,conc,negs,pick,term,
     #                   storage,slprbc,datri,weight,tstr,tend,tint,ehm,
     #                 msz,nsep,seps,cnc,score,betas,ssize,dcph,ordrs,
     #                   nfcnt,penalty,resp,mtm,mcmc,hyperpars,iotree,
     #                   iocoef,ioscores,rd4,bout)
          ioscores(2)=score(3)
          CALL annealing(n1,n2,mdl,nkn,nsp,ntrup,conc,negs,pick,term,
     #                   storage,slprbc,datri,weight,tstr,tend,tint,ehm,
     #                   0,nsep,seps,cnc,score,betas,ssize,dcph,ordrs,
     #                   nfcnt,penalty,resp,mtm,mcmc,hyperpars,iotree,
     #                   iocoef,ioscores,rd4,bout)
          ioscores(1)=score(1)
          ltree=0
          DO ntr=ntrlo,ntrup
            IF (ehm.GE.0) THEN
              CALL stringprint('  ',2)
              astring(1:39)="The number of trees in these models is "
              CALL makeistring(40,42,astring,ntr,3)
              CALL stringprint(astring,42)
            ENDIF
            DO msz=mszlo,mszup
              IF (msz.GE.0.AND.msz.LT.ntr.AND.ntr.GT.ntrlo) THEN
                IF(ehm.GE.0)THEN
                  CALL stringprint('  ',2)
                 astring(1:40)='The size for this model is smaller than'
                  astring(41:76)=' the number of trees you requested.'
                  CALL stringprint(astring,76)
                  astring(1:1)="("
                  CALL makeistring(2,5,astring,msz,2)
                  astring(4:10)=" versus"
                  CALL makeistring(11,13,astring,ntr,3)
                  astring(14:14)=")"
                  CALL stringprint(astring,14)
                  CALL stringprint(
     #               'To save CPU time, we will skip this run.',40)
                  CALL stringprint("On to the next model...",23)
                END IF
              ELSE
                ltree=ltree+1
                ntrnew=MIN(ntr,msz)
                IF (ehm.GE.0) THEN
                  CALL stringprint('  ',2)
                  astring(1:18)="The model size is "
                  CALL makeistring(19,21,astring,msz,3)
                  CALL stringprint(astring,21)
                ENDIF
                ioscores((ltree-1)*(nrep+2)+3)=ntrnew
                ioscores((ltree-1)*(nrep+2)+4)=msz
                DO j=1,nrep
                  IF(ntrnew.eq.0)THEN
                    mszmax=ntrup*nkn
                    CALL nullmodelx(n1,n2,mdl,nkn,nsp,ntrup,conc,negs,
     #                            pick,term,storage,slprbc,datri,weight,
     #                            tstr,tend,tint,ehm,mszmax,nsep,seps,
     #                            cnc,score,betas,ssize,dcph,ordrs,
     #                            nfcnt,resp,mtm)
                  ELSE
                    CALL randomizationx(n1,n2,mdl,nkn,nsp,ntrnew,
     #                                conc,negs,pick,term,storage,
     #                                slprbc,datri,weight,tstr,tend,
     #                                tint,ehm,msz,nsep,seps,cnc,score,
     #                                betas,ssize,dcph,ordrs,iotree,
     #                                nfcnt,ntrup,error,resp,mtm)
                  END IF
                  IF(error.eq.1)THEN
                    astring(1:23)="No data for model size "
                    CALL makeistring(24,27,astring,msz,4)
                    astring(28:33)=" with "
                     CALL makeistring(34,36,astring,ntrnew,3)
                     astring(37:44)=" tree(s)"
                     CALL stringprint(astring,44)
                     ioscores((ltree-1)*(nrep+2)+3)= -1
                     GO TO 1299
                  END IF
                  IF (ehm.GE.0) THEN  
                     astring(1:18)="Permutation number"
                     CALL makeistring(19,23,astring,j,5)
                     astring(24:31)="  out of"
                     CALL makeistring(32,36,astring,nrep,5)
                     astring(37:47)="  has score"
                     CALL makerstring(48,61,astring,score(3),9,4)
                     astring(61:73)="   model size"
                     CALL makeistring(74,76,astring,msz,4)
                     astring(77:81)=" with"
                     CALL makeistring(82,84,astring,ntrnew,3)
                     astring(85:92)=" tree(s)"
                     CALL stringprint(astring,92)
                  END IF
                  ioscores((ltree-1)*(nrep+2)+j+4)=score(3)
                END DO
              END IF
1299          CONTINUE
            END DO
          END DO
        END IF

 105    FORMAT(A18,I5,2X,A6,I5,A11,F14.5,10X,A10,I3,A5,I2,A8)
 108    FORMAT(A22,I3,A5,I2,A8)

      END 

      ! *****************************************************************
      ! *****************************************************************

      ! this subroutine selects a move for a given tree
      ! last modification 10/16/02

      SUBROUTINE moving(n2,nkn,nsp,ntr,conc,negs,pick,term,slprbc,cnc,
     #                mcmc,npckmv,pickmv,msz,ssize,nop,wh,knt,mtp,mctry)
      IMPLICIT NONE

        ! arguments in
          INTEGER n2,msz,nkn,nsp,ntr,nop,ssize,cnc(3),npckmv(6,ntr)
          INTEGER pickmv(6,nkn,ntr),mcmc,mctry
          INTEGER conc(nkn,ntr,3)
          INTEGER negs(nkn,ntr,3)
          INTEGER pick(nkn,ntr,3)
          INTEGER term(nkn,ntr,3)
          REAL slprbc(25)
        ! local
          REAL rnum,myrand,mtpx,whx
          INTEGER i,j
        ! arguments out
          INTEGER knt,mtp,wh

      ! randomly select tree 
 1000   CONTINUE
        CALL copytree(ntr,nkn,conc,negs,pick,term,-1,2,1)
        rnum=myrand(0)
        IF(mctry.GT.0)THEN
          IF(mtp.NE.0)THEN
            DO j=0,(nop)
              DO i=1,2
                IF(i+j.GT.1)THEN
                  whx=wh+j
                  IF(whx.GT.nop)whx=whx-nop
                  mtpx=mtp
                  IF(mtp.LT.3.and.i.EQ.2)mtpx=3-mtp
                  IF(mtp.GE.3.and.i.EQ.2)mtpx=9-mtp
                  IF(npckmv(mtpx,whx).GT.0)THEN
                    mtp=mtpx
                    wh=whx
                    CALL decision(n2,nkn,nsp,ntr,wh,conc,negs,pick,
     #              term,slprbc,cnc,npckmv,pickmv,msz,ssize,nop,
     #              knt,mtp,mcmc,mctry)
                    GOTO 169
                  END IF
                END IF
              END DO
            END DO
          ELSE
            CALL firstknot(n2,nkn,ntr,wh,conc,negs,pick,term,-1,-1)
            knt=0
            mtp=0
            GOTO 169
          END IF
        END IF
        mtp=-1
        mctry=0
      ! wh=INT(REAL(MIN(ntr,nop+1))*rnum)+1
        wh=INT(REAL(ntr)*rnum)+1
        IF(wh.GT.nop+1)wh=nop+1
        IF (pick(1,wh,1).EQ.0) THEN
          IF (ssize.EQ.msz) GOTO 1000 
          IF(mcmc.GT.0)mcmc=2
          CALL firstknot(n2,nkn,ntr,wh,conc,negs,pick,term,-1,-1)
          knt=0
          mtp=0
        ELSE
          CALL decision(n2,nkn,nsp,ntr,wh,conc,negs,pick,term,slprbc,
     #                  cnc,npckmv,pickmv,msz,ssize,nop,
     #                  knt,mtp,mcmc,mctry)
        END IF
169     CONTINUE
      END 

      ! *****************************************************************
      ! *****************************************************************

      ! this subroutine starts a new tree
      ! last modification 10/16/02

      SUBROUTINE firstknot(n2,nkn,ntr,wh,conc,negs,pick,term,r1,r2)
      IMPLICIT NONE

        ! arguments in
          INTEGER n2,nkn,ntr,wh,r1,r2
          INTEGER conc(nkn,ntr,3)
          INTEGER negs(nkn,ntr,3)
          INTEGER pick(nkn,ntr,3)
          INTEGER term(nkn,ntr,3)
        ! local
          INTEGER letter,neg
          REAL rnum,myrand
        ! arguments out

      ! randomly select predictor
        IF(r1.LT.0)THEN
           rnum=myrand(0)
           letter=INT(REAL(n2)*rnum)+1
           rnum=myrand(0)
           neg=INT(2.0*rnum)
        ELSE
           letter=r1
           neg=r2
        END IF

        conc(1,wh,1)=3
        term(1,wh,1)=letter
        negs(1,wh,1)=neg
        pick(1,wh,1)=1

      END 

      ! *****************************************************************
      ! *****************************************************************



      ! this subroutine makes a decision which move to carry out
      ! last modification 10/16/02

      SUBROUTINE decision(n2,nkn,nsp,ntr,wh,conc,negs,pick,term,slprbc,
     #                    cnc,npckmv,pickmv,msz,ssize,nop,
     #                    knt,mtp,mcmc,mctry)
      IMPLICIT NONE

        ! arguments in
          INTEGER msz,n2,nkn,nop,nsp,ntr,ssize,wh,cnc(3)
          INTEGER npckmv(6,ntr),pickmv(6,nkn,ntr)
          INTEGER conc(nkn,ntr,3)
          INTEGER negs(nkn,ntr,3)
          INTEGER pick(nkn,ntr,3)
          INTEGER term(nkn,ntr,3),mcmc
          REAL slprbc(25)
        ! local
          INTEGER sibling,sng,dbl,mctry
          REAL rnum,myrand
        ! arguments out
          INTEGER knt,mtp,rnd1,rnd2,rnd3

      ! select move from selection scheme
        IF(mcmc.GT.0)mcmc=1
        IF(mcmc.GT.0)mcmc=mcmc+1
        IF(mcmc.EQ.0)mctry=0
 1000   CONTINUE
        rnum=myrand(0)
        IF(mctry.GT.0)rnum=3.
 2000   CONTINUE
        IF (rnum.LE.slprbc(1).OR.(mctry.GT.0.AND.mtp.EQ.1)) THEN
          rnum=myrand(0)
          knt=pickmv(1,INT(REAL(npckmv(1,wh))*rnum)+1,wh)
          rnd1= -1
          CALL altlf(knt,n2,nkn,ntr,wh,conc,negs,term,rnd1,rnd2)
          mtp=1
        ELSE IF (rnum.LT.slprbc(2).OR.(mctry.GT.0.AND.mtp.EQ.2)) THEN
          IF (npckmv(2,wh).EQ.0) GOTO 1000
          rnum=myrand(0)
          knt=pickmv(2,INT(REAL(npckmv(2,wh))*rnum)+1,wh)
          CALL altop(knt,nkn,ntr,wh,conc)
          mtp=2
        ELSE IF (rnum.LT.slprbc(3).OR.(mctry.GT.0.AND.mtp.EQ.3)) THEN
          rnum=myrand(0)
          IF (npckmv(3,wh).EQ.0) GOTO 1000
          knt=pickmv(3,INT(REAL(npckmv(3,wh))*rnum)+1,wh)
          CALL xdelete(knt,nkn,ntr,wh,conc,negs,pick,term)
          mtp=3
        ELSE IF (rnum.LT.slprbc(4).OR.(mctry.GT.0.AND.mtp.EQ.4)) THEN
          IF (ssize.EQ.msz) GOTO 1000
          rnum=myrand(0)
          IF (npckmv(4,wh).EQ.0) GOTO 1000
          knt=pickmv(4,INT(REAL(npckmv(4,wh))*rnum)+1,wh)
          rnd1= -1
          CALL xsplit(knt,n2,nkn,ntr,wh,cnc,conc,negs,pick,term,
     #                rnd1,rnd2,rnd3)
          mtp=4
        ELSE IF (rnum.LT.slprbc(5).OR.(mctry.GT.0.AND.mtp.EQ.5)) THEN
          IF (ssize.EQ.msz) GOTO 1000
          rnum=myrand(0)
          IF (npckmv(5,wh).EQ.0) GOTO 1000
          knt=pickmv(5,INT(REAL(npckmv(5,wh))*rnum)+1,wh)
          rnd1= -1
          CALL branch(knt,n2,nkn,ntr,wh,cnc,conc,negs,pick,term,
     #                rnd1,rnd2,rnd3)
          mtp=5
        ELSE IF (rnum.LE.slprbc(6).OR.(mctry.GT.0.AND.mtp.EQ.6)) THEN
          IF (npckmv(6,wh).EQ.0) GOTO 1000
          rnum=myrand(0)
          knt=pickmv(6,INT(REAL(npckmv(6,wh))*rnum)+1,wh)
          IF(knt.GT.0)THEN
            sng=2*knt
            dbl=2*knt+1
          ELSE
            knt= -knt
            sng=2*knt+1
            dbl=2*knt
          END IF
          CALL prune(knt,dbl,sng,nkn,ntr,wh,conc,negs,pick,term)
          mtp=6
        ELSE
          CALL stringprint("error: this move is not defined!",32)
          STOP
        END IF
      END 

      ! *****************************************************************
      ! *****************************************************************

      ! this subroutine alternates a leaf
      ! last modification 10/17/02

      SUBROUTINE altlf(knt,n2,nkn,ntr,wh,conc,negs,term,rnd1,rnd2)
      IMPLICIT NONE

        ! arguments in
          INTEGER knt,n2,nkn,ntr,wh,rnd1,rnd2
          INTEGER conc(nkn,ntr,3)
          INTEGER negs(nkn,ntr,3)
          INTEGER term(nkn,ntr,3)
        ! local
          INTEGER letter,neg
          REAL rnum,myrand
        ! arguments out

      ! alternate leaf
 300    CONTINUE
        IF(rnd1.GT.0)THEN
          letter=rnd1
        ELSE
          rnum=myrand(0)
          letter=INT(REAL(n2)*rnum)+1
        END IF
      ! avoid creating redundancies such as (X and X^c)
        IF (knt.GT.1) THEN
          IF (MOD(knt,2).EQ.0) THEN
            IF (letter.EQ.term(knt+1,wh,1)) THEN
              IF(rnd1.GT.0)GOTO 388
              GOTO 300
            END IF
          ELSE
            IF (letter.EQ.term(knt-1,wh,1)) THEN
              IF(rnd1.GT.0)THEN
                 rnd1=-1
                 GOTO 388
              END IF
              GOTO 300
            END IF
          ENDIF
        END IF
        IF(rnd1.GT.0)THEN
          neg=rnd2
        ELSE
          rnum=myrand(0)
          neg=INT(2.0*rnum)
        END IF
        term(knt,wh,1)=letter         
        negs(knt,wh,1)=neg
388     CONTINUE
      END 

      ! *****************************************************************
      ! *****************************************************************

      ! this subroutine alternates an operator
      ! last modification 10/17/02

      SUBROUTINE altop(knt,nkn,ntr,wh,conc)
      IMPLICIT NONE

        ! arguments in
          INTEGER knt,nkn,ntr,wh
          INTEGER conc(nkn,ntr,3)
        ! arguments out

      ! exchange operator
        conc(knt,wh,1)=3-conc(knt,wh,1)       

      END 

      ! *****************************************************************
      ! *****************************************************************



      ! this subroutine deletes a leaf
      ! last modification 10/17/02

      SUBROUTINE xdelete(knt,nkn,ntr,wh,conc,negs,pick,term)
      IMPLICIT NONE

        ! arguments in
          INTEGER knt,nkn,ntr,wh
          INTEGER conc(nkn,ntr,3)
          INTEGER negs(nkn,ntr,3)
          INTEGER pick(nkn,ntr,3)
          INTEGER term(nkn,ntr,3)
        ! local
          INTEGER head,sibling,token
        ! arguments out

      ! delete the leaf 
        IF (knt.EQ.1) THEN
          conc(knt,wh,1)=0
          term(knt,wh,1)=0
          negs(knt,wh,1)=0
          pick(knt,wh,1)=0
        ELSE
          head=INT(knt/2)
          IF (MOD(knt,2).EQ.0) THEN
            sibling=knt+1
          ELSE
            sibling=knt-1
          END IF
          conc(head,wh,1)=conc(sibling,wh,1)
          term(head,wh,1)=term(sibling,wh,1)
          negs(head,wh,1)=negs(sibling,wh,1)
          token=knt
          conc(token,wh,1)=0
          term(token,wh,1)=0
          negs(token,wh,1)=0
          pick(token,wh,1)=0
          token=sibling
          conc(token,wh,1)=0
          term(token,wh,1)=0
          negs(token,wh,1)=0
          pick(token,wh,1)=0 
        END IF

      END 

      ! *****************************************************************
      ! *****************************************************************

      ! this subroutine splits a leaf
      ! last modification 10/17/02

      SUBROUTINE xsplit(knt,n2,nkn,ntr,wh,cnc,conc,negs,pick,term,r1,r2,
     #                  r3)
      IMPLICIT NONE


        ! arguments in
          INTEGER knt,n2,nkn,ntr,wh,r1,r2,r3
          INTEGER cnc(3)
          INTEGER conc(nkn,ntr,3)
          INTEGER negs(nkn,ntr,3)
          INTEGER pick(nkn,ntr,3)
          INTEGER term(nkn,ntr,3)
        ! local
          INTEGER letter,neg
          REAL rnum,myrand
        ! arguments out

      ! move leaf and create operator
        conc(2*knt,wh,1)=3
        term(2*knt,wh,1)=term(knt,wh,1)
        negs(2*knt,wh,1)=negs(knt,wh,1)
        pick(2*knt,wh,1)=1
        term(knt,wh,1)=0
        negs(knt,wh,1)=0        
      ! create the new leaf 
        IF(r1.GE.0)THEN
          conc(knt,wh,1)=cnc(r2)
          letter=r1
          IF (letter.EQ.term(2*knt,wh,1))THEN
            r1= -1
            GOTO 599
          END IF
          neg=r3
        ELSE
          rnum=myrand(0)
          conc(knt,wh,1)=cnc(INT(2*rnum)+1)
 500      CONTINUE
          rnum=myrand(0)
          letter=INT(n2*rnum)+1
          IF (letter.EQ.term(2*knt,wh,1)) GOTO 500
          rnum=myrand(0)
          neg=INT(2*rnum)
        END IF
        conc(2*knt+1,wh,1)=3
        term(2*knt+1,wh,1)=letter
        negs(2*knt+1,wh,1)=neg
        pick(2*knt+1,wh,1)=1
599     CONTINUE
      END 

      ! *****************************************************************
      ! *****************************************************************


      ! this subroutine grows a branch
      ! last modification 10/17/02

      SUBROUTINE branch(knt,n2,nkn,ntr,wh,cnc,conc,negs,pick,term,
     #                  r1,r2,r3)
      IMPLICIT NONE

        ! arguments in
          INTEGER knt,n2,nkn,ntr,wh,r1,r2,r3
          INTEGER cnc(3)
          INTEGER conc(nkn,ntr,3)
          INTEGER negs(nkn,ntr,3)
          INTEGER pick(nkn,ntr,3)
          INTEGER term(nkn,ntr,3)
        ! local
          INTEGER letter,neg,loc1,loc2
          REAL rnum,myrand
        ! arguments out

      ! pass old substructure
        loc1=2*knt
        loc2=4*knt
        conc(loc2,wh,1)=3
        term(loc2,wh,1)=term(loc1,wh,1)
        negs(loc2,wh,1)=negs(loc1,wh,1)
        pick(loc2,wh,1)=1
        loc1=2*knt+1
        loc2=4*knt+1
        conc(loc2,wh,1)=3
        term(loc2,wh,1)=term(loc1,wh,1)
        negs(loc2,wh,1)=negs(loc1,wh,1)
        pick(loc2,wh,1)=1
        loc1=knt
        loc2=2*knt
        conc(loc2,wh,1)=conc(loc1,wh,1)
        term(loc2,wh,1)=term(loc1,wh,1)
        negs(loc2,wh,1)=0
      ! create the new leaf 
        conc(2*knt+1,wh,1)=3
        pick(2*knt+1,wh,1)=1
        term(knt,wh,1)=0
        negs(knt,wh,1)=0        
        IF(r1.LT.0)THEN
 500      CONTINUE
          rnum=myrand(0)
          letter=INT(n2*rnum)+1
          IF (letter.EQ.term(4*knt,wh,1).OR.
     #      letter.EQ.term(4*knt+1,wh,1)) GOTO 500
          rnum=myrand(0)
          neg=INT(2*rnum)
      ! create operator in knt
          rnum=myrand(0)
          conc(knt,wh,1)=cnc(INT(2*rnum)+1)
        ELSE
          letter=r1
          IF (letter.EQ.term(4*knt,wh,1).OR.
     #      letter.EQ.term(4*knt+1,wh,1)) THEN
            r1 = -1
            GOTO 592
          END IF
          neg=r2
          conc(knt,wh,1)=cnc(r3)
        END IF
        term(2*knt+1,wh,1)=letter
        negs(2*knt+1,wh,1)=neg
592     CONTINUE
      END

      ! *****************************************************************
      ! *****************************************************************

      ! this subroutine prunes a leaf
      ! last modification 10/17/02

      SUBROUTINE prune(knt,dbl,sng,nkn,ntr,wh,conc,negs,pick,term)
      IMPLICIT NONE

        ! arguments in
          INTEGER knt,nkn,ntr,wh,dbl,sng,loc1,loc2
          INTEGER conc(nkn,ntr,3)
          INTEGER negs(nkn,ntr,3)
          INTEGER pick(nkn,ntr,3)
          INTEGER term(nkn,ntr,3)
        ! arguments out

      ! move operator 
        loc1=dbl
        loc2=knt
        conc(loc2,wh,1)=conc(loc1,wh,1)
        term(loc2,wh,1)=0
        negs(loc2,wh,1)=0
      ! move leaves
        loc1=2*dbl
        loc2=sng
        conc(loc2,wh,1)=conc(loc1,wh,1)
        term(loc2,wh,1)=term(loc1,wh,1)
        negs(loc2,wh,1)=negs(loc1,wh,1)
        pick(loc2,wh,1)=pick(loc1,wh,1)
        conc(loc1,wh,1)=0
        term(loc1,wh,1)=0
        negs(loc1,wh,1)=0
        pick(loc1,wh,1)=0
        loc1=2*dbl+1
        loc2=dbl
        conc(loc2,wh,1)=conc(loc1,wh,1)
        term(loc2,wh,1)=term(loc1,wh,1)
        negs(loc2,wh,1)=negs(loc1,wh,1)
        pick(loc2,wh,1)=pick(loc1,wh,1)
        conc(loc1,wh,1)=0
        term(loc1,wh,1)=0
        negs(loc1,wh,1)=0
        pick(loc1,wh,1)=0

      END 

      ! *****************************************************************
      ! *****************************************************************

      ! this subroutine does the randomization test
      ! last modification 10/17/02

      SUBROUTINE nullmodelx(n1,n2,mdl,nkn,nsp,ntr,conc,negs,pick,term,
     #                     storage,slprbc,datri,weight,
     #                     tstr,tend,tint,ehm,msz,nsep,seps,cnc,
     #                     score,betas,ssize,dcph,ordrs,nfcnt,resp,mtm)
      IMPLICIT NONE

        ! parameters
          INTEGER LGCnsepMAX,LGCn1MAX
          PARAMETER (LGCn1MAX   = 10000)
          PARAMETER (LGCnsepMAX =    50)
        ! arguments in
          INTEGER ehm,n1,n2,mdl,msz,nkn,nsp,ntr,nsep,nfcnt,cnc(3)
          INTEGER dcph(n1),ordrs(n1),mtm
          REAL tstr,tend,tint,weight(n1),slprbc(25)
          INTEGER datri(n2,n1)
          REAL seps(nsep,n1),resp(n1)
        ! local
          REAL rseps(LGCnsepMAX,LGCn1MAX)
          INTEGER xstop
        ! arguments out
          INTEGER ssize,conc(nkn,ntr,3)
          INTEGER negs(nkn,ntr,3)
          INTEGER pick(nkn,ntr,3)
          INTEGER term(nkn,ntr,3)
          INTEGER storage(2*ntr*nkn*n1)
          REAL score(3),betas(0:(nsep+ntr),3)

        xstop=0
        CALL stopper(LGCn1MAX,n1,"LGCn1MAX","nullmodelx()",8,12,xstop,0)
        CALL stopper(LGCnsepMAX,nsep,"LGCnsepMAX","nullmodelx()",10,12,
     #          xstop,1)
        CALL nullmodel(n1,n2,mdl,nkn,nsp,ntr,conc,negs,pick,term,
     #                     storage,slprbc,datri,weight,
     #                     tstr,tend,tint,ehm,msz,nsep,seps,cnc,score,
     #                     betas,ssize,dcph,ordrs,nfcnt,rseps,resp,mtm)

        END

      ! *****************************************************************
      ! *****************************************************************

      ! *****************************************************************
      ! *****************************************************************

      ! this subroutine does the randomization test
      ! last modification 10/17/02

      SUBROUTINE nullmodel(n1,n2,mdl,nkn,nsp,ntr,conc,negs,pick,term,
     #                     storage,slprbc,datri,weight,
     #                     tstr,tend,tint,ehm,msz,nsep,seps,cnc,score,
     #                     betas,ssize,dcph,ordrs,nfcnt,rseps,resp,mtm)
      IMPLICIT NONE

        ! parameters
          INTEGER LGCn1MAX
          PARAMETER (LGCn1MAX   = 10000)
        ! random number generator
          REAL myrand               ! Declare the type of the rand() function
        ! arguments in
          INTEGER ehm,n1,n2,mdl,msz,nkn,nsp,ntr,nsep,nfcnt,cnc(3)
          INTEGER dcph(n1),ordrs(n1),mtm
          REAL tstr,tend,tint,weight(n1),slprbc(25)
          INTEGER datri(n2,n1)
          REAL seps(nsep,n1),resp(n1)
        ! local
          INTEGER i,j,dummy(LGCn1MAX),inums(LGCn1MAX),rdcph(LGCn1MAX)
          INTEGER rnumsi(LGCn1MAX),sordrs(LGCn1MAX)
          REAL datrtoken(LGCn1MAX),rnumsr(LGCn1MAX),rrsp(LGCn1MAX)
          REAL rwgt(LGCn1MAX),rseps(nsep,n1),penalty
        ! arguments out
          INTEGER ssize,conc(nkn,ntr,3)
          INTEGER negs(nkn,ntr,3)
          INTEGER pick(nkn,ntr,3)
          INTEGER term(nkn,ntr,3)
          INTEGER storage(2*ntr*nkn*n1),xstop,mcmc,rd1(2)
          REAL score(3),betas(0:(nsep+ntr),3),rdummy(10),rd2(2),rd3(2)
          INTEGER rd4(2),bout

        mcmc=0
        xstop=0
        CALL stopper(LGCn1MAX,n1,"LGCn1MAX","nullmodel()",8,11,xstop,1)

      ! randomize all cases, save responses and randomize
        DO j=1,n1
          dummy(j)=j
          rnumsi(j)=j
          rnumsr(j)=myrand(0)
        END DO
        CALL psort2(rnumsr,n1,dummy,n1,rnumsi)
        DO j=1,n1
          rrsp(j)=resp(j)
          sordrs(j)=ordrs(j)
        END DO
        DO j=1,n1
          resp(j)=rrsp(rnumsi(j))
          rwgt(j)=weight(rnumsi(j))
          rdcph(j)=dcph(rnumsi(j))
          inums(j)=j
          ordrs(j)=j
          datrtoken(j)=resp(j)
        END DO
        IF(nsep.GT.0)THEN
          DO j=1,n1
            DO i=1,nsep
               rseps(i,j)=seps(i,rnumsi(j))
            END DO
          END DO
        END IF
        CALL psort2(datrtoken,n1,inums,n1,ordrs)
 1000   CONTINUE
        penalty=0.
        CALL annealing(n1,n2,mdl,nkn,nsp,ntr,conc,negs,pick,term,
     #                 storage,slprbc,datri,rwgt,tstr,tend,tint,ehm,
     #                 msz,nsep,rseps,cnc,score,betas,ssize,rdcph,ordrs,
     #                 nfcnt,penalty,resp,mtm,mcmc,rdummy,rd1,rd2,rd3,
     #                 rd4,bout)
        DO j=1,n1
          resp(j)=rrsp(j)
          ordrs(j)=sordrs(j)
        END DO

      END 

      ! *****************************************************************
      ! *****************************************************************


      ! this subroutine does the randomization test
      ! last modification 10/17/02

      SUBROUTINE randomizationx(n1,n2,mdl,nkn,nsp,ntr,conc,negs,pick,
     #                         term,storage,slprbc,datri,weight,tstr,
     #                         tend,tint,ehm,msz,nsep,seps,cnc,score,
     #                         betas,ssize,dcph,ordrs,iotrees,nfcnt,
     #                         ntrup,error,resp,mtm)
      IMPLICIT NONE

        ! parameters
          INTEGER LGCnsepMAX,LGCn1MAX,LGCntrMAX,LGCn2MAX
          PARAMETER (LGCn1MAX   = 10000)
          PARAMETER (LGCn2MAX   =  1000)
          PARAMETER (LGCnsepMAX =    50)
          PARAMETER (LGCntrMAX  =     5)
        ! arguments in
          INTEGER ehm,n1,n2,mdl,msz,nkn,nsp,ntr,nsep,nfcnt
          INTEGER cnc(3),dcph(n1),ordrs(n1),ntrup,mtm,iotrees(1)
          REAL tstr,tend,tint,slprbc(25),weight(n1)
          INTEGER datri(n2,n1)
          REAL seps(nsep,n1),resp(n1)
        ! local
          REAL rseps(LGCnsepMAX,LGCn1MAX)
          INTEGER prtr(LGCn1MAX,LGCntrMAX)
          INTEGER rdatri(LGCn2MAX,LGCn1MAX)
          INTEGER xstop,i,j
        ! arguments out
          INTEGER ssize,conc(nkn,ntr,3)
          INTEGER negs(nkn,ntr,3),error
          INTEGER pick(nkn,ntr,3)
          INTEGER term(nkn,ntr,3)
          INTEGER storage(2*ntr*nkn*n1)
          REAL score(3),betas(0:(nsep+ntr),3)

          xstop=0
          CALL stopper(LGCn1MAX,n1,"LGCn1MAX","randomizationx()",8,16,
     #            xstop,0)
          CALL stopper(LGCn2MAX,n2,"LGCn2MAX","randomizationx()",8,16,
     #            xstop,0)
          CALL stopper(LGCntrMAX,ntr,"LGCntrMAX","randomizationx()",9,
     #            16,xstop,0)
          CALL stopper(LGCnsepMAX,nsep,"LGCnsepMAX","randomizationx()",
     #            10,16,xstop,1)
        DO i=1,LGCn1MAX
           DO j=1,LGCntrMAX
              prtr(i,j)=0
           END DO
        END DO
            CALL randomization(n1,n2,mdl,nkn,nsp,ntr,conc,negs,pick,
     #                         term,storage,slprbc,datri,weight,tstr,
     #                         tend,tint,ehm,msz,nsep,seps,cnc,score,
     #                         betas,ssize,dcph,ordrs,iotrees,nfcnt,
     #                         ntrup,error,rseps,prtr,rdatri,resp,mtm)

        END

      ! *****************************************************************
      ! *****************************************************************


      ! this subroutine does the randomization test
      ! last modification 10/17/02

      SUBROUTINE randomization(n1,n2,mdl,nkn,nsp,ntr,conc,negs,pick,
     #                         term,storage,slprbc,datri,weight,tstr,
     #                         tend,tint,ehm,msz,nsep,seps,cnc,score,
     #                         betas,ssize,dcph,ordrs,iotrees,nfcnt,
     #                         ntrup,error,rseps,prtr,rdatri,resp,mtm)
      IMPLICIT NONE

        ! parameters
          INTEGER LGCn1MAX,LGCntrMAX
          PARAMETER (LGCn1MAX   = 10000)
          PARAMETER (LGCntrMAX  =     5)
        ! arguments in
          INTEGER ehm,n1,n2,mdl,msz,nkn,nsp,ntr,nsep,iotrees(1),nfcnt
          INTEGER cnc(3),dcph(n1),ordrs(n1),ntrup,mtm
          REAL tstr,tend,tint,slprbc(25),weight(n1)
          INTEGER datri(n2,n1)
          REAL seps(nsep,n1),resp(n1)
        ! local
          INTEGER j,k,l,wh,ncl,msznew,rdcp(LGCn1MAX)
          INTEGER prtr(n1,ntr),xstop,rordrs(LGCn1MAX)
          INTEGER nprdcl(2**LGCntrMAX),prdcl(LGCn1MAX,2**LGCntrMAX)
          REAL rwgt(LGCn1MAX),rseps(nsep,n1)
          INTEGER rdatri(n2,n1)
          REAL penalty,rresp(LGCn1MAX)
        ! arguments out
          INTEGER ssize,conc(nkn,ntr,3)
          INTEGER negs(nkn,ntr,3),error
          INTEGER pick(nkn,ntr,3)
          INTEGER term(nkn,ntr,3)
          INTEGER storage(2*ntr*nkn*n1),mcmc,rd1(2)
          REAL score(3),betas(0:(nsep+ntr),3),dummy(10),rd2(2),rd3(2)
          INTEGER rd4(2),bout

          xstop=0
          mcmc=0
          CALL stopper(LGCn1MAX,n1,"LGCn1MAX","randomization()",8,15,
     #                 xstop,0)
          CALL stopper(LGCntrMAX,ntr,"LGCntrMAX","randomization()",9,15,
     #            xstop,1)

      ! read and evaluate
        DO j=1,nkn
          DO k=1,ntr
            DO l=1,3
              conc(j,k,l)=0
              negs(j,k,l)=0
              pick(j,k,l)=0
              term(j,k,l)=0
            END DO
          END DO
        END DO

        DO wh=1,ntr
          CALL read_treex(wh,1,msz,nkn,ntr,conc,term,negs,pick,
     #                    iotrees,error)
          IF(error.eq.1) RETURN
        END DO
        DO wh=1,ntr
          CALL evaluate_first(wh,n1,n2,nkn,ntr,conc,term,negs,pick,
     #                        datri,prtr)
        END DO
        DO j=1,nkn
          DO k=1,ntr
            DO l=1,3
              conc(j,k,l)=0
              negs(j,k,l)=0
              pick(j,k,l)=0
              term(j,k,l)=0
            END DO
          END DO
        END DO

      ! reset to unlimited size
        msznew=ntrup*nkn

      ! pass on
        DO j=1,n1
          rwgt(j)=weight(j)
          rdcp(j)=dcph(j)
          rordrs(j)=ordrs(j)
          DO k=1,n2
            rdatri(k,j)=datri(k,j)
          END DO
          rresp(j)=resp(j)
          DO k=1,nsep
            rseps(k,j)=seps(k,j)
          END DO
        END DO
        ncl=2**ntr
        penalty=0
         

      ! randomization
        CALL ident_prdcl(n1,n2,ntr,prtr,ncl,nprdcl,prdcl)
        CALL rand_prdcl(n1,n2,nsep,rresp,rwgt,rseps,
     #                  ncl,nprdcl,prdcl,rdcp,rordrs)
        CALL annealing(n1,n2,mdl,nkn,nsp,ntrup,conc,negs,pick,term,
     #                storage,slprbc,rdatri,rwgt,tstr,tend,tint,ehm,
     #                msznew,nsep,rseps,cnc,score,betas,ssize,rdcp,
     #                rordrs,nfcnt,penalty,rresp,mtm,mcmc,dummy,rd1,rd2,
     #                rd3,rd4,bout)

      END 

      ! *****************************************************************
      ! *****************************************************************

      ! this subroutine identifies the prediction classes
      ! last modification 10/17/02

      SUBROUTINE ident_prdcl(n1,n2,ntr,prtr,ncl,nprdcl,prdcl)
      IMPLICIT NONE

        ! arguments in
          INTEGER n1,n2,ntr
          INTEGER prtr(n1,ntr)
        ! local
          INTEGER j,k
          INTEGER cl
        ! arguments out
          INTEGER ncl
          INTEGER nprdcl(ncl)
          INTEGER prdcl(n1,ncl)

      ! initialize and identify the prediction classes

        DO j=1,ncl
          nprdcl(j)=0
          DO k=1,n1
            prdcl(k,j)=0
          END DO
        END DO
        DO k=1,n1
          cl=1
          DO j=1,ntr
            cl=cl+2**(j-1)*prtr(k,j)
          END DO
          nprdcl(cl)=nprdcl(cl)+1
          prdcl(nprdcl(cl),cl)=k
        END DO

      END 

      ! *****************************************************************
      ! *****************************************************************

      ! this subroutine randomizes the prediction classes
      ! last modification 10/17/02

      SUBROUTINE rand_prdcl(n1,n2,nsep,resp,rwgt,rseps,
     #                      ncl,nprdcl,prdcl,rdcp,ordrs)
      IMPLICIT NONE

        ! parameters
          INTEGER LGCn1MAX
          PARAMETER (LGCn1MAX   = 10000)
        ! random number generator
          REAL myrand               ! declare the type of the rand() function
        ! arguments in
          INTEGER n1,n2,ncl,nsep
          INTEGER rdcp(n1),ordrs(n1)
          INTEGER nprdcl(ncl)
          INTEGER prdcl(n1,ncl)
          REAL rwgt(n1),rseps(nsep,n1),resp(n1)
        ! local
          INTEGER j,k,l
          INTEGER nn,xstop
          INTEGER rnumsi(LGCn1MAX),wk2(LGCn1MAX)
          REAL wk1(LGCn1MAX)
        ! arguments out
        xstop=0
        CALL stopper(LGCn1MAX,n1,"LGCn1MAX","rand_prdcl()",8,12,xstop,1)

      ! randomize the prediction classes
        DO j=1,ncl
          nn=nprdcl(j)
          IF (nn.GT.0) THEN
            DO k=1,nn
              wk2(k)=k
              rnumsi(k)=k
              wk1(k)=myrand(0) 
            END DO
            CALL psort2(wk1,nn,wk2,nn,rnumsi)
            DO k=1,nn
              wk1(k)=resp(prdcl(k,j))
            END DO
            DO k=1,nn
              resp(prdcl(k,j))=wk1(rnumsi(k))
            END DO 
            DO k=1,nn
              wk1(k)=rwgt(prdcl(k,j))
              wk2(k)=rdcp(prdcl(k,j))
            END DO
            DO k=1,nn
              rwgt(prdcl(k,j))=wk1(rnumsi(k))
              rdcp(prdcl(k,j))=wk2(rnumsi(k))
            END DO 
            IF (nsep.GT.0) THEN
              DO l=1,nsep
                DO k=1,nn
                  wk1(k)=rseps(l,prdcl(k,j))
                END DO
                DO k=1,nn
                  rseps(l,prdcl(k,j))=wk1(rnumsi(k))
                END DO
              END DO
            END IF 
          END IF
        END DO

      ! calculate order
        DO j=1,n1
          wk2(j)=j
          ordrs(j)=j
          wk1(j)=resp(j)
        END DO
        CALL psort2(wk1,n1,wk2,n1,ordrs)

      END 

      ! *****************************************************************
      ! *****************************************************************

      ! this subroutine reads a tree
      ! last modification 10/17/02

      SUBROUTINE read_treex(wh,stst,msz,nkn,ntr,conc,term,negs,pick,
     #                      iotrees,error)
      IMPLICIT NONE

        ! arguments in
          INTEGER msz,nkn,ntr,wh,stst,iotrees(1)
        ! local
          INTEGER k,ltree,j
        ! arguments out
          INTEGER conc(nkn,ntr,3)
          INTEGER negs(nkn,ntr,3)
          INTEGER error
          INTEGER pick(nkn,ntr,3)
          INTEGER term(nkn,ntr,3)

        error=-1
        ltree=4*nkn+3
        DO k=0,1000
           IF(error.LT.0)THEN
              IF(iotrees(k*ltree+1).LT.0)THEN
                 error=1
              ELSE
                IF(iotrees(k*ltree+1).EQ.msz.AND.
     #             iotrees(k*ltree+2).EQ.ntr.AND.
     #             iotrees(k*ltree+3).EQ.wh)THEN
                  error=0
                  DO j=1,nkn
                    conc(j,wh,stst)=iotrees(k*ltree+(j-1)*4+4)
                    term(j,wh,stst)=iotrees(k*ltree+(j-1)*4+5)
                    negs(j,wh,stst)=iotrees(k*ltree+(j-1)*4+6)
                    pick(j,wh,stst)=iotrees(k*ltree+(j-1)*4+7)
                  END DO
               END IF
             END IF
           END IF
        END DO
        IF(error.EQ.-1)error=1

      END 

      ! *****************************************************************
      ! *****************************************************************

      ! this subroutine keeps track of the scores, trees etc.
      ! last modification 10/17/02

      SUBROUTINE recording(accept,wh,nkn,ntr,nsep,score,betas,
     #                     conc,negs,pick,term,mcmc)
      IMPLICIT NONE

        ! arguments in
          INTEGER accept,nkn,nsep,ntr,wh,mcmc
          INTEGER conc(nkn,ntr,3)
          INTEGER negs(nkn,ntr,3)
          INTEGER pick(nkn,ntr,3)
          INTEGER term(nkn,ntr,3)
          REAL score(3)
          REAL betas(3,0:(nkn+ntr))
        ! local
          INTEGER j,k
        ! arguments out

      ! record trees
        IF (score(1).LT.score(3)) THEN
          IF(accept.GT.0)THEN
          CALL copytree(ntr,nkn,conc,negs,pick,term,-1,1,3)
          DO k=0,(nsep+ntr)
            betas(3,k)=betas(1,k)
          END DO
          END IF
        END IF
        IF(mcmc.GE.1)THEN
          CALL copytree(ntr,nkn,conc,negs,pick,term,-1,2,3)
        END IF
        IF (accept.EQ.1) THEN
          CALL copytree(ntr,nkn,conc,negs,pick,term,wh,1,2)
          DO k=0,(nsep+ntr)
            betas(2,k)=betas(1,k)
          END DO
        ELSE
          CALL copytree(ntr,nkn,conc,negs,pick,term,wh,2,1)
          DO k=0,(nsep+ntr)
            betas(1,k)=betas(2,k)
          END DO
        END IF

      ! record scores
        IF (score(1).LT.score(3).and.accept.eq.1) score(3)=score(1)
        IF (accept.EQ.1) THEN
          score(2)=score(1)
        ELSE
          score(1)=score(2)
        END IF

      END 

      ! *****************************************************************
      ! *****************************************************************

      ! storing or restoring data after mtp 0
      ! last modification 10/17/02

      SUBROUTINE restoring(accept,wh,n1,nkn,ntr,storage,nwkv,wkv)
      IMPLICIT NONE

        ! arguments in
          INTEGER accept,n1,nkn,ntr,nwkv,wh
          INTEGER wkv(nkn)
        ! local
          INTEGER j,k,wkvj,pp,pp2
        ! arguments out
          INTEGER storage(2*ntr*nkn*n1)

      ! storing/restoring
        pp2=ntr*nkn*n1
        IF (accept.EQ.1) THEN
          DO j=1,nwkv
            wkvj=wkv(j)
            pp=(wh-1)*nkn*n1+(wkvj-1)*n1
            DO k=1,n1
              storage(pp2+pp+k)=storage(pp+k)
            END DO
          END DO
        ELSE
          DO j=1,nwkv
            wkvj=wkv(j)
            pp=(wh-1)*nkn*n1+(wkvj-1)*n1
            DO k=1,n1
              storage(pp+k)=storage(pp2+pp+k)
            END DO
          END DO
        END IF

      END

      ! *****************************************************************
      ! *****************************************************************

      ! this subroutine calculates the score of the model under
      ! consideration
      ! last modification 10/17/02

      SUBROUTINE scoring(prtr,rsp,dcph,ordrs,weight,n1,ntr,mdl,nop,wh,
     #                   nsep,seps,score,betas,reject,xtxsep,mtm,nopold)
      IMPLICIT NONE

        ! parameters
          INTEGER LGCbetaMAX,mtm
          PARAMETER (LGCbetaMAX = 55)
        ! arguments in
          INTEGER mdl,n1,nop,nsep,ntr,wh
          INTEGER dcph(n1),ordrs(n1),nopold
          INTEGER prtr(n1,ntr)
          REAL rsp(n1),weight(n1)
          REAL seps(nsep,n1),xtxsep(0:nsep,0:nsep)
        ! local
          INTEGER j,oops,xstop
          REAL smbetas(0:LGCbetaMAX)
        ! arguments out
          INTEGER reject
          REAL score(3)
          REAL betas(0:(nsep+ntr))
          character *100 astring

        xstop=0
        CALL stopper(LGCbetaMAX,nsep+ntr,"LGCbetaMAX","scoring()",10,9,
     #          xstop,1)


      ! initialize
        DO j=0,nsep+ntr
          betas(j)=0.0
        END DO
        reject=0

      ! choose scoring function for model type
        IF (mdl.NE.1.AND.nopold.LE.nop)
     #    CALL singularities(n1,nop,ntr,wh,prtr,nsep,seps,reject,mtm)
        IF (reject.EQ.0) THEN
          IF (mdl.EQ.0)THEN
            CALL My_own_fitting(prtr,rsp,dcph,ordrs,weight,n1,ntr,
     #                        nop,wh,nsep,seps,score(1),smbetas,reject)
            DO j=0,(nsep+ntr)
              betas(j)=smbetas(j)
            END DO 
          ELSE IF (mdl.EQ.5)THEN
            CALL expofit(prtr,rsp,dcph,weight,n1,ntr,
     #                        nop,nsep,seps,score(1),smbetas,reject)
            DO j=0,(nsep+ntr)
              betas(j)=smbetas(j)
            END DO 
          ELSE IF (mdl.EQ.1) THEN
            score(1)=0.0
            DO j=1,n1
              score(1)=score(1)+weight(j)*(REAL(prtr(j,1))-rsp(j))**2.0
            END DO
          ELSE IF (mdl.EQ.2) THEN
            CALL calcbetarss(n1,nop,ntr,prtr,nsep,seps,rsp,weight,
     #                       smbetas,oops,xtxsep)
            IF(oops.EQ.1)THEN
               reject=1   
            ELSE
              CALL calcrss(nop,n1,ntr,smbetas,prtr,nsep,seps,rsp,
     #                     weight,score)
              DO j=0,(nsep+ntr)
                betas(j)=smbetas(j)
              END DO 
            END IF
          ELSE IF (mdl.EQ.3) THEN
            CALL calcdev(n1,nop,ntr,prtr,nsep,seps,rsp,weight,
     #                   betas,score,reject)
          ELSE IF (mdl.EQ.4) THEN
            CALL calcplcph(nop,n1,ntr,betas,prtr,nsep,seps,rsp,weight,
     #                   dcph,ordrs,score,oops)
          END IF
        END IF

      END 

      ! *****************************************************************
      ! *****************************************************************

      ! this subroutine checks for singularities in the proposed model
      ! last modification 10/16/02

      SUBROUTINE singularities(n1,nop,ntr,wh,prtr,nsep,seps,reject,mtm)
      IMPLICIT NONE

        ! arguments in
          INTEGER n1,nop,nsep,ntr,wh,prtr(n1,ntr),mtm
          REAL seps(nsep,n1)
        ! local
          INTEGER hmuch,sum1,SUM2I,j,k,l,m
          REAL nr1
        ! arguments out
          INTEGER reject,n4

      ! check for singularities in the predicted trees 
        reject=0
        n4=0
        sum1=0
        IF (nop.GT.0) THEN
          sum1=SUM2I(prtr,n1,ntr,1,wh,1,n1)
          nr1=REAL(n1)
          n4=INT(0.05*nr1)
          if(n4.gt.15)n4=15
          if(mtm.gt.0)n4=mtm
          IF (sum1.LT.n4.OR.sum1.GT.(n1-n4)) THEN
            reject=1
          END IF
        END IF

      ! check for singularities with other predictors 
        IF (reject.EQ.0.AND.nop.GT.1) THEN
          DO j=1,nop
            IF (j.NE.wh) THEN
              hmuch=0
              l=1
              m=-1
              IF(prtr(1,wh).EQ.prtr(1,j))l=0
              IF(prtr(1,wh).EQ.prtr(1,j))m= 1
              DO k=1,n1
                IF(prtr(k,wh).NE.(l+m*prtr(k,j)))GOTO 4000
              END DO 
              reject=1
              GOTO 1000
            END IF
 4000       CONTINUE
          END DO
        END IF
 1000   CONTINUE

        IF (reject.EQ.0.AND.nop.GE.1.AND.nsep.GT.0) THEN
          DO j=1,nsep
            DO k=1,n1
              IF((seps(j,k).NE.0.0).AND.(seps(j,k).NE.1.0))GOTO 3000
            END DO 
            l=1
            m=-1
            IF(prtr(1,wh).EQ.seps(j,1))l=0
            IF(prtr(1,wh).EQ.seps(j,1))m= 1
            DO k=1,n1
              IF(prtr(k,wh).NE.(l+m*seps(j,k))) GOTO 3000
            END DO 
            reject=1
            GOTO 2000
 3000       CONTINUE
          END DO
        END IF
 2000   CONTINUE
      END 

      ! *****************************************************************
      ! *****************************************************************

      ! this subroutine estimates the parameters in the linear 
      ! regression model
      ! last modification 10/17/02

      SUBROUTINE calcbetarss(n1,nop,ntr,prtr,nsep,seps,rsp,weight,betas,
     #                       oops,xtxsep)
      IMPLICIT NONE

        ! parameters
          INTEGER LGCbetaMAX
          PARAMETER (LGCbetaMAX =    55)
        ! arguments in
          INTEGER n1,nop,nsep,ntr
          INTEGER prtr(n1,ntr),oops
          REAL rsp(n1),weight(n1)
          REAL seps(nsep,n1),xtxsep(0:nsep,0:nsep)
        ! local
          INTEGER i,j,k,l,nj,xstop
          INTEGER*4 dummy_vec(LGCbetaMAX+1)
          DOUBLE PRECISION dummy_scal
          DOUBLE PRECISION xtx(0:LGCbetaMAX,0:LGCbetaMAX)
          DOUBLE PRECISION xty(0:LGCbetaMAX)
          REAL tweight,tmp
        ! arguments out
          REAL betas(0:(nsep+ntr))

        xstop=0
        CALL stopper(LGCbetaMAX,nsep+ntr,"LGCbetaMAX","calcbetarss()",
     #          10,13,xstop,1)

        ! estimate the parameters for linear regression 
        xtx(0,0)=xtxsep(0,0)
        xty(0)=0.0
        DO k=1,n1
          xty(0)=xty(0)+rsp(k)*weight(k)
        END DO
        IF (nsep.GT.0) THEN
          DO j=1,nsep
            xtx(0,j)=xtxsep(0,j)
            DO k=1,nsep
             xtx(j,k)=xtxsep(j,k)
            END DO        
            xty(j)=0.0
            DO i=1,n1
              tmp=weight(i)*seps(j,i)
              xty(j)=xty(j)+rsp(i)*tmp
            END DO        
            xtx(j,0)=xtx(0,j)
          END DO
        END IF
        IF (nop.GT.0) THEN
          DO j=1,nop
            nj=nsep+j
            xtx(0,nj)=0.0
            xty(nj)=0.0
            DO k=1,nop
              xtx(nj,nsep+k)=0.0
            END DO
            DO k=1,nsep
              xtx(nj,k)=0
            END DO
            DO i=1,n1
              IF(prtr(i,j).EQ.1)THEN
                xtx(0,nj)=xtx(0,nj)+weight(i)
                xty(nj)=xty(nj)+weight(i)*rsp(i)
                DO k=1,j
                  xtx(nj,nsep+k)=xtx(nj,nsep+k)+
     #                                weight(i)*REAL(prtr(i,k))
                END DO        
                DO k=1,nsep
                  xtx(nj,k)=xtx(nj,k)+weight(i)*seps(k,i)
                END DO
              END IF
            END DO        
            xtx(nj,0)=xtx(0,nj)
            DO k=1,nsep
              xtx(k,nj)=xtx(nj,k)
            END DO        
          END DO        
          DO j=1,nop
            DO k=j+1,nop
              xtx(nsep+j,nsep+k)=xtx(nsep+k,nsep+j)
            END DO        
          END DO        
        END IF
        tweight=xtx(0,0)
        DO k=0,(nsep+nop)
          IF(xty(k).lt.1.0D-10)xty(k)=0
          xty(k)=xty(k)/tweight
          DO l=0,(nsep+nop)
            IF(xtx(k,l).lt.1.0D-10)xtx(k,l)=0
            xtx(k,l)=xtx(k,l)/tweight
          END DO        
        END DO
        DO j=1,(nsep+nop+1)
          dummy_vec(j)=0
        END DO
        dummy_scal=0.0
        CALL LUDCMP(xtx,(nsep+nop+1),dummy_vec,dummy_scal,oops,
     #              LGCbetaMAX+1)
        IF(oops.EQ.1)THEN
           RETURN
        END IF
        CALL LUBKSB(xtx,(nsep+nop+1),dummy_vec,xty,LGCbetaMAX+1)
        DO j=0,(nsep+nop)
          betas(j)=xty(j)
        END DO

      END 

      ! *****************************************************************
      ! *****************************************************************

      ! this subroutine calculates the residual standard error for
      ! the linear regression model
      ! last modification 10/17/02

      SUBROUTINE calcrss(nop,n1,ntr,betas,prtr,nsep,seps,rsp,weight,
     #                   score)
      IMPLICIT NONE

        ! arguments in
          INTEGER n1,nop,nsep,ntr
          INTEGER prtr(n1,ntr)
          REAL betas(0:(nsep+ntr))
          REAL rsp(n1),weight(n1)
          REAL seps(nsep,n1)
        ! local
          INTEGER i,j
          REAL prediction
        ! arguments out
          REAL score(3)

      ! calculate residual standard error  
        score(1)=0.0
        DO i=1,n1
          prediction=betas(0)
          IF (nsep.GT.0) THEN
            DO j=1,nsep
              prediction=prediction+betas(j)*seps(j,i)
            END DO
          END IF
          IF (nop.GT.0) THEN
            DO j=1,nop
              prediction=prediction+betas(nsep+j)*REAL(prtr(i,j))
            END DO
          END IF
          score(1)=score(1)+weight(i)*((prediction-rsp(i))**2.0)
        END DO
        score(1)=SQRT(score(1)/REAL(n1-1-nsep-nop))

      END 

      ! *****************************************************************
      ! *****************************************************************

      ! this subroutine estimates the parameters in the logistic 
      ! regression model and calculates the deviance
      ! last modification 10/17/02

      SUBROUTINE calcdev(n1,nop,ntr,prtr,nsep,seps,rsp,weight,betas,
     #                   score,reject)
      IMPLICIT NONE

        ! parameters
          INTEGER LGCbetaMAX,LGCn1MAX
          PARAMETER (LGCn1MAX   = 10000)
          PARAMETER (LGCbetaMAX =    55)
        ! arguments in
          INTEGER n1,nop,nsep,ntr,prtr(n1,ntr)
          REAL rsp(n1),weight(n1),seps(nsep,n1)
        ! local
          INTEGER j,k,conv,iter,info,job,ipvt(LGCbetaMAX+1)
          INTEGER xstop,n1tmp,nu
          REAL eps
          DOUBLE PRECISION dldmu(LGCn1MAX),eta(LGCn1MAX),n(LGCn1MAX)
          DOUBLE PRECISION p(LGCn1MAX),w(LGCn1MAX),winv(LGCn1MAX)
          DOUBLE PRECISION beta(LGCbetaMAX+1),beta_new(LGCbetaMAX+1),det
          DOUBLE PRECISION work(LGCbetaMAX+1),sco(LGCbetaMAX+1)
          DOUBLE PRECISION y(LGCn1MAX),x(LGCn1MAX,LGCbetaMAX+1)
          DOUBLE PRECISION dmudb(LGCn1MAX,LGCbetaMAX+1),loglik
          DOUBLE PRECISION d2mat(LGCbetaMAX+1,LGCbetaMAX+1),loglik_old
          DOUBLE PRECISION prb,mylog,prb2
        ! arguments out
          INTEGER reject
          REAL betas(0:(nsep+ntr)),score(3)

        xstop=0
        CALL stopper(LGCn1MAX,n1,"LGCn1MAX","calcdev()",8,9,xstop,0)
        CALL stopper(LGCbetaMAX,nsep+ntr,"LGCbetaMAX","calcdev()",10,9,
     #          xstop,1)
          
        n1tmp=n1
        IF(nsep+ntr.LT.10)THEN
        CALL redater(nu,x,y,n,prtr,seps,rsp,weight,n1,nop,nsep,LGCn1MAX)
        IF(nu.NE.0)n1tmp=nu
        END IF

      ! initialize
        IF(nu.EQ.0)THEN
        DO j=1,n1
          x(j,1)=1.0
        END DO
        IF (nsep.GT.0) THEN
          DO j=1,n1
            DO k=1,nsep
              x(j,k+1)=seps(k,j)
            END DO
          END DO
        END IF 
        IF (nop.GT.0) THEN
          DO j=1,n1
            DO k=1,nop
              x(j,k+nsep+1)=REAL(prtr(j,k))
            END DO
          END DO
        END IF 
        DO j=1,n1
          y(j)=rsp(j)*weight(j)
          n(j)=weight(j)
        END DO
        END IF
        job=1
        iter=0
        conv=0
        reject=0
        eps=0.000001
        loglik_old=(-100000.0)
        DO j=1,(nsep+ntr+1)
          beta(j)=0.D0
        END DO

      ! estimate the parameters for logistic regression 
        DO WHILE (iter.LT.20.AND.conv.EQ.0)
          iter=iter+1
          CALL lgtderiv(n1tmp,nsep+nop+1,n,x,y,beta,sco,d2mat,eta,p,w,
     #                  winv,dldmu,dmudb,loglik,LGCn1MAX,LGCbetaMAX+1)
          CALL dgefa(d2mat,LGCbetaMAX+1,nsep+nop+1,ipvt,info)
          IF(info.GT.0)THEN
            reject=1
            GOTO 1000
          END IF
          CALL dgedi(d2mat,LGCbetaMAX+1,nsep+nop+1,ipvt,det,work,job)
          DO j=1,nsep+nop+1
            beta_new(j)=beta(j)
            DO k=1,nsep+nop+1
              beta_new(j)=beta_new(j)+d2mat(j,k)*sco(k)
            END DO
            conv=1
            IF (ABS(beta_new(j)-beta(j)).GT.eps*ABS(beta(j))) THEN
              conv=0
            END IF
            IF (ABS(loglik_old-loglik).GT.eps) conv=0
            beta(j)=beta_new(j)
            betas(j-1)=beta_new(j)
          END DO
          loglik_old=loglik
        END DO
        IF (conv.EQ.0) THEN
          reject=1
          GOTO 1000
        END IF

      ! calculate the deviance
        score(1)=0.0
        DO k=1,n1tmp
          IF (p(k).GT.0.D0.AND.p(k).LT.1.D0) THEN
            prb=p(k)
            prb2=1.-prb
            score(1)=score(1)-2.0*(y(k)*mylog(prb)+
     #                             (n(k)-y(k))*mylog(prb2))
          ELSE
            reject=1
            GOTO 1000
          ENDIF
        ENDDO

 1000   CONTINUE

      END 

      ! *****************************************************************
      ! *****************************************************************


      ! this subroutine calculates the partial likelihood for the
      ! cox proportional hazards model
      ! last modification 10/17/02

      SUBROUTINE calcplcph(nop,n1,ntr,betas,prtr,nsep,seps,rsp,weight,
     #                     dcph,ordrs,score,oops)
      IMPLICIT NONE

        ! parameters
          INTEGER LGCbetaMAX,LGCn1MAX
          PARAMETER (LGCn1MAX   = 10000)
          PARAMETER (LGCbetaMAX =    55)
        ! arguments in
          INTEGER n1,nop,nsep,ntr,dcph(n1),ordrs(n1),prtr(n1,ntr)
          REAL rsp(n1),weight(n1),seps(nsep,n1)
        ! local
          INTEGER i,j,k,nnf(2),xstop
          DOUBLE PRECISION loglf,betaf(LGCbetaMAX)
          DOUBLE PRECISION covsf(LGCn1MAX*LGCbetaMAX)
        ! arguments out
          REAL score(3),betas(0:(nsep+ntr)),r
          INTEGER oops
          CHARACTER *80 astring

        xstop=0
        CALL stopper(LGCn1MAX,n1,"LGCn1MAX","calcplcph()",8,11,xstop,0)
        CALL stopper(LGCbetaMAX,nsep+ntr,"LGCbetaMAX","calcplcph()",10,
     #          11,xstop,1)

        DO i=1,n1
          IF((dcph(i).NE.0).and.(dcph(i).NE.1)) THEN
            astring(1:15)="censoring case "
            CALL makeistring(16,23,astring,i,8)
            astring(23:42)="not 0 or 1 -- sorry"
            CALL stringprint(astring,42)
            STOP
          END IF
        END DO

      ! calculate predictors
        nnf(1)=nop+nsep
        nnf(2)=n1
        DO j=1,(n1*(nsep+ntr))
          covsf(j)=0.D0
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
        
      ! calculate partial likelihood
        CALL myphxx(dcph,ordrs,covsf,nnf(1),n1,nsep,ntr,loglf,betaf,
     #              oops,weight)
        r=n1
        score(1)=-loglf
        betas(0)=0.0
        DO j=1,(nsep+nop)
          betas(j)=betaf(j)
        END DO

      END 

      ! *****************************************************************
      ! *****************************************************************



      ! this subroutine stores some variables needed in other routines
      ! last modification 10/17/02

      SUBROUTINE storing(nkn,ntr,conc,pick,npckmv,pickmv,ssize,nop)
      IMPLICIT NONE

        ! parameters 
        ! arguments in
          INTEGER nkn,ntr
          INTEGER conc(nkn,ntr,3)
          INTEGER pick(nkn,ntr,3)
        ! local
          INTEGER i,j,sibling
        ! arguments out
          INTEGER nop,ssize
          INTEGER npckmv(6,ntr)
          INTEGER pickmv(6,nkn,ntr)

      ! determine which moves we can make
        ssize=0
        nop=0
        DO j=1,ntr
          DO i=1,6
            npckmv(i,j)=0
          END DO
          DO i=1,nkn
            IF (pick(i,j,1).EQ.1) THEN
              nop=j
              IF (conc(i,j,1).EQ.3) THEN
                ssize=ssize+1
                npckmv(1,j)=npckmv(1,j)+1
                pickmv(1,npckmv(1,j),j)=i
                IF(i.EQ.1)THEN
                   npckmv(3,j)=npckmv(3,j)+1
                   pickmv(3,npckmv(3,j),j)=i
                ELSE
                   IF (MOD(i,2).EQ.0) THEN
                     sibling=i+1
                   ELSE
                      sibling=i-1
                   END IF
                   IF (conc(sibling,j,1).EQ.3) THEN
                     npckmv(3,j)=npckmv(3,j)+1
                     pickmv(3,npckmv(3,j),j)=i
                   END IF
                END IF
                IF(2*i.LE.nkn) THEN
                   npckmv(4,j)=npckmv(4,j)+1
                   pickmv(4,npckmv(4,j),j)=i
                END IF
              ELSE
                npckmv(2,j)=npckmv(2,j)+1
                pickmv(2,npckmv(2,j),j)=i
                IF(4*i.LE.nkn)THEN
                 IF((conc(2*i,j,1).EQ.3).AND.(conc(2*i+1,j,1).EQ.3))THEN
                     npckmv(5,j)=npckmv(5,j)+1
                     pickmv(5,npckmv(5,j),j)=i
                  END IF
                  IF(conc(2*i,j,1).EQ.3.AND.conc(4*i+2,j,1).EQ.3.AND.
     #               conc(4*i+3,j,1).EQ.3) THEN
                     npckmv(6,j)=npckmv(6,j)+1
                     pickmv(6,npckmv(6,j),j)=i
                  ELSE IF (conc(2*i+1,j,1).EQ.3.AND.conc(4*i,j,1).EQ.3
     #                      .AND.conc(4*i+1,j,1).EQ.3) THEN
                     npckmv(6,j)=npckmv(6,j)+1
                     pickmv(6,npckmv(6,j),j)= -i
                  END IF
                END IF
              END IF
            END IF
          END DO
        END DO
      END 

      ! *****************************************************************
      ! *****************************************************************
      ! this subroutine scores a test set for a given model
      ! last modification 10/17/02

      SUBROUTINE testsetx(n1,n2,mdl,nkn,ntr,conc,negs,pick,term,betas,
     #                   datri,weight,dcph,ordrs,nsep,seps,score,resp)
      IMPLICIT NONE

        ! parameters 
          INTEGER LGCnknMAX,LGCntrMAX,LGCn1MAX,LGCbetaMAX
          PARAMETER (LGCn1MAX   = 10000)
          PARAMETER (LGCnknMAX  =   128)
          PARAMETER (LGCntrMAX  =     5)
          PARAMETER (LGCbetaMAX =    55)
        ! arguments in
          INTEGER n1,n2,mdl,nkn,ntr,nsep,dcph(n1),ordrs(n1)
          INTEGER conc(nkn,ntr,3)
          INTEGER negs(nkn,ntr,3)
          INTEGER pick(nkn,ntr,3)
          INTEGER term(nkn,ntr,3)
          REAL betas(3,0:(nsep+ntr)),weight(n1)
          INTEGER datri(n2,n1)
          REAL seps(nsep,n1),resp(n1)
        ! local
          INTEGER pickmv(6,LGCnknMAX,LGCntrMAX),xstop
          INTEGER prtr(LGCn1MAX,LGCntrMAX)
          INTEGER npckmv(6,LGCntrMAX),i,j
          REAL rsp(LGCn1MAX),smbetas(0:LGCbetaMAX)
        ! arguments out
          REAL score(3)
        xstop=0
        CALL stopper(LGCn1MAX,n1,"LGCn1MAX","testsetx()",8,10,xstop,0)
        CALL stopper(LGCnknMAX,nkn,"LGCnknMAX","testsetx()",9,10,
     #               xstop,0)
        CALL stopper(LGCntrMAX,ntr,"LGCntrMAX","testsetx()",9,10,
     #               xstop,0)
        CALL stopper(LGCbetaMAX,ntr+nsep,"LGCbetaMAX","testsetx()",10,
     #               10,xstop,1)
        DO i=1,LGCn1MAX
           DO j=1,LGCntrMAX
              prtr(i,j)=0
           END DO
        END DO
        DO i=0,(nsep+ntr)
           smbetas(i)=betas(3,i)
        END DO
        CALL testset(n1,n2,mdl,nkn,ntr,conc,negs,pick,term,smbetas,
     #                   datri,weight,dcph,ordrs,nsep,seps,score,
     #                   pickmv,prtr,rsp,npckmv,resp)

        END

      ! *****************************************************************
      ! *****************************************************************
      ! this subroutine scores a test set for a given model
      ! last modification 10/17/02

      SUBROUTINE testset(n1,n2,mdl,nkn,ntr,conc,negs,pick,term,betas,
     #                   datri,weight,dcph,ordrs,nsep,seps,score,
     #                   pickmv,prtr,rsp,npckmv,resp)
      IMPLICIT NONE

        ! arguments in
          INTEGER n1,n2,mdl,nkn,ntr,nsep,dcph(n1),ordrs(n1)
          INTEGER conc(nkn,ntr,3)
          INTEGER negs(nkn,ntr,3)
          INTEGER pick(nkn,ntr,3)
          INTEGER term(nkn,ntr,3)
          REAL betas(0:(nsep+ntr)),weight(n1)
          INTEGER datri(n2,n1)
          REAL seps(nsep,n1),resp(n1)
        ! local
          INTEGER wh,i,j,k,nop,ssize,npckmv(6,ntr)
          INTEGER pickmv(6,nkn,ntr),prtr(n1,ntr)
          REAL rsp(n1)
        ! arguments out
          REAL score(3)

      ! calculate response for scoring function
        DO j=1,n1
          rsp(j)=resp(j) 
        END DO

      ! store and evaluate variables
        CALL copytree(ntr,nkn,conc,negs,pick,term,-1,3,1)
        CALL storing(nkn,ntr,conc,pick,npckmv,pickmv,ssize,nop)
        DO wh=1,ntr
          CALL evaluate_first(wh,n1,n2,nkn,ntr,conc,term,negs,pick,
     #                        datri,prtr)
        END DO

      ! choose scoring function for model type
        IF (mdl.EQ.1) THEN
          score(1)=0.0
          DO i=1,n1
            score(1)=score(1)+
     #               weight(i)*((REAL(prtr(i,1))-resp(i))**2.0)
          END DO
        ELSE IF (mdl.EQ.2) THEN
          DO j=1,3
            score(j)=0        
          END DO
          CALL calcrss(nop,n1,ntr,betas,prtr,nsep,seps,
     #                 rsp,weight,score)
        ELSE IF (mdl.EQ.3) THEN
          CALL scoredev(n1,nop,ntr,prtr,nsep,seps,rsp,
     #                  weight,betas,score)
        ELSE IF (mdl.EQ.4) THEN
          CALL scorepll(n1,nop,ntr,nsep,seps,prtr,betas,dcph,ordrs,
     #                  score,weight)
        ELSE IF (mdl.EQ.5) THEN 
           CALL exposcore(prtr,rsp,dcph,weight,n1,ntr,nop,
     #                         nsep,seps,score(1),betas)
        ELSE IF (mdl.EQ.0) THEN 
           CALL My_own_scoring(prtr,rsp,dcph,ordrs,weight,n1,ntr,nop,
     #                         wh,nsep,seps,score(1),betas)
        ELSE 
          CALL stringprint('not done yet!',13)
          STOP
        END IF

      END 

      ! *****************************************************************
      ! *****************************************************************



      ! this subroutine calculates the deviance for logistic regression
      ! last modification 10/17/02

      SUBROUTINE scoredev(n1,nop,ntr,prtr,nsep,seps,rsp,weight,betas,
     #                    score)
      IMPLICIT NONE

        ! arguments in
          INTEGER n1,nop,nsep,ntr,prtr(n1,ntr)
          REAL rsp(n1),weight(n1),betas(0:(nsep+ntr))
          REAL seps(nsep,n1)
        ! local
          INTEGER i,j
          DOUBLE PRECISION prd,myexp,mylog
        ! arguments out
          REAL score(3)

      ! calculate the predictions
        score(1)=0.0
        DO i=1,n1
          prd=betas(0)
          IF (nsep.GT.0) THEN
            DO j=1,nsep
              prd=prd+betas(j)*seps(j,i)
            END DO
          END IF
          IF (nop.GT.0) THEN
            DO j=1,nop
              prd=prd+betas(nsep+j)*REAL(prtr(i,j))
            END DO
          END IF
          prd=myexp(prd)
          prd=prd/(1.0+prd)

      ! calculate the deviance
          IF (prd.GT.0.D0.AND.prd.LT.1.D0) THEN
            IF (rsp(i).EQ.0) prd=1-prd
            score(1)=score(1)-2.0*weight(i)*mylog(prd)
          ELSE
            STOP " * Fitted probabilities of 0 or 1 * "
          ENDIF
        ENDDO

      END 

      ! *****************************************************************
      ! *****************************************************************

      ! this subroutine calculates the partial likelihood for the 
      ! proportional hazards model
      ! last modification 10/17/02

      SUBROUTINE scorepll(n1,nop,ntr,nsep,seps,prtr,betas,dcph,ordrs,
     #                    score,weight)
      IMPLICIT NONE

        ! parameters 
          INTEGER LGCbetaMAX,LGCn1MAX
          PARAMETER (LGCn1MAX   = 10000)
          PARAMETER (LGCbetaMAX=    55)
        ! arguments in
          INTEGER n1,nop,ntr,nsep,dcph(n1),ordrs(n1),prtr(n1,ntr)
          REAL betas(0:(nsep+ntr)),seps(nsep,n1),weight(n1)
        ! local
          INTEGER j,k,nnf(2),xstop
          DOUBLE PRECISION loglf,betaf(LGCbetaMAX)
          DOUBLE PRECISION covsf(LGCn1MAX*LGCbetaMAX)
        ! arguments out
          REAL score(3),r

        xstop=0
        CALL stopper(LGCn1MAX,n1,"LGCn1MAX","scorepll()",8,10,xstop,0)
        CALL stopper(LGCbetaMAX,nsep+ntr,"LGCbetaMAX","scorepll()",10,
     #          10,xstop,1)

      ! calculate the partial likleihood
        nnf(1)=nop+nsep
        nnf(2)=n1
        DO j=1,(n1*(nsep+ntr))
          covsf(j)=0.D0
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
        DO j=1,(nsep+nop)
          betaf(j)=betas(j)
        END DO
        j=nsep+nop
        CALL mypllxx(loglf,betaf,dcph,ordrs,covsf,j,n1,weight)
        r=n1
        score(1)=-loglf

      END 

      ! *****************************************************************
      ! *****************************************************************

      ! This routines does the optimization for a Cox PH model
      ! Last modified 10/17/02

      SUBROUTINE myphxx(delta,idx,covs,np,n1,nsep,ntr,logl,beta,
     #                   oops,weight)
      IMPLICIT none
        ! parameters 
          INTEGER LGCbetaMAX
          PARAMETER (LGCbetaMAX=    55)
        ! i/o
          INTEGER n1,nsep,ntr,zolala,delta(n1),idx(n1),oops
          DOUBLE PRECISION beta(nsep+ntr),covs(n1*(nsep+ntr))
          DOUBLE PRECISION logl
          REAL weight(n1)
        ! local
          DOUBLE PRECISION ologl,nlogl,alpha,pp,alphap1,alphap2,prec
          DOUBLE PRECISION nbeta(LGCbetaMAX),grad(LGCbetaMAX)
          DOUBLE PRECISION hess(LGCbetaMAX,LGCbetaMAX)   
          INTEGER iter,i,np,xstop

        xstop=0
        CALL stopper(LGCbetaMAX,nsep+ntr,"LGCbetaMAX","myphxx()",10,8,
     #          xstop,1)

          DO i=1,np
            beta(i)=0.
          END DO
          oops=0
          alphap1=0.001
          alphap2=0.00001
          alpha=1
          prec=0.00001
          iter=0
          pp=prec+10
          DO WHILE((iter.LT.10).AND.(pp.GT.prec).AND.(alpha.GT.alphap2))
            iter=iter+1
            CALL mygradph(grad,hess,beta,delta,idx,covs,np,n1,ologl,
     #                    LGCbetaMAX,weight)
            DO i=1,np
              IF((hess(i,i).LT.1.0e-10).AND.(hess(i,i).GT.-1.0e-10))THEN
                CALL mypllxx(logl,beta,delta,idx,covs,np,n1,weight)
                GOTO 1234
              END IF
            END DO
            CALL lusolveph(hess,grad,np,oops,LGCbetaMAX)
            IF(oops.EQ.1)RETURN
            alpha=1
            zolala = 0
            DO WHILE (((alpha.GT.alphap2) .AND. (nlogl.LT.ologl))
     #                               .OR.(zolala.EQ.0))
              zolala = 1
              DO i=1,np
                nbeta(i)=beta(i)+alpha*grad(i)
              END DO
              CALL mypllxx(nlogl,nbeta,delta,idx,covs,np,n1,weight)
              IF(nlogl.LT.ologl) alpha=alpha/2
            END DO
            IF(alpha.GT.alphap1) THEN
              pp=0
              DO i=1,np
                pp=pp+(nbeta(i)-beta(i))*(nbeta(i)-beta(i))
                beta(i)=nbeta(i)
              END DO
              pp=SQRT(pp)
              IF(iter.LT.3) pp=prec+10.
            END IF
          END DO
          CALL mygradph(grad,hess,beta,delta,idx,covs,np,n1,logl,
     #                  LGCbetaMAX,weight)
1234      CONTINUE

      END 

      ! *****************************************************************
      ! *****************************************************************

      ! This routines computes the gradient for a Cox PH model
      ! Last modified 10/17/02

      SUBROUTINE mygradph(grad,hess,beta,delta,idx,covs,np,n1,logl,lda,
     #                    weight)
      IMPLICIT none
        ! parameters 
          INTEGER LGCbetaMAX,LGCn1MAX
          PARAMETER (LGCn1MAX=   10000)
          PARAMETER (LGCbetaMAX=    55)
        ! i/o
          INTEGER n1,np,delta(n1),idx(n1),lda
          DOUBLE PRECISION beta(np),grad(np),covs(n1*np),hess(lda,np)
          DOUBLE PRECISION logl
          REAL weight(n1)
        ! local
          DOUBLE PRECISION ff(LGCn1MAX),s1(LGCn1MAX),s1s(LGCn1MAX)
          DOUBLE PRECISION gg(LGCn1MAX),ff2(LGCn1MAX),myexp,mylog
          DOUBLE PRECISION s0,s1r,u,z,s2(LGCbetaMAX*LGCbetaMAX)
          INTEGER i,i2,j,k,r,s,it,xstop
          xstop=0
          CALL stopper(LGCn1MAX,n1,"LGCn1MAX","mygradph()",8,10,xstop,0)
          CALL stopper(LGCbetaMAX,np,"LGCbetaMAX","mygradph()",10,10,
     #                 xstop,1)

          u=0
          DO i=1,n1
            ff(i)=0
            DO k=1,np
              ff(i)=ff(i)+beta(k)*covs(i+n1*(k-1))
            END DO
          END DO
          s0=0
          DO i=1,np
            grad(i)=0
            s1(i)=0
            DO k=1,np
              s2((i-1)*np+k)=0
              hess(i,k)=0
            END DO
          END DO
          DO i=1,n1
             gg(i)=ff(idx(i))
             ff2(i)=myexp(gg(i))
          END DO
          DO i2=1,n1
            i=n1+1-i2
            j=idx(i)
            s0=s0+ff2(i)*weight(j)
            DO r=1,np
              s1r=ff2(i)*covs(j+n1*(r-1))*weight(j)
              s1(r)=s1(r)+s1r
              it=(r-1)*np
              DO s=r,np
                s2(it+s)=s2(it+s)+s1r*covs(j+n1*(s-1))
              END DO
            END DO
            IF(delta(j).EQ.1) THEN
              DO r=1,np
                 s1s(r)=s1(r)/s0
              END DO
              DO r=1,np
                it=(r-1)*np
                grad(r)=grad(r)+
     #                  weight(j)*(covs(j+n1*(r-1))-s1s(r))
                DO s=r,np
                   hess(r,s)=hess(r,s)-
     #                  weight(j)*(s1s(r)*s1s(s)-s2(it+s)/s0)
                END DO
              END DO
              z=ff2(i)/s0
              z=mylog(z)
              u=u+z*weight(j)
            END IF
          END DO
          DO r=1,np
            DO s=1,r
              hess(r,s)=hess(s,r)
            END DO
          END DO
          logl=u

      END 

      ! *****************************************************************
      ! *****************************************************************

      ! This routine computes the partial likelihood
      ! Last modified 10/17/02

      SUBROUTINE mypllxx(logl,beta,delta,idx,covs,np,n1,weight)
      IMPLICIT none
        ! i/o
          INTEGER n1,np,delta(n1),idx(n1),LGCn1MAX
          PARAMETER (LGCn1MAX=10000)
          DOUBLE PRECISION beta(np),covs(n1*np),logl,myexp,mylog
          REAL weight(n1)
        ! local
          INTEGER i,k,xstop
          DOUBLE PRECISION z,s0,ff(LGCn1MAX),ff2(LGCn1MAX),gg(LGCn1MAX)
          xstop=0
          CALL stopper(LGCn1MAX,n1,"LGCn1MAX","mypllxx()",8,9,xstop,1)

          logl=0.
          DO i=1,n1
            ff(i)=0.
            DO k=1,np
              ff(i)=ff(i)+beta(k)*covs(i+n1*(k-1))
            END DO
          END DO
          DO i=1,n1
             gg(i)=ff(idx(i))
             ff2(i)=myexp(gg(i))
          END DO
          s0=0.
          DO k=1,n1
            i=n1+1-k
            s0=s0+ff2(i)*weight(idx(i))
            IF(delta(idx(i)).EQ.1) THEN
              z=ff2(i)/s0
              z=mylog(z)
              logl=logl+z*weight(idx(i))
            END IF
          END DO

      END 

      ! *****************************************************************
      ! *****************************************************************

      ! This routine solves a system
      ! Last modified 10/17/02

      SUBROUTINE lusolveph(a,b,n,k,lda)
      IMPLICIT none

          ! parameters
          INTEGER lda,LGCbetaMAX
          PARAMETER (LGCbetaMAX   = 55)
          INTEGER n,k
          DOUBLE PRECISION a(lda,n)
          DOUBLE PRECISION b(n)
          INTEGER job,info
          INTEGER ipvt(LGCbetaMAX),xstop
          xstop=0
          CALL stopper(LGCbetaMAX,n,"LGCbetaMAX","lusolveph()",10,11,
     #                 xstop,1)

          job=0
          k=0
          if(n.GT.0)THEN
             CALL dgefa(a,lda,n,ipvt,info)
             IF(info.NE.0)THEN
               k=1
             ELSE
               CALL dgesl(a,lda,n,ipvt,b,job)
             END IF
          END IF

      END 

      ! *****************************************************************
      ! *****************************************************************

      ! This routine glues the directory to the filename
      ! Last modified 10/17/02

      SUBROUTINE STRINGCOM(s1,s2,ls1,ls2)
      IMPLICIT none

          CHARACTER *80 s1,s2,s3
          INTEGER ls1,ls2

          s3(1:ls1)=s1(1:ls1)
          s3((ls1+1):(ls1+ls2))=s2(1:ls2)
          s2=s3
          ls2=ls2+ls1

      END 
      
      ! *****************************************************************
      ! *****************************************************************



      ! This function calculates the sum of an array (integers)
      ! Last modified 10/17/02

      INTEGER FUNCTION SUM2I(x,n1,n2,iin,f1,fr,tto)

        ! arguments in
          INTEGER n1,n2,iin,f1,fr,tto
          INTEGER x(n1,n2)
        ! local
          INTEGER j

      ! here it goes
        SUM2I=0   

        IF (iin.GT.2) STOP "index in SUM2I bigger than 2"
        IF (iin.EQ.1) THEN
          DO j=fr,tto
            SUM2I=SUM2I+x(j,f1)
          ENDDO
        ELSE IF (iin.EQ.2) THEN
          DO j=fr,tto
            SUM2I=SUM2I+x(f1,j)
          ENDDO
        ENDIF

      END

      ! *****************************************************************
      ! *****************************************************************

      SUBROUTINE clearly(iearly,ntr,nkn,n2)
      
      ! parameters
          INTEGER ntr,nkn,n2
          INTEGER iearly(0:6,ntr,0:nkn,n2,0:1,2)
          INTEGER i,j,k,l
          DO i=1,ntr
            DO j=0,nkn
              DO k=1,n2
                DO l=0,6
                  iearly(l,i,j,k,0,1)=0
                  iearly(l,i,j,k,1,1)=0
                  iearly(l,i,j,k,0,2)=0
                  iearly(l,i,j,k,1,2)=0
                END DO
              END DO
            END DO
          END DO
          END 
      ! *****************************************************************
      ! *****************************************************************
      SUBROUTINE  stopper(i,j,s1,s2,m,n,vv,ww)

        INTEGER i,j,m,n,vv,ww
        CHARACTER*20 s1,s2
        CHARACTER*80 s3

      IF(i.LT.j)THEN
        CALL stringprint("Insufficient declaration",24)
             
        s3(1:m)=s1(1:m)
        s3((m+1):(m+4))=" in "
        s3((m+5):(m+4+n))=s2(1:n)
        s3((m+5+n):(m+8+n))=" is "
        CALL makeistring((m+9+n),(m+16+n),s3,i,8)
        CALL stringprint(s3,m+16+n)
        s3(1:m)=s1(1:m)
        s3((m+1):(m+21))=" should be at least "
        CALL makeistring((m+22),(m+29),s3,j,8)
        CALL stringprint(s3,m+29)
        vv=vv+1
      END IF
      IF(ww.GT.0.AND.vv.GT.0)THEN
        CALL stringprint('Please fix and recompile....',28)
        STOP
      END IF
      END
      ! *****************************************************************
      ! *****************************************************************
      SUBROUTINE makeistring(k1,k2,astring,i,j)
      INTEGER i,j,k,k1,k2
      CHARACTER *40 aa
      CHARACTER *80 astring
      CALL makeiistring(aa,i,j,k,0)
      astring(k1:k2)=aa(1:(k2-k1+1))
      END
      SUBROUTINE makerstring(k1,k2,astring,rr,i,j)
      CHARACTER *20 aa,bb
      CHARACTER *80 astring
      REAL r,rr
      INTEGER i,j,p,k,l,s2,k3,ww
      ww=1
      IF(rr.LT.0.AND.rr.GT.(-9))ww=0
      r=rr
      s2=1
      IF(r.LT.0)THEN
         s= -r
      ELSE
         s=r
      END IF
      k=s
      if(r.LT.0)k=-k
      k3=k
      IF(j.GT.0)THEN
        k=s
        s=s-k
        DO k=1,j
           s2=s2*10
           s=s*10
        END DO
        k=(s+0.5)
      END IF
      IF(k.eq.s2)THEN
         k=0
         IF(r.GE.0)THEN
            k3=k3+1
         ELSE
            k3=k3-1
         END IF
      END IF
      CALL makeiistring(bb,k3,i,p,0)
      aa(1:i)=bb(1:i)
      aa((i+1):(i+1))="."
      IF(j.GT.0)THEN
        CALL makeiistring(bb,k,j,l,1)
        aa((i+2):(i+1+j))=bb(1:j)
      END IF
      IF(r.lt.0.and.k2.eq.0.and.i.gt.1)aa((i-1):(i-1))="-"
      IF(p.eq.1)THEN
        DO k=1,(i+j+1)
          aa((i+j):(i+j))="*"
        END DO
      END IF
      IF(ww.EQ.0.and.i.gt.1)aa((i-1):(i-1))="-"
      astring(k1:k2)=aa(1:(k2-k1+1))
      END
      
      SUBROUTINE makeiistring(aa,i,j,p,f)
      CHARACTER *20 aa
      INTEGER i,j,p,f
      INTEGER i2,k,l,m,n
      i2=i
      p=0
      DO k=1,20
         aa(k:k)=' '
      END DO
      IF(i.EQ.0) THEN
        IF(f.EQ.0) THEN
          aa(j:j)='0'
        ELSE
          DO k=1,20
            aa(k:k)='0'
          END DO
        END If
      ELSE
        IF(i.GE.0) THEN
          l=0
        ELSE
          l=1
          i2=-i2
        END IF
        m=0
        DO k=1,j
          i3=i2
          m=i2/10
          n=i2-10*m
          i2=m
          IF(i3.GT.0)THEN
            aa((j+1-k):(j+1-k))=char(48+n)
            IF(i2.EQ.0.AND.l.eq.1)THEN
              IF(k.NE.J)THEN
                aa((j-k):(j-k))='-'
              ELSE
                DO i4=1,j
                  aa(i4:i4)='*'
                END DO
                p=1
              END IF
            END IF
          ELSE
            IF(f.EQ.1)THEN
              aa((j+1-k):(j+1-k))="0"
            END IF
          END IF
        END DO
        IF(i2.ne.0)then
          DO i4=1,j
            aa(i4:i4)='*'
          END DO
          p=1
        END IF
      END IF
      END
      SUBROUTINE stringprint(aa,i)
      CHARACTER *80 aa
      INTEGER i
      REAL k
      CALL realpr(aa,i,k,0)
      END
      SUBROUTINE reorder(iotree,conc,term,negs,pick,ltree,k,j,nkn,
     #                   ntrnew)
      INTEGER iotree(1),ltree,k,j,nkn,ntrnew
      INTEGER conc(nkn,ntrnew,3)
      INTEGER negs(nkn,ntrnew,3)
      INTEGER pick(nkn,ntrnew,3)
      INTEGER term(nkn,ntrnew,3)
                    iotree((ltree-1)*(4*nkn+3)+(k-1)*4+4)=conc(k,j,3)
                    iotree((ltree-1)*(4*nkn+3)+(k-1)*4+5)=term(k,j,3)
                    iotree((ltree-1)*(4*nkn+3)+(k-1)*4+6)=negs(k,j,3)
                    iotree((ltree-1)*(4*nkn+3)+(k-1)*4+7)=pick(k,j,3)
      END
      SUBROUTINE smackonprior(score,ssize,ntr,nkn,conc,term,negs,pick,
     #           hyperpars,n2,mtp,slprbc,rrat,nopdiff)
      IMPLICIT none
      INTEGER ssize,ntr,nkn,n2,term(nkn,ntr,3),mtp,nopdiff
      REAL ll,hyperpars(10),score(3),slprbc(6),rrat
      INTEGER conc(nkn,ntr,3),negs(nkn,ntr,3),pick(nkn,ntr,3)
      DOUBLE PRECISION mylog,postrat,rr,zz
      CHARACTER *100 astring
      DOUBLE PRECISION u1,u2,u3
      CALL getv(rr,ssize,ntr,nkn,n2)
      ! a simple hypergeometric prior
      score(1)=0.5*score(1)*exp(hyperpars(2))+hyperpars(1)*ssize
      hyperpars(10)=score(1)
      score(1)=score(1)+rr
      hyperpars(9)=score(1)
      hyperpars(8)=0.
      zz=(slprbc(3)-slprbc(2))/(slprbc(4)-slprbc(2)+slprbc(1))
      IF(mtp.GE.0)THEN
         postrat=1.
         IF(mtp.EQ.3)postrat=rrat/(8*n2)
         IF(mtp.EQ.4)postrat=8*n2*rrat
         IF(mtp.EQ.5)postrat=8*n2*rrat
         IF(mtp.EQ.6)postrat=rrat/(8*n2)
         IF(mtp.EQ.3.AND.ssize.EQ.0)postrat=1./(zz*(2*n2))
         IF(mtp.EQ.3.AND.nopdiff.GT.0)postrat=1./(zz*(2*n2))
         IF(mtp.EQ.0)postrat=(2*n2)*zz
         hyperpars(8)=mylog(postrat)
      END IF
      END 
      SUBROUTINE getv2(nn,ssize,nkn,n2)
      DOUBLE PRECISION nn
      INTEGER i,nkn,n2,ssize
      DOUBLE PRECISION r1,r2,myexp,mylog
      nn=0.
      DO i=0,ssize
         CALL getv1(r1,i,nkn,n2)
         CALL getv1(r2,ssize-i,nkn,n2)
         nn=nn+myexp(r1+r2)
      END DO
      nn=mylog(nn)
      END
      SUBROUTINE getv3(nn,ssize,nkn,n2)
      DOUBLE PRECISION nn
      INTEGER i,nkn,n2,ssize
      DOUBLE PRECISION r1,r2,myexp,mylog
      nn=0.
      DO i=0,ssize
         CALL getv1(r2,i,nkn,n2)
         CALL getv2(r1,ssize-i,nkn,n2)
         nn=nn+myexp(r1+r2)
      END DO
      nn=mylog(nn)
      END
      SUBROUTINE getv4(nn,ssize,nkn,n2)
      DOUBLE PRECISION nn
      INTEGER i,nkn,n2,ssize
      DOUBLE PRECISION r1,r2,myexp,mylog
      nn=0.
      DO i=0,ssize
         CALL getv1(r1,i,nkn,n2)
         CALL getv3(r2,ssize-i,nkn,n2)
         nn=nn+myexp(r1+r2)
      END DO
      nn=mylog(nn)
      END
      SUBROUTINE getv5(nn,ssize,nkn,n2)
      DOUBLE PRECISION nn
      INTEGER i,nkn,n2,ssize
      DOUBLE PRECISION r1,r2,myexp,mylog
      DO i=0,ssize
         CALL getv1(r1,i,nkn,n2)
         CALL getv4(r2,ssize-i,nkn,n2)
         nn=nn+myexp(r1+r2)
      END DO
      nn=mylog(nn)
      END
      SUBROUTINE  getv(ll,ssize,ntr,nkn,n2) 
      INTEGER ssize,ntr,nkn,n2
      DOUBLE PRECISION ll
      IF(ntr.EQ.1)THEN
         CALL getv1(ll,ssize,nkn,n2)
         goto 567
      END IF
      IF(ntr.EQ.2)THEN
         CALL getv2(ll,ssize,nkn,n2)
         goto 567
      END IF
      IF(ntr.EQ.3)THEN
         CALL getv3(ll,ssize,nkn,n2)
         GOTO 567
      END IF
      IF(ntr.EQ.4)THEN
         CALL getv4(ll,ssize,nkn,n2)
         GOTO 567
      END IF
      IF(ntr.EQ.3)THEN
         CALL getv5(ll,ssize,nkn,n2)
         GOTO 567
      END IF
567   CONTINUE
      END 
      SUBROUTINE getv1(l,ssize,nb,n2)
      INTEGER ll,ssize,nb,precomp(25),n2
      DOUBLE PRECISION l,m,mylog,r
      IF(ssize.GT.20)THEN
        ll=0
        GOTO 568
      END IF
      IF(ssize.LE.2)THEN
        ll=1
        GOTO 568
      END IF
      IF((ssize*2-1).GT.nb)THEN
        ll=0
        GOTO 568
      END IF
      precomp(3)=2
      IF(nb.LT.8)THEN
        precomp(4)=1
        ll=precomp(ssize)
        GOTO 568
      END IF
      precomp(4)=5
      IF(nb.LT.16)THEN
        precomp(5)=6
        precomp(6)=6
        precomp(7)=4
        precomp(8)=1
        ll=precomp(ssize)
        GOTO 568
      END IF
      precomp(5)=14
      IF(nb.LT.32)THEN
        precomp(6)=26
        precomp(7)=44
        precomp(8)=69
        precomp(9)=94
        precomp(10)=114
        precomp(11)=116
        precomp(12)=94
        precomp(13)=60
        precomp(14)=28
        precomp(15)=8
        precomp(16)=1
        ll=precomp(ssize)
        GOTO 568
      END IF
      precomp(6)=42
      IF(nb.LT.64)THEN
        precomp(7)=100
        precomp(8)=221
        precomp(9)=470
        precomp(10)=958
        precomp(11)=1860
        precomp(12)=3434
        precomp(13)=6036
        precomp(14)=10068
        precomp(15)=15864
        precomp(16)=23461
        precomp(17)=32398
        precomp(18)=41658
        precomp(19)=49700
        precomp(20)=54746
        ll=precomp(ssize)
        goto 568
      END IF
      precomp(7)=132
      IF(nb.LT.128)THEN
        precomp(8)=365
        precomp(9)=950
        precomp(10)=2398
        precomp(11)=5916
        precomp(12)=14290
        precomp(13)=33708
        precomp(14)=77684
        precomp(15)=175048
        precomp(16)=385741
        precomp(17)=831014
        precomp(18)=1749654
        precomp(19)=3598964
        precomp(20)=7228014
        ll=precomp(ssize)
        goto 568
      END IF
      precomp(8)=429
      IF(nb.LT.256)THEN
        precomp(9)=1302
        precomp(10)=3774
        precomp(11)=10652
        precomp(12)=29538
        precomp(13)=80812
        precomp(14)=218324
        precomp(15)=582408
        precomp(16)=1534301
        precomp(17)=3993030
        precomp(18)=10269590
        precomp(19)=26108844
        precomp(20)=65626918
        ll=precomp(ssize)
        GOTO 568
      END IF
      precomp(9)=1430
      precomp(10)=4862
      precomp(11)=16796
      precomp(12)=58786
      precomp(13)=208012
      precomp(14)=742900
      precomp(15)=2674440
      precomp(16)=9694845
      precomp(17)=35357670
      precomp(18)=129644790
      precomp(19)=477638700
      precomp(20)=1767263190
      ll=precomp(ssize)
568   CONTINUE
      l=ll
      l=mylog(l)
      m=2.
      m=mylog(m)
      r=n2
      IF(ssize.GT.0)THEN
        l=l+m*(2*ssize-1)+ssize*mylog(r)
      END IF
      END 
      SUBROUTINE storeone(mcmc,new,hyperpars,lvisit,visit,
     #           ntr,nkn,conc,negs,pick,term,nac,rd1,rd2,rd3,rd4,
     #           bout,n2)
      IMPLICIT NONE
      INTEGER LGCn2MAX,LGCntrMAX
      PARAMETER(LGCn2MAX=1000,LGCntrMAX=5)
      INTEGER mcmc,new,nkn,ntr,i,j,k,ssize,n2,used1(LGCn2MAX)
      INTEGER conc(nkn,ntr,3),negs(nkn,ntr,3),pick(nkn,ntr,3)
      INTEGER visit(2+ntr*nkn),term(nkn,ntr,3),nac,rd1(1)
      REAL lvisit(2),rd2(1),rd3(1)
      REAL hyperpars(10)
      INTEGER iz2,iz,jz,i2,j2,xstop,rd4(1),zused
      INTEGER iz3,j3,bout,xused(LGCntrMAX,LGCn2MAX),yused
      INTEGER vused(LGCntrMAX)
      xstop=0
      CALL stopper(LGCn2MAX,n2,"LGCn2MAX","storeone",9,8,xstop,1)
      xstop=0
      CALL stopper(LGCntrMAX,ntr,"LGCntrMAX","storeone",9,9,xstop,1)
      IF(new.eq.nac)THEN
         new=0
      ELSE
         new=1
      END IF
      IF(nac.LT.-10) new=2
      IF(new.GT.0)THEN
        IF(new.NE.2)THEN
          DO i=1,n2
            used1(i)=0
          END DO
          visit(2)=visit(2)+mcmc-2
          DO i2=1,n2
            used1(i2)=0
          END DO
          jz=0
          DO i=1,ntr
            DO j=1,nkn
              IF(conc(j,i,3).EQ.3)THEN
                jz=jz+1
              END IF
            END DO
            DO i2=1,n2
              xused(i,i2)=0
            END DO
            DO i2=1,nkn
              IF(conc(i2,i,3).EQ.3)THEN
                xused(i,term(i2,i,3))=1
                used1(term(i2,i,3))=1
              END IF
            END DO
          END DO
          DO iz=1,n2
            IF(used1(iz).EQ.1)THEN
              rd2(iz)=rd2(iz)+visit(2)
            END IF
          END DO
          IF(bout.GE.2.or.bout.LE.(-2))THEN
            DO iz=1,(n2-1)
              IF(used1(iz).EQ.1)THEN
                DO iz2=(iz+1),n2
                  IF(used1(iz2).EQ.1)THEN
                    yused=0
                    DO i=1,ntr
                      vused(i)=xused(i,iz)*xused(i,iz2)
                      yused=yused+vused(i)
                    END DO
                    IF(yused.GT.0)THEN
                      rd3(iz2+(iz-1)*n2)=rd3(iz2+(iz-1)*n2)+visit(2)
                    END IF
                    IF(bout.GE.3.or.bout.LE.(-3).and.yused.GT.0)THEN
                      DO iz3=(iz2+1),n2
                        IF(used1(iz3).EQ.1)THEN
                          zused=0
                          DO i=1,ntr
                            zused=zused+vused(i)*xused(i,iz3)
                          END DO
                          IF(zused.GT.0)THEN
                            i2=iz3+(iz2-1)*n2+(iz-1)*n2*n2
                            rd4(i2)=rd4(i2)+visit(2)
                          END IF
                        END IF
                      END DO
                    END IF
                  END IF
                END DO
              END IF
            END DO
          END IF
          rd1(jz+1)=rd1(jz+1)+visit(2)
        END IF
        IF(new.EQ.2)THEN
          DO i=1,128
            rd1(i)=0
          END DO
          DO i=1,n2
            rd2(i)=0.
            DO j=1,n2
              rd3(i+(j-1)*n2)=0.
            END DO
          END DO
          lvisit(1)=0
          lvisit(2)=0
        END IF
        IF(bout.GT.0)THEN
          DO i=1,ntr
            DO j=1,nkn
              IF(conc(j,i,3).LT.3)THEN
                visit(2+(i-1)*nkn+j)=1000*conc(j,i,3)
              ELSE
                IF(negs(j,i,3).EQ.0)THEN
                  visit(2+(i-1)*nkn+j)=term(j,i,3)
                ELSE
                  visit(2+(i-1)*nkn+j)= -term(j,i,3)
                END IF
              END IF
            END DO
          END DO
          k=2+ntr*nkn
          CALL cwrite(lvisit,visit,k)
        END IF
        lvisit(1)=hyperpars(9)
        lvisit(2)=hyperpars(10)
        visit(1)=1
        visit(2)=1
      ELSE
        visit(2)=visit(2)+mcmc-2
        visit(1)=visit(1)+1
        visit(2)=visit(2)+1
      END IF
      END 
      SUBROUTINE redater(ii,x,y,w,prtr,seps,resp,wt,n1,ntr,nsep,ln)
      INTEGER n1,ntr,nsep,ii,prtr(n1,ntr)
      DOUBLE PRECISION y(ln),w(ln),x(ln,nsep+n1+1)
      REAL wt(n1),seps(nsep,n1),resp(n1)
      INTEGER i,j,k,l,i2,j2
      REAL rt
      ii=1
      DO i=1,nsep
        DO j=1,n1
          IF(seps(i,j).GT.1.000001.OR.seps(i,j).LT.-0.000001)THEN
            ii=0
            GOTO 117
          ELSE
            IF(seps(i,j).GT.0.000001.AND.seps(i,j).LE.0.999999)THEN
              ii=0
              GOTO 117
            END IF
          END IF
        END DO
      END DO
      rt=nsep+ntr
      i2=nsep+ntr
      k=2.0**rt
      DO i=1,k
         DO j=1,i2
            x(i,j+1)=0.
         END DO
         y(i)=0.
         w(i)=0.
         x(i,1)=1.
      END DO
      DO i=1,k
         j=k/2
         l=i-1
         DO j2=1,i2
            IF(l.GE.j)THEN
               x(i,i2+2-j2)=1.
               l=l-j
            END IF
            j=j/2
         END DO
      END DO
      DO i=1,n1
         j=0
         k=1
         DO l=1,nsep
            j=j+k*seps(l,i)
            k=k*2
         END DO
         DO l=1,ntr
            j=j+k*prtr(i,l)
            k=k*2
         END DO
         IF(resp(i).EQ.1)y(j+1)=y(j+1)+wt(i)
         w(j+1)=w(j+1)+wt(i)
      END DO
      ii=0
      DO i=1,k
         IF(w(i).GT.0)THEN
            ii=ii+1
            w(ii)=w(i)
            y(ii)=y(i)
            DO j=1,(1+nsep+ntr)
               x(ii,j)=x(i,j)
            END DO
         END IF
      END DO
117   CONTINUE 
      END
      SUBROUTINE copytree(ntr,nkn,conc,negs,pick,term,one,iin,iout)
      IMPLICIT NONE 
      INTEGER ntr,nkn,one,i1,i2,i3,i4,iin,iout
      INTEGER conc(nkn,ntr,3)
      INTEGER negs(nkn,ntr,3)
      INTEGER pick(nkn,ntr,3)
      INTEGER term(nkn,ntr,3)
      i3=one
      i4=one
      IF(one.LT.0)THEN
         i3=1
         i4=ntr
      END IF
      DO i1=1,nkn
        DO i2=1,i4
          conc(i1,i2,iout)=conc(i1,i2,iin)
          negs(i1,i2,iout)=negs(i1,i2,iin)
          pick(i1,i2,iout)=pick(i1,i2,iin)
          term(i1,i2,iout)=term(i1,i2,iin)
        END DO
      END DO
      END

      SUBROUTINE greedyonestep(nkn,ntr,conc,negs,pick,term,n1,n2,
     #           storage,prtr,datri,nsep,seps,cbetas,xtxsep,mdl,
     #           nop,orders,dcph,rsp,weight,bestscore,mtm,cnc,
     #           penalty)
      IMPLICIT NONE
      INTEGER nkn,ntr,n1,n2,nsep,mdl,nop,mtm,cnc(3)
      INTEGER conc(nkn,ntr,3)
      INTEGER negs(nkn,ntr,3)
      INTEGER pick(nkn,ntr,3)
      INTEGER term(nkn,ntr,3)
      INTEGER storage(2*ntr*nkn*n1),orders(n1),dcph(n1)
      INTEGER prtr(n1,ntr),datri(n2,n1)
      REAL seps(nsep,n1),cbetas(0:(nsep+ntr)),xtxsep(0:nsep,0:nsep)
      REAL rsp(n1),weight(n1),bestscore,penalty

      INTEGER mtp,wh,knt,letter,op,neg,sng,dbl,improve,bestmove(3)
      INTEGER yes,l2
      improve=-1
      CALL copytree(ntr,nkn,conc,negs,pick,term,-1,1,2)
      CALL copytree(ntr,nkn,conc,negs,pick,term,-1,1,3)
      
      DO wh=1,nop
        DO mtp=1,6
          DO knt=1,nkn
            CALL isallowed(wh,mtp,knt,negs,pick,term,conc,nkn,ntr,yes)
            IF(yes.GE.1)THEN
              CALL copytree(ntr,nkn,conc,negs,pick,term,-1,2,1)
              letter=1
              op=1
              neg=1
              IF(mtp.EQ.1)THEN
                DO letter=1,n2
                  DO neg=0,1
                    IF(letter.NE.term(knt,wh,1).OR.
     #                     neg.NE.negs(knt,wh,1))THEN
                      l2=letter
                      CALL altlf(knt,n2,nkn,ntr,wh,conc,negs,term,l2,
     #                                        neg)
                      IF(l2.LT.0)GOTO 391
                    
                      CALL evalgreed(nkn,ntr,conc,pick,negs,term,mtp,n1,
     #                     n2,wh,knt,storage,datri,prtr,mdl,nsep,dcph,
     #                     orders,mtm,rsp,weight,seps,cbetas,xtxsep,nop,
     #                     bestscore,improve,0,bestmove,penalty)
                    END IF
                    CALL copytree(ntr,nkn,conc,negs,pick,term,-1,2,1)
                  END DO
391               CONTINUE
                END DO
              END IF
              IF(mtp.EQ.4.OR.mtp.EQ.5)THEN
                DO letter=1,n2
                  DO neg=0,1
                    DO op=1,2
                      l2=letter
                      IF(mtp.EQ.4)CALL xsplit(knt,n2,nkn,ntr,wh,cnc,
     #                            conc,negs,pick,term,l2,op,neg)
                      IF(mtp.EQ.5)CALL branch(knt,n2,nkn,ntr,wh,cnc,
     #                            conc,negs,pick,term,l2,op,neg)
                      IF(l2.LT.0)GOTO 394
                      CALL evalgreed(nkn,ntr,conc,pick,negs,term,mtp,n1,
     #                     n2,wh,knt,storage,datri,prtr,mdl,nsep,dcph,
     #                     orders,mtm,rsp,weight,seps,cbetas,xtxsep,nop,
     #                     bestscore,improve,0,bestmove,penalty)
                      CALL copytree(ntr,nkn,conc,negs,pick,term,-1,2,1)
                    END DO
                  END DO
394               CONTINUE
                END DO
              END IF
              IF(mtp.EQ.2.OR.mtp.EQ.3.OR.mtp.EQ.6)THEN
                IF(mtp.EQ.2)conc(knt,wh,1)=3-conc(knt,wh,1)
                IF(mtp.EQ.3)
     #              CALL xdelete(knt,nkn,ntr,wh,conc,negs,pick,term)
                IF(mtp.EQ.6)THEN
                  IF(yes.EQ.1)THEN
                    sng=2*knt
                    dbl=2*knt+1
                  ELSE
                    sng=2*knt+1
                    dbl=2*knt
                  END IF
                  CALL prune(knt,dbl,sng,nkn,ntr,wh,conc,negs,pick,term)
                END IF
                CALL evalgreed(nkn,ntr,conc,pick,negs,term,mtp,n1,
     #                     n2,wh,knt,storage,datri,prtr,mdl,nsep,dcph,
     #                     orders,mtm,rsp,weight,seps,cbetas,xtxsep,nop,
     #                     bestscore,improve,0,bestmove,penalty)
                CALL copytree(ntr,nkn,conc,negs,pick,term,-1,2,1)
              END IF
            END IF
          END DO
        END DO
      END DO
      IF(nop.LT.ntr)THEN
         wh=nop+1
         knt=1
         mtp=0
         neg=0
         DO letter=1,n2
            CALL firstknot(n2,nkn,ntr,wh,conc,negs,pick,term,letter,neg)
            CALL evalgreed(nkn,ntr,conc,pick,negs,term,mtp,n1,
     #                     n2,wh,knt,storage,datri,prtr,mdl,nsep,dcph,
     #                     orders,mtm,rsp,weight,seps,cbetas,xtxsep,nop,
     #                     bestscore,improve,0,bestmove,penalty)
            CALL copytree(ntr,nkn,conc,negs,pick,term,-1,2,1)
         END DO
      END IF
      IF(improve.GE.0)THEN
        CALL copytree(ntr,nkn,conc,negs,pick,term,-1,3,1)
        mtp=bestmove(1)
        wh=bestmove(2)
        knt=bestmove(3)
        CALL evalgreed(nkn,ntr,conc,pick,negs,term,mtp,n1,
     #                     n2,wh,knt,storage,datri,prtr,mdl,nsep,dcph,
     #                     orders,mtm,rsp,weight,seps,cbetas,xtxsep,nop,
     #                     bestscore,improve,1,bestmove,penalty)
        CALL copytree(ntr,nkn,conc,negs,pick,term,-1,1,2)
        CALL copytree(ntr,nkn,conc,negs,pick,term,-1,1,3)
      ELSE
         ntr=-10
      END IF
      END
      SUBROUTINE isallowed(wh,mtp,knt,negs,pick,term,conc,nkn,ntr,yes)
      IMPLICIT NONE
          INTEGER wh,mtp,knt,nkn,ntr,yes
          INTEGER conc(nkn,ntr,3)
          INTEGER negs(nkn,ntr,3)
          INTEGER pick(nkn,ntr,3)
          INTEGER term(nkn,ntr,3)
          INTEGER sibling
          yes=0
          IF(pick(knt,wh,1).EQ.1)THEN
          IF(conc(knt,wh,1).EQ.3)THEN
            IF(mtp.EQ.1)yes=1
            IF(mtp.EQ.3)THEN
              IF(nkn.EQ.1)THEN
                yes=1
              ELSE
                IF (MOD(knt,2).EQ.0) THEN
                  sibling=knt+1
                ELSE
                  sibling=knt-1
                END IF
                IF (conc(sibling,wh,1).EQ.3) yes=1
              END IF
            END IF
            IF(mtp.EQ.4)THEN
              IF(2*knt.LE.nkn)yes=1
            END IF
          ELSE
            IF(mtp.EQ.2)yes=1
            IF(4*knt.LE.nkn)THEN
              IF(mtp.EQ.5)THEN
                IF((conc(2*knt,wh,1).EQ.3).AND.
     #             (conc(2*knt+1,wh,1).EQ.3))yes=1
              END IF
              IF(mtp.EQ.6)THEN
                IF(conc(2*knt,wh,1).EQ.3.AND.conc(4*knt+2,wh,1).EQ.3
     #               .AND.conc(4*knt+3,wh,1).EQ.3) yes=1
                IF (conc(2*knt+1,wh,1).EQ.3.AND.conc(4*knt,wh,1).EQ.3
     #               .AND.conc(4*knt+1,wh,1).EQ.3) yes=2
              END IF
            END IF
          END IF
          END IF
          END
 
        SUBROUTINE evalgreed(nkn,ntr,conc,pick,negs,term,mtp,n1,n2,wh,
     #                      knt,storage,datri,prtr,mdl,nsep,dcph,ordrs,
     #                      mtm,rsp,weight,seps,cbetas,xtxsep,nop,
     #                      bestscore,improve,oncemore,bestmove,penalty)

      IMPLICIT NONE
        INTEGER LGCnknMAX
        PARAMETER (LGCnknMAX  =   128)
        INTEGER nkn,ntr,mtp,n1,n2,wh,knt,nop
        INTEGER conc(nkn,ntr,3)
        INTEGER negs(nkn,ntr,3)
        INTEGER pick(nkn,ntr,3)
        INTEGER term(nkn,ntr,3)
        INTEGER storage(2*ntr*nkn*n1)
        INTEGER prtr(n1,ntr),datri(n2,n1),bestmove(3)
        INTEGER mdl,nsep,dcph(n1),ordrs(n1),mtm,improve,oncemore
        REAL seps(nsep,n1),cbetas(0:(nsep+ntr)),xtxsep(0:nsep,0:nsep)
        REAL weight(n1),bestscore,rsp(n1),penalty
      ! local

        REAL score(3)
        INTEGER ssize,nopold,xstop,reject
        INTEGER nwkv,wkv(LGCnknMAX),l1,l2
        xstop=0
        nopold=nop
        CALL stopper(LGCnknMAX,nkn,"LGCnknMAX","evalgreed()",9,11,
     #          xstop,1)
        CALL gstoring(nkn,ntr,conc,pick,ssize,nop)
        CALL evaluating(wh,knt,mtp,n1,n2,nkn,ntr,conc,term,negs,
     #                  datri,prtr,storage,nwkv,wkv)
        reject=0
        CALL scoring(prtr,rsp,dcph,ordrs,weight,n1,ntr,mdl,nop,wh,
     #             nsep,seps,score,cbetas,reject,xtxsep,mtm,nopold)

       !   take care that the there is real convergence
        l1=0 
        l2=0
        
        IF(cbetas(1).GT. -10000.)l1=1
        IF(cbetas(1).LT. 10000.)l2=1
        IF(l1+l2.EQ.0)reject=1
        IF (reject.NE.1) THEN
        IF(mdl.EQ.2)THEN
          score(1)=score(1)+penalty/(REAL(n1))*ssize
        ELSE 
          score(1)=score(1)+penalty*ssize
        END IF
        IF(oncemore.EQ.0)THEN
          IF(reject.NE.1.AND.score(1).LT.bestscore)THEN
            CALL copytree(ntr,nkn,conc,negs,pick,term,-1,1,3)
            bestscore=score(1)
            bestmove(1)=mtp
            bestmove(2)=wh
            bestmove(3)=knt
            improve=1
          END IF
          CALL restoring(0,wh,n1,nkn,ntr,storage,nwkv,wkv)
          CALL gstoring(nkn,ntr,conc,pick,ssize,nop)
        ELSE
          CALL restoring(1,wh,n1,nkn,ntr,storage,nwkv,wkv)
          CALL gstoring(nkn,ntr,conc,pick,ssize,nop)
        END IF
        ELSE
          CALL copytree(ntr,nkn,conc,negs,pick,term,-1,3,1)
        END IF
        END
      ! *****************************************************************
      ! *****************************************************************

      ! last modification 10/17/02

      SUBROUTINE gstoring(nkn,ntr,conc,pick,ssize,nop)
      IMPLICIT NONE

        ! parameters 
        ! arguments in
          INTEGER nkn,ntr
          INTEGER conc(nkn,ntr,3)
          INTEGER pick(nkn,ntr,3)
        ! local
          INTEGER i,j
        ! arguments out
          INTEGER nop,ssize

      ! determine which moves we can make
        ssize=0
        nop=0
        DO j=1,ntr
          DO i=1,nkn
            IF (pick(i,j,1).EQ.1) THEN
              nop=j
              IF (conc(i,j,1).EQ.3) THEN
                ssize=ssize+1
              END IF
            END IF
          END DO
        END DO
      END 
      ! *****************************************************************
      ! *****************************************************************
      ! this routine scores exponential survival models

      SUBROUTINE exposcore(prtr,rsp,dcph,weight,n1,ntr,
     #                          nop,nsep,seps,score,betas)

      IMPLICIT NONE

        ! arguments in
          INTEGER n1,nop,nsep,ntr,dcph(n1),prtr(n1,ntr)
          REAL rsp(n1),weight(n1),seps(nsep,n1),betas(0:(nsep+ntr))
        ! local
          REAL myrand
          INTEGER i,j
        ! arguments out
          REAL score
          DOUBLE PRECISION myexp,prd

      ! for now, the score is just a random number
c       CALL RANDOM_NUMBER(score)
        score=0.
        DO i=1,n1
          prd=betas(0)
          IF (nsep.GT.0) THEN
            DO j=1,nsep
              prd=prd+betas(j)*seps(j,i)
            END DO
          END IF
          IF (nop.GT.0) THEN
            DO j=1,nop
              prd=prd+betas(nsep+j)*REAL(prtr(i,j))
            END DO
          END IF
          score=score-myexp(prd)*rsp(i)*weight(i)
          IF(dcph(i).EQ.1)score=score+prd*weight(i)
             
        END DO
        score = - score
      END 
      ! *****************************************************************
      ! *****************************************************************
      ! this routine fits exponential survival models
      SUBROUTINE expofit(prtr,rsp,dcph,weight,n1,ntr,nop,
     #                          nsep,seps,score,betas,reject)

      IMPLICIT NONE

        ! arguments in
          INTEGER nop,nsep,ntr,n1
          INTEGER dcph(n1)
          INTEGER prtr(n1,ntr)
          REAL rsp(n1),weight(n1)
          REAL seps(nsep,n1)
        ! local
          INTEGER j,vv,pow2(20),uu,i,zz,k,i2,l,m
          DOUBLE PRECISION s1,s2, myrand,nn(16384),phi(2,14)
          DOUBLE PRECISION dd(2,14),eps,phio(2)
        ! arguments out
          INTEGER reject
          REAL score
          REAL betas(0:(nsep+ntr))
          DOUBLE PRECISION mylog,myexp,pp,vvv

      ! for now, the score is just a random number
        reject=0
        DO j=0,(nsep+ntr)
          betas(j)=0.0
        END DO
c       CALL RANDOM_NUMBER(score)
        pow2(1)=1
        DO j=1,19
           pow2(j+1)=pow2(j)*2
        END DO
        uu=nop+nsep
        vv=pow2(uu+1)
        DO j=1,vv
           nn(j)=0
        END DO
        DO j=1,uu+1
           dd(1,j)=0.0000000000000000000000000000001
           dd(2,j)=0.0000000000000000000000000000001
           phi(1,j)=1.
           phi(2,j)=1.
        END DO
        vvv=0.
        DO i=1,n1
           zz = 0
           IF(nsep.GT.0)THEN
              DO j=1,nsep
                 zz=2*zz+seps(j,i)
              END DO
           END IF
           IF(nop.GT.0)THEN
              DO j=1,nop
                 zz=2*zz+prtr(i,j)
              END DO
           END IF
           nn(zz+1) = nn(zz+1)+weight(i)*rsp(i)
           IF(dcph(i).EQ.1)THEN
              vvv=vvv+weight(i)
              IF(nsep.GT.0)THEN
                 DO j=1,nsep
                    i2=seps(j,i)+1
                    dd(i2,j)=dd(i2,j)+weight(i)
                 END DO
              END IF
              IF(nop.GT.0)THEN
                 DO j=1,nop
                    i2=prtr(i,j)+1
                    dd(i2,j+nsep)=dd(i2,j+nsep)+weight(i)
                 END DO
              END IF
           END IF
        END DO
        score=0.
        IF(uu.EQ.0)THEN
           s1 = 0
           s2 = 0
           DO i=1,n1
              s1 = s1 + weight(i)*rsp(i)
              s2 = s2 + weight(i)*dcph(i)
           END DO
           betas(0)=mylog(s1/s2)
           GOTO 117
        END IF
        DO i=1,100
           eps=0
           DO j=1,uu
              phio(1)=phi(1,j)
              phio(2)=phi(2,j)
              CALL upphi(phi,dd,nn,j,uu,pow2)
              s1 = phi(1,j)-phio(1)
              s2 = phi(2,j)-phio(2)
              eps=eps+SQRT(s1*s1)+SQRT(s2*s2)
           END DO
           IF(eps.LT.0.000006)GOTO 118
        END DO
118     IF(uu.GT.1)THEN
           DO i=2,uu
              phi(1,1) = phi(1,1)*phi(1,i)
              phi(2,1) = phi(2,1)*phi(1,i)
              phi(2,i) = phi(2,i)/phi(1,i)
           END DO
           phi(2,1) = phi(2,1)/phi(1,1)
           betas(0)=mylog(phi(1,1))
           DO i=1,uu
              betas(i)=mylog(phi(2,i))
           END DO
        END IF
        reject=0
        j=pow2(uu+1)
        score=0
        DO i=1,j
          l=i
          pp=betas(0)
          DO k=1,uu
             m=(l+1)/2
             IF(2*m.EQ.l)pp=pp+betas(uu+1-k)
             l=m
          END DO
          score=score-myexp(pp)*nn(i)
        END DO
        DO i=1,uu
          score=score+betas(i)*dd(2,i)
        END DO
        score=score+betas(0)*vvv
        score=-score
117     CONTINUE
      END
      ! *****************************************************************
      ! *****************************************************************
      ! utility for fiting exponential survival models
      SUBROUTINE upphi(phi,dd,nn,j,uu,pow2)
      IMPLICIT NONE
      DOUBLE PRECISION nn(16384),dd(2,14),phi(2,14),mm(16384)
      INTEGER j,uu,pow2(20)
      INTEGER i,k,vv,i2,k2
      vv=pow2(uu+1)
      DO i=1,vv
         mm(i)=nn(i)
      END DO
      IF(j.GT.1)THEN
         DO i=1,(j-1)
            k=pow2(uu-i+1)
            DO i2=1,k
               mm(i2) = mm(i2)*phi(1,i)+mm(i2+k)*phi(2,i)
            END DO
         END DO
      END IF
      k2=pow2(uu-j+1)
      IF(j.LT.uu)THEN
         DO i=(j+1),uu
            k=pow2(uu-i+1)
            DO i2=1,k
               mm(i2) = mm(i2)*phi(1,i)+mm(i2+k)*phi(2,i)
               mm(i2+k2) = mm(i2+k2)*phi(1,i)+mm(i2+k+k2)*phi(2,i)
            END DO
         END DO
      END IF
      phi(1,j)=dd(1,j)/mm(1)
      phi(2,j)=dd(2,j)/mm(k2+1)
      END
