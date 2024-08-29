      ! This particular routine is intended to be rewritten by users

      ! this subroutine calculates the score under my own model
      ! last modification 07/02/02


      ! you can write your own scoring function!
      !
      ! essentially you need to provide two routines -
      ! (1) a routine which fits your model: it provides a coefficient
      ! (beta) for each of the logic trees and provides a score
      ! of how good the model is. Low scores are good. (So add a
      ! minus sign if your score is a log-likelihood.)
      ! (2) a routine which - given the betas - provides the score
      ! of your model. [If you don't use cross-validation, this
      ! second routine is not needed.]
      !
      ! We pass on the variables!
      !
      ! Below is a list of variables that are passed on. Most of them
      ! are as you expect - response, predictors (binary ones and
      ! continuous ones), number of cases, number of predictors.
      ! In addition there are two columns - dcph and weight - that
      ! can either be used to pass on an auxilliary variable for
      ! each case (discrete for dcph and continous for weight), or
      ! even some overall auxilliary variables - as these numbers are
      ! not used anywhere else.
      !
      ! If you do not need any of the variables - just ignore them!

        ! prtr     the predictions of the logic trees in the current
        !          model: this is an integer matrix of size
        !          n1 times ntr - although only the first nop columns
        !          contain usefule information
        ! rsp      the response variable: this is a real (single
        !          precision) vector of length n1
        ! dcph     censor times: this is an integer vector of length n1
        !          this could be used as an auxilliary (integer)
        !          vector - as it is just passed on.  (There is no
        !          check that this is a 0/1 variable, unless you use
        !          a proportional hazards model.) (For example,
        !          you could use this to pass on cluster membership.)
        ! weight   weights for the cases
        !          this is a real vector of length n1.
        !          this could be used as an auxilliary (real)
        !          vector - as it is just passed on.  (There is no
        !          check that these numbers are positive.)
        ! ordrs    the order (by response size) of the cases
        !          This is an integer vector of length n1. For
        !          the case with the smallest response this one is
        !          1, for the second smallest 2, and so on. Ties
        !          are resolved arbitrary. Use it if you wish!
        ! n1       the total number of cases in the data
        ! ntr      the number of logic trees ALLOWED in the tree
        !          (you probably don't need this one)
        ! nop      the number of logic trees in the CURRENT model
        ! wh       the index of the tree that has been edited in the
        !          last move - i.e. the column of prtr that has
        !          changes since the last call
        ! nsep     number of varaibles that get fit a separate
        !          parameter
        ! seps     array of the above variables - this is a single
        !          precision matrix of size n1 times nsep

      ! For Myownfitting:

        ! You are allowed to change the values of dcph, weight and
        ! ordrs     not of any of the other variables
        ! You should return
        ! betas     a vector of parameters of the model that you
        !           fit.
        !           betas(0) should be the parameter for the intercept
        !           betas(1:nsep) should be the parameters for the
        !                   continuous variables in seps
        !           betas((nsep+1):(nsep+nop)) should be the
        !                   parameters for the binary trees in prtr
        !           if you have more parameters, use dcph, weight or
        !           orders ==> these variables will not be printed
        !           however.
        ! score     whatever score you assign to your model
        !           small should be good (i.e. deviance or
        !           -log.likelihood)
        ! reject    an indicator whether or not to reject the proposed
        !           move *regardless* of the score (for example when an
        !           iteration necessary to determine the score failed
        !           to converge (0 = move is OK ; 1 = reject move)
        !           set this one to 0 if there is no such condition

      ! For Myownscoring:

        ! Additional input is
        !
        ! betas     the coefficients
        !
        ! You should return
        !
        ! score     whatever score you assign to your model
        !           small should be good (i.e. deviance or
        !           -log.likelihood)
        !           If the model "crashes", you should simply return a
        !           very large number

      ! Note: both subroutines should work if nop and or nsep is zero.
      ! while we try to prevent that models are singular, it is possible
      ! that for your model a single or degenerate model is passed on
      ! for evaluation. For Myownfitting you can pass the model
      ! back with reject=1, for Myownscoring you can pass it on with
      ! a very large value for score.

        ! Below are the frames to put these scoring functions in.
        ! currently they just yield random numbers as the fitted models.


      SUBROUTINE Myownscoring(prtr,rsp,dcph,ordrs,weight,n1,ntr,
     #                          nop,wh,nsep,seps,score,betas)

      IMPLICIT NONE

        ! arguments in
          INTEGER n1,nop,nsep,ntr,wh
          INTEGER dcph(n1),ordrs(n1)
          INTEGER prtr(n1,ntr)
          REAL rsp(n1),weight(n1)
          REAL seps(nsep,n1)
          REAL betas(0:(nsep+ntr))
        ! arguments out
          REAL score

      ! line below only to prevent warnings during compilation
        score=prtr(1,1)*rsp(1)*dcph(1)*ordrs(1)*weight(1)
        score=score+nop+wh+seps(1,1)+betas(0)
      ! for now, the score is 0
        score=0.

      END
      SUBROUTINE Myownfitting(prtr,rsp,dcph,ordrs,weight,n1,ntr,nop,
     #                          wh,nsep,seps,score,betas,reject)

      IMPLICIT NONE

        ! arguments in
          INTEGER nop,nsep,ntr,wh,n1
          INTEGER dcph(n1),ordrs(n1)
          INTEGER prtr(n1,ntr)
          REAL rsp(n1),weight(n1)
          REAL seps(nsep,n1)
        ! local
          INTEGER j
        ! arguments out
          INTEGER reject
          REAL score
          REAL betas(0:(nsep+ntr))

      ! line below only to prevent warnings during compilation
        score=prtr(1,1)*rsp(1)*dcph(1)*ordrs(1)*weight(1)
        score=score+nop+wh+seps(1,1)+betas(0)
      ! for now, the score is 0
        reject=0
        DO j=0,(nsep+ntr)
          betas(j)=0.0
        END DO
        score=0.

      END
