/*
 * Function for selecting good sets of frequencies for range estimation.
 * @author Robby McKilliam 12/1/2015
 */

package rangeestimation

import numbers.Rational
import numbers.Integer
import Util.sqr
import numbers.EuclideanDomain.gcd

/** Calculates an optimised set of N wavelengths for range estimation */
class WavelengthOptimiser(N : Int, rmin : Rational, lambdamin : Rational, lambdamax : Rational, beta : Rational) {
  
  val B = sqr(lambdamin*lambdamax)/(sqr(lambdamin) + sqr(lambdamax)*(N-1))
  val d = (List(rmin/lambdamax, lambdamin/(lambdamax-lambdamin)).max).ceil
  protected var lambdahat = (1 to N-1).map{ i => lambdamax*d/(d+1)}.toArray //stores the currently best found wavelength
  protected var Ltilde = L1( lambdahat ) //initial bound on objective function
  protected var pm = Array.fill(N-1)(Integer.one) //memory for p's
  protected var qm = Array.fill(N-1)(Integer.one) //memory of q's
  psearch(0) //start the search
  
  /** Return optimised wavelengths */
  val wavelengths = lambdamax :: lambdahat.toList
    
  /** First recursive function enumerates over candidate numerators, i.e., the p's */
  protected def psearch(n : Int) : Unit = {
    while( sqr(pm(n)) <= (Ltilde - beta*B)/N ) {
      qm(n) = pm(n)
      qsearch(n)
      pm(n) = pm(n) + 1
    }
  }
  
  /** Second recurive function enumerates over the denominators, i.e., the q's */
  protected def qsearch(n : Int) : Unit = {
    while( qm(n) <= pm(n)*lambdamax/lambdamin ) {
      if(gcd(pm(n),qm(n)) == Integer.one) { //ignore if this p(n) and q(n) are not relatively prime
        if(n < N-2) { //we aren't at the end yet so keep recursing
          pm(n+1) = pm(n)
          psearch(n+1)
        }
        else if( Rational(numbers.EuclideanDomain.lcm(pm)) >= rmin/lambdamax && L2(pm,qm) < Ltilde ) {
          //this is the currently best found set of wavelengths so update Ltilde and lambdahat
          Ltilde = L2(pm,qm)
          for( i <- pm.indices) lambdahat(i) = lambdamax*Rational(pm(i),qm(i))
        }
      }
      qm(n) = qm(n) + 1
    }
  }
  
  /** Objective function assuming that the first wavelength lambda1 is lambdamax */
  def L1(w : Seq[Rational]) = {
    if(w.size != N-1) throw new ArrayIndexOutOfBoundsException("Number of wavelengths must be N = " + (N-1))
    val P = Util.lcm(lambdamax :: w.toList) //lcm of a sequence of rationals
    val sumw = w.foldLeft(1/sqr(lambdamax)){ (sum, v) => sum + 1/sqr(v) }
    P*P*sumw + beta/sumw
  }
   
  /** Transformed objective function with one real/rational parameter w1 and 2(N-1) integer parameters */
  def L2(p : Seq[Integer], q : Seq[Integer]) = {
    if(p.size != N-1 || q.size != N-1) throw new ArrayIndexOutOfBoundsException("p and q must have size N-1 = " + (N-1))
    val Q = numbers.EuclideanDomain.lcm(p) //this is the lcm of a sequence of integers
    val D = (0 until N-1).foldLeft(Rational.one){ (sum, n) => sum + Rational(q(n)*q(n),p(n)*p(n)) }
    Q*Q*D + beta*lambdamax*lambdamax/D
  }
  
  /** Compute objective function directly from wavelengths. Just here for testing purposes */
  def L(w : Seq[Rational]) : Rational = {
    if(w.size != N) throw new ArrayIndexOutOfBoundsException("Number of wavelengths must be N = " + N)
    val P = Util.lcm(w) //lcm of a sequence of rationals
    val sumw = w.foldLeft(Rational.zero){ (sum, v) => sum + 1/v/v }
    P*P*sumw + beta/sumw
  }
  
  
  
}
