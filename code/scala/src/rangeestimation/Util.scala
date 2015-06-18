/*
 * Miscellaeous useful functions.
 * @author Robby McKilliam
 */

package rangeestimation

import numbers.RationalMatrix
import numbers.Rational
import numbers.Integer
import org.mckilliam.lattices.Lattice

object Util {
  
  def sqr(x : Rational) = x*x
  
  /** Computes the centered fractional part of x */
  def fracpart(x : Double) = x - x.round
  
  /** Least common multiple of two rational numbers */
  def lcm(a : Rational, b : Rational) : Rational = {
    val g = a.n * b.d
    val f = b.n * a.d
    val k = numbers.EuclideanDomain.lcm(g,f)
    return Rational(k,a.d*b.d)
  }
  
  /** Least common multiple of a sequence of rational numbers */
  def lcm(s : Seq[Rational]) : Rational = s.reduceLeft { (g,v) => lcm(g,v) }
  
  /** Return the (a) multiplicative inverse of a mod b */
  def inverse_mod(a : Integer, b : Integer) : Integer = {
    if( numbers.EuclideanDomain.gcd(a,b) != Integer.one ) throw new RuntimeException("a and b must be relatively prime for inverse of a mod b to exist")
    val (s,t) = numbers.EuclideanDomain.extended_gcd(a,b)
    return s
  }
  
  /** 
   *Solves x = a(i) mod m(i) for all i.  Uses the algorithm implied by Theorem 1 of
   * Oystein Ore, "The general Chinese remainder theorem", The American Mathematical Monthly, Vol. 59, No. 6 (Jun. - Jul., 1952), pp. 365-370 
   */ 
  def chinese_remainder(a : Seq[Integer], m : Seq[Integer]) : Integer = {
    if(a.size != m.size) throw new RuntimeException("a and m must contain the same number of elements.")
    val M = numbers.EuclideanDomain.lcm(m)
    val B = m.map( x => M/x )
    if( numbers.EuclideanDomain.gcd(B) != Integer.one ) throw new RuntimeException("Somehow the gcd of B is not relatively prime.  This shouldn't be possible.")
    val c = numbers.EuclideanDomain.extended_gcd(B)
    return a.indices.foldLeft(Integer.zero)( (s, i) => s + a(i)*c(i)*B(i) )
  }
  
  /** Convert ScalarNumber RationalMatrix to Jama matrix */
  def RationalMatrixToJama(B : RationalMatrix) : Jama.Matrix = {
    val J = new Jama.Matrix(B.M,B.N)
    for( m <- 0 until B.M ) for( n <- 0 until B.N ) J.set(m,n,B(m,n).toDouble)
    return J
  }
   
  trait ClosestPoint {
    val cvp : org.mckilliam.lattices.ClosestVectorInterface
    def closestPointIntegers(y : Seq[Double]) : Seq[Integer] = {
      cvp.closestPoint(y.toArray)
      val w = cvp.getIndex()
      w.map( i => Integer(i.toInt) )
    }
  }
  
  class SphereDecoder(B : RationalMatrix) extends ClosestPoint {
    override val cvp = new org.mckilliam.lattices.cvp.SphereDecoderSchnorrEuchner(new Lattice(RationalMatrixToJama(B))) 
  }
  
  class Babai(B : RationalMatrix) extends ClosestPoint {
    override val cvp = new org.mckilliam.lattices.cvp.Babai(new Lattice(RationalMatrixToJama(B))) 
  }

}
