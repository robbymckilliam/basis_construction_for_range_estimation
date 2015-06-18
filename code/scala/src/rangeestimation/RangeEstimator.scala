/*
 * @author Robby McKilliam and Assad Akhlaq
 */

package rangeestimation

import numbers.EuclideanDomain.gcd
import numbers.EuclideanDomain
import numbers.EuclideanDomain.extended_gcd
import numbers.Rational
import numbers.Integer
import numbers.IntegerMatrix
import numbers.RationalMatrix

/** Base trait for estimator of range from multiple frequencies */
trait RangeEstimator {
  
  /** Return range estimate given phase differences y */
  def estimateRange(y : Seq[Double]) : Double
  
}

/** Least squares range estimator. */
class LeastSquaresRangeEstimator(val wavelengths : Seq[Rational]) extends RangeEstimator {
  
  ///The number of wavelengths
  val N = wavelengths.length
  
  ///Least common multiple of the wavelength
  val P = Util.lcm(wavelengths)
  
  ///Integers computed from wavelength after dividing the least common multple P
  val v : Seq[Integer] = wavelengths.map( w => (P/w).n )
  ///Norm of v, i.e., sum of elements squared
  val vnorm = v.foldLeft(Integer.zero) { (s, vi) => s + vi*vi }

  ///The sequence g of gcd's. This is a somewhat obscure way way to write this, but it does 
  ///replicate the recursive definition in the paper.
  val g : Seq[Integer] = v.slice(0,N-1).reverse.foldLeft(List(v.last)) { (list, v) => gcd(list.head,v) :: list }
  
  /// Return the N by N unimodular matrix A_k from the paper
  def A(k : Int) = {
    if( k < 0 || k >= N-1 ) throw new ArrayIndexOutOfBoundsException("The matrix A(k) is only defined for 0 <= k < N-1")
    val c = v(k)/g(k)
    val d = g(k+1)/g(k)
    val (b,minusa) = extended_gcd(c,d)
    val a = -minusa
    def f(m : Int, n : Int) : Integer = { //function returning the elements of the matrix A
      if(m==k && n==k) return c
      else if(m==k+1 && n==k) return d
      else if(m==k && n==k+1) return a
      else if(m==k+1 && n==k+1) return b
      else if(m==n) return Integer.one
      else return Integer.zero
    }
    new IntegerMatrix(f, N,N)
  }
  
  ///The unimodular matrix U given by the product of A matrices
  val U = (1 until N-1).foldLeft(A(0)) { (prod, k) => A(k)*prod } 
  ///The N by N-1 matrix formed by removing the first column from U
  val U2 = U.submatrix(0 until N, 1 until N)
  ///The projection orthogonal to v
  val Q = {
    val vv = RationalMatrix( (m,n) => Rational(v(m)*v(n),Integer.one)/vnorm, N, N)
    RationalMatrix.identity(N) - vv
  }
  ///The generator of the lattice for this range estimator
  val B = Q*U2
  ///The lattice and closest point algorithm for this range estimator
  val Lambda : Util.ClosestPoint = new Util.SphereDecoder(B)
  
  /// Return the range estimate given unwrapping variables z
  def beta(y : Seq[Double], z : IntegerMatrix) = {
    val num = (0 until N).foldLeft(0.0){ (s, n) => s + (y(n) - z(n,0).toDouble)*v(n).toDouble } 
    num / vnorm.toDouble
  }
  
  def estimateRange(y : Seq[Double]) : Double = {
    if( y.size != N ) throw new RuntimeException("Length of y should be N = " + N + " but was " + y.size)
    //paper computes closest point to Qy but the projection isn't necessary because this lattice lies in the space orthogonal to y
    val wseq = Lambda.closestPointIntegers(y)
    val w  = IntegerMatrix.constructColumn(n => wseq(n), N-1)
    val z = U2*w
    val betahat = beta(y,z)
    return P.toDouble*(betahat - betahat.floor)
  }
  
  /** 
   *Variance of this estimator under the assumption that the wrapping variables z are "correct"
   *and the phase noise is i.i.d. with variance sigma2
   */
  final def variance_unwrapping_variables_correct(sigma2 : Double) = {
    val sw = wavelengths.foldLeft(Rational.zero) { (s, w) => s + Rational.one/w/w }
    sigma2/sw.toDouble
  }
    
}

/** 
 *Same as the least squares estimator, but uses Babai's nearest plane algorithm rather than
 *computing the least square estimator by Phost enumerations/sphere decoders.  It's a bit
 *faster, not that it really matters.
 */
class BabaiRangeEstimator(override val wavelengths : Seq[Rational]) extends LeastSquaresRangeEstimator(wavelengths) {
  override val Lambda : Util.ClosestPoint = new Util.Babai(B)
}

/**
 * Single stage CRT.  This is the single stage of the Chinese remainder algorithm described on page 4775 of
 * 
 * Xiao et. al., "Multistage robust Chinese Remainder Theorem", IEEE Trans. on Signal Processing, vol. 62, no. 18, Sep. 15, 2014
 * 
 */
class CRTRangeEstimator(val wavelengths : Seq[Rational]) extends RangeEstimator {
  
  ///The number of wavelengths
  val N = wavelengths.length
  
  ///Least common multiple of the wavelength.  Maximum range of this estimator
  val P = Util.lcm(wavelengths)
  
  /// Numerators and denominators in the rational relationship between the first wavelength and others
  val p = (0 until N).map( i => (wavelengths(i)/wavelengths(0)).n )
  val q = (0 until N).map( i => (wavelengths(i)/wavelengths(0)).d )
  
  /// Least common multiple of the denominators of the wavelength
  val Q = EuclideanDomain.lcm(q)
  
  /// This estimator supposes all the wavelength are integers so scale them to integers
  val M = {
    val rationalM = wavelengths.map(w => Q*w/wavelengths(0))
    rationalM.foreach( r => if(!r.isInteger) throw new RuntimeException("Elements of M should be integers") )
    rationalM.map(_.n)
  }
  
  /// Lowercase m1 from CRT paper (the subscript 1 seemed redundant in the paper)
  val m = (0 until N).map( n=> EuclideanDomain.gcd(M(0),M(n)) )
  
  //These are Gamma1i, Gammai1, and barGammai in Xiao et. al.
  val A = (0 until N).map( n => M(0)/m(n) )
  val B = (0 until N).map( n => M(n)/m(n) )
  val C = (0 until N).map( n => Util.inverse_mod(A(n), B(n)) )
  
  def estimateRange(y : Seq[Double]) : Double = {
    if( y.size != N ) throw new RuntimeException("Length of y should be N = " + N + " but was " + y.size)
    if( N == 1) return P.toDouble*(y(0) - (y(0)/P.toDouble).floor); //special case when N is 1
    val r = (0 until N).map( i => M(i).toDouble * y(i) ) //these are the observations normalised as in the CRT paper
    val q = (0 until N).map( i => ((r(i) - r(0))/m(i).toDouble).round ) //step 2 of CRT paper
    val xi = (0 until N).map( i => (q(i)*C(i)) mod B(i) ) //step 3 of CRT paper
    val n1 = Util.chinese_remainder(xi.slice(1,N), B.slice(1,N)) //step 4, solve generalised Chinese Remainder theorem
    val n = (0 until N).map( i => Rational(n1*A(i) - q(i),B(i)) )
    //if( Rational(n1,1) != n(0) ) throw new RuntimeException("n1 and computed n are not the same. This shouldn't be possible.")
    val hatN = (0 until N).foldLeft(0.0)( (s, i) => s + (n(i)*M(i)).toDouble + r(i) )
    val rhat = wavelengths(0).toDouble*hatN/N/Q.toDouble
    return rhat - P.toDouble*(rhat/P.toDouble).floor;
  }
  
}

/**
 * The two stage CRT range estimator from
 * 
 * Xiao et. al., "Multistage robust Chinese Remainder Theorem", IEEE Trans. on Signal Processing, vol. 62, no. 18, Sep. 15, 2014
 * 
 */
class TwoStageCRTRangeEstimator(val wavelengths : Seq[Rational], val subsets : Seq[Seq[Int]])  extends RangeEstimator {
  
  //the total number of wavelengths
  val N = wavelengths.size
  if(N != subsets.foldLeft(0)( (s, w) => s + w.size )) throw new ArrayIndexOutOfBoundsException("Number of indices in subsets must be equal to the number of wavelengths N = " + N)
  
  //first level of estimators run
  val crtestimators = subsets.map( subset => new CRTRangeEstimator( subset.map( i => wavelengths(i)) ) )
  
  //psuedowavelength of the first level estimator
  val psuedowavelengths = crtestimators.map( e => e.P )
  //top level estimator using the psuedowavelengths
  val toplevelestimator = new CRTRangeEstimator(psuedowavelengths)
  
  /// Return estimates of the range  from each of the CRT estimators in the first stage
  def pseudoranges(y : Seq[Double]) : Seq[Double] = {
    return subsets.indices.map( si => crtestimators(si).estimateRange( subsets(si).map( i => y(i)) ) )
  }
  
  /// Returns psudeo phases from each of the CRT estimators in the first stage
  def psuedophases(y : Seq[Double]) : Seq[Double] = {
    val r = pseudoranges(y)
    return r.indices.map( i => Util.fracpart( r(i)/psuedowavelengths(i).toDouble ) )
  }
  
  def estimateRange(y : Seq[Double]) : Double = {
    if( y.size != N ) throw new RuntimeException("Length of y should be N = " + N + " but was " + y.size)
    val p = psuedophases(y)
    return toplevelestimator.estimateRange(p)
  }
  
}
