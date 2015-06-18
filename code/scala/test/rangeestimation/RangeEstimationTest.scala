/*
 * @author Robby McKilliam.
 */

package rangeestimation

import org.junit._
import Assert._
import scala.math.round
import numbers.Rational
import numbers.RationalMatrix
import numbers.Integer
import numbers.IntegerMatrix
import numbers.EuclideanDomain.gcd
import org.mckilliam.distributions.Gaussian
import org.mckilliam.distributions.Uniform

class RangeEstimationTest {

  def fracpart(x : Double) = x - round(x)
  val tol = 1e-7
  
  @Test
  def testP = {
    println("Testing P")
    val l1 = List(Rational(1,5),Rational(1,3),Rational(1,1),Rational(1,2))
    val e1 = new LeastSquaresRangeEstimator( l1 )
    assertTrue(e1.P == Rational(1,1))
  }
  
  @Test
  def testg = {
    println("Testing the sequence g")
    val l1 = List(Rational(1,5),Rational(1,3),Rational(1,1),Rational(1,2), Rational(2,3), Rational(3,7))
    val e1 = new LeastSquaresRangeEstimator( l1 )
    for( k <- 0 until e1.N ) assertTrue( e1.g(k) == gcd(e1.v.slice(k,e1.N)) )
    for(  k <- 0 until e1.N-1 ) assertTrue( Rational( e1.g(k+1), e1.g(k)).d == Integer.one ) //check that g(k) divides g(k+1)
    for(  k <- 0 until e1.N ) assertTrue( Rational( e1.v(k), e1.g(k)).d == Integer.one ) //check that g(k) divides v(k)
    for( k <- 0 until e1.N-1 ) {
      val c = Rational( e1.g(k+1), e1.g(k)).n
      val d = Rational( e1.v(k), e1.g(k)).n
      assertTrue( c == e1.g(k+1)/e1.g(k) ) //check that direct integer division works
      assertTrue( d == e1.v(k)/e1.g(k) )
      assertTrue( gcd(c,d) == Integer.one ) //check that these integers are relatively prime
    }
  }
  
  @Test
  def testAdet = {
    println("Test that A(k) has determinant 1, i.e. that it is unimodular")
    val l1 = List(Rational(1,5),Rational(1,3),Rational(1,1),Rational(1,2), Rational(2,3), Rational(3,7))
    val e1 = new LeastSquaresRangeEstimator( l1 )
    for( k <- 0 until e1.N-1 ) {
      println(e1.A(k).det)
      assertTrue( e1.A(k).det == Integer.one ) 
    }
  }
  
  @Test 
  def testArecursivevFormula = {
    println("Test the recursive formula v_k = A_k v_{k-1}")
    val l1 = List(Rational(1,5),Rational(1,3),Rational(1,1),Rational(1,2), Rational(2,3), Rational(3,7))
    val e1 = new LeastSquaresRangeEstimator( l1 )
    
    //returns the vector v_k from the paper
    def v(k : Int) : IntegerMatrix = {
      def f(n : Int) : Integer = {
        if(n <= k) return e1.v(n)
        else if(n==k+1) return e1.g(n)
        else return Integer.zero
      }
      return IntegerMatrix.constructColumn(f, e1.N)
    }
    
    { //assert that v(N-1) is v
      val vNm1 = v(e1.N-1)
      for( i <- 0 until e1.N) assertTrue(e1.v(i) == vNm1(i)) 
    }
    
    //asert that v(k+1) = A(k+1) v(k)
    for( k <- 0 until e1.N-2){ 
      val vk = v(k)
      val vk1 = v(k+1)
      val Ak1 = e1.A(k+1)
      val comp = Ak1*vk
      for( i <- 0 until e1.N) assertTrue(vk1(i) == comp(i)) 
    }
    
  }
  
  @Test
  def testU = {
    println("Test that the matrix U is unimoular and has first column equal to v")
    val l1 = List(Rational(1,5),Rational(1,3),Rational(1,1),Rational(1,2), Rational(2,3), Rational(3,7))
    val e1 = new LeastSquaresRangeEstimator( l1 )
    println(e1.U)
    println(e1.U.det)
    for( i <- 0 until e1.N) assertTrue(e1.v(i) == e1.U(i,0) ) //test first column is v
    assertTrue( e1.U.det == Integer.one ) //test determinant is 1
  }
  
  @Test
  def testQ = {
    println("Test the projection matrix Q")
    val l1 = List(Rational(1,5),Rational(1,3),Rational(1,1),Rational(1,2), Rational(2,3), Rational(3,7))
    val e1 = new LeastSquaresRangeEstimator( l1 )
    println(e1.Q)
    val v = RationalMatrix( (m,n) => e1.v(m), e1.N, 1)
    assertTrue( (e1.Q * v).squaredFrobeniusNorm == Rational.zero )
  }
  
  @Test
  def testEstimateNoiseLess = {
    println("Test estimate method of least squares with no noise")
    val l1 = List(Rational(1,5),Rational(1,3),Rational(1,1),Rational(1,2), Rational(2,3), Rational(3,7))
    val e1 = new LeastSquaresRangeEstimator(l1)
    val P = e1.P.toDouble
    val r0 = Uniform.constructFromMinMax(0,P).noise
    val y = l1.map( l => Util.fracpart(r0/l.toDouble) )
    val rhat = e1.estimateRange(y)
    println(r0,rhat)
    assertTrue( (Util.fracpart((r0 - rhat)/P)*P).abs < tol  )
  }
  
    @Test
  def testEstimateWithNoise = {
    println("Test estimate method of least squares with a small amount of noise")
    val l1 = List(Rational(1,5),Rational(1,3),Rational(1,1),Rational(1,2), Rational(2,3), Rational(3,7))
    val e1 = new LeastSquaresRangeEstimator(l1)
    val P = e1.P.toDouble
    val r0 = Uniform.constructFromMinMax(0,P).noise
    val sigma = 0.001
    val Phi = new Gaussian(0,sigma*sigma).noise
    val y = l1.map( l => Util.fracpart(r0/l.toDouble + Phi))
    val rhat = e1.estimateRange(y)
    println(r0,rhat)
    assertTrue( (Util.fracpart((r0 - rhat)/P)*P).abs < 30*sigma  )
  }
  
  @Test
  def testEstimateBabaiNoiseLess = {
    println("Test estimate method of Babai least squares with no noise")
    val l1 = List(Rational(1,5),Rational(1,3),Rational(1,1),Rational(1,2), Rational(2,3), Rational(3,7))
    val e1 = new BabaiRangeEstimator(l1)
    val P = e1.P.toDouble
    val r0 = Uniform.constructFromMinMax(0,P).noise
    val y = l1.map( l => Util.fracpart(r0/l.toDouble) )
    val rhat = e1.estimateRange(y)
    println(r0,rhat)
    assertTrue( (Util.fracpart((r0 - rhat)/P)*P).abs < tol  )
  }
  
  @Test
  def testEstimateBabaiWithNoise = {
    println("Test estimate method of Babai least squares with a small amount of noise")
    val l1 = List(Rational(1,5),Rational(1,3),Rational(1,1),Rational(1,2), Rational(2,3), Rational(3,7))
    val e1 = new BabaiRangeEstimator(l1)
    val P = e1.P.toDouble
    val r0 = Uniform.constructFromMinMax(0,P).noise
    val sigma = 0.001
    val Phi = new Gaussian(0,sigma*sigma).noise
    val y = l1.map( l => Util.fracpart(r0/l.toDouble + Phi))
    val rhat = e1.estimateRange(y)
    println(r0,rhat)
    assertTrue( (Util.fracpart((r0 - rhat)/P)*P).abs < 30*sigma  )
  }
  
  @Test 
  def testBabaiAndSphereDecoderDifferent = {
    println("Check that Babai's algorithm and least squares are different even when integers are relatively prime.")
    val tol = 1e-6
    val lambda = List(Rational(2),Rational(3),Rational(5),Rational(7))
    val babai = new BabaiRangeEstimator(lambda)
    val leastsq = new LeastSquaresRangeEstimator(lambda)
    val P = babai.P.toDouble
    val maxiters = 10000
    def testf(iters : Int) : Boolean = {
      val r0 = Uniform.constructFromMinMax(0,P).noise
      val Phi = new Gaussian(0,1000*1000).noise
      val y = lambda.map( l => Util.fracpart(r0/l.toDouble + Phi))
      val rbabai = babai.estimateRange(y)
      val rls = leastsq.estimateRange(y)
      if( (rbabai - rls).abs > tol) return true
      else if(iters > maxiters) return false
      else return testf(iters+1)
    }
    assertTrue(testf(1))
  }
  
  @Test
  def testCRTMareIntegers = {
    println("Test that the normalised wavelengths M for the CRT estimator turn out to be integers");
    {
      val l1 = List(Rational(1,5),Rational(1,3),Rational(1,1),Rational(1,2), Rational(2,3), Rational(3,7))
      val crtest = new CRTRangeEstimator(l1) //Constructor will throw an exception if any M's aren't integers
    }
    {
      val l1 = List(Rational(2),Rational(3),Rational(5),Rational(7))
      val crtest = new CRTRangeEstimator(l1) //Constructor will throw an exception if any M's aren't integers
    }
  }
  
  @Test
  def testCRTNoiseLess = {
    println("Test CRT estimate with no noise");
        {
      //test with one
      val l1 = List(Rational(1,3))
      val e1 = new CRTRangeEstimator(l1)
      val P = Util.lcm(l1).toDouble
      val r0 = P/2.0
      val y = l1.map( l => Util.fracpart(r0/l.toDouble) )
      val rhat = e1.estimateRange(y)
      println(r0,rhat)
      assertTrue( (Util.fracpart((r0 - rhat)/P)*P).abs < tol  )
    }
    {
      //test with two
      val l1 = List(Rational(1,5),Rational(1,3))
      val e1 = new CRTRangeEstimator(l1)
      val P = Util.lcm(l1).toDouble
      val r0 = P/2.0
      val y = l1.map( l => Util.fracpart(r0/l.toDouble) )
      val rhat = e1.estimateRange(y)
      println(r0,rhat)
      assertTrue( (Util.fracpart((r0 - rhat)/P)*P).abs < tol  )
    }
    {
      //test with two reordered
      val l1 = List(Rational(1,3), Rational(1,5))
      val e1 = new CRTRangeEstimator(l1)
      val P = Util.lcm(l1).toDouble
      val r0 = P/2.0
      val y = l1.map( l => Util.fracpart(r0/l.toDouble) )
      val rhat = e1.estimateRange(y)
      println(r0,rhat)
      assertTrue( (Util.fracpart((r0 - rhat)/P)*P).abs < tol  )
    }
    {
      val l1 = List(Rational(1,5),Rational(1,3),Rational(1,1),Rational(1,2), Rational(2,3), Rational(3,7))
      val e1 = new CRTRangeEstimator(l1)
      val P = Util.lcm(l1).toDouble
      val r0 = P/2.0
      val y = l1.map( l => Util.fracpart(r0/l.toDouble) )
      val rhat = e1.estimateRange(y)
      println(r0,rhat)
      assertTrue( (Util.fracpart((r0 - rhat)/P)*P).abs < tol  )
    }
    {
      val l1 = List(Rational(2,3), Rational(3,7),Rational(1,5),Rational(1,3),Rational(1,1),Rational(1,2))
      val e1 = new CRTRangeEstimator(l1)
      val P = Util.lcm(l1).toDouble
      val r0 = P/2.0
      val y = l1.map( l => Util.fracpart(r0/l.toDouble) )
      val rhat = e1.estimateRange(y)
      println(r0,rhat)
      assertTrue( (Util.fracpart((r0 - rhat)/P)*P).abs < tol  )
    }
    {
      val l1 = List(Rational(2),Rational(3),Rational(5),Rational(7)) //list of wavelenghts
      val e1 = new CRTRangeEstimator(l1)
      val P = Util.lcm(l1).toDouble
      val r0 = 20*scala.math.Pi
      val y = l1.map( l => Util.fracpart(r0/l.toDouble) )
      val rhat = e1.estimateRange(y)
      println(r0,rhat)
      assertTrue( (Util.fracpart((r0 - rhat)/P)*P).abs < tol  )
    }
    {
      //test with two integers
      val l1 = List(Rational(1),Rational(6)) //list of wavelenghts
      val e1 = new CRTRangeEstimator(l1)
      val P = Util.lcm(l1).toDouble
      val r0 = 20*scala.math.Pi
      val y = l1.map( l => Util.fracpart(r0/l.toDouble) )
      val rhat = e1.estimateRange(y)
      println(r0,rhat)
      assertTrue( (Util.fracpart((r0 - rhat)/P)*P).abs < tol  )
    }
    {
      val l1 = List(Rational(2),Rational(42,13),Rational(210,37),Rational(7)) //list of wavelenghts
      val e1 = new CRTRangeEstimator(l1)
      val P = Util.lcm(l1).toDouble
      val r0 = 20*scala.math.Pi
      val y = l1.map( l => Util.fracpart(r0/l.toDouble) )
      val rhat = e1.estimateRange(y)
      println(r0,rhat)
      assertTrue( (Util.fracpart((r0 - rhat)/P)*P).abs < tol  )
    }
    {
      val l1 = List(Rational(210,79),Rational(210,61),Rational(210,41),Rational(210,31)) //list of wavelenghts
      val e1 = new CRTRangeEstimator(l1)
      val P = Util.lcm(l1).toDouble
      val r0 = 20*scala.math.Pi
      val y = l1.map( l => Util.fracpart(r0/l.toDouble) )
      val rhat = e1.estimateRange(y)
      println(r0,rhat)
      assertTrue( (Util.fracpart((r0 - rhat)/P)*P).abs < tol  )
    }
  }
  
  @Test
  def testCRTWithNoise = {
    println("Test CRT estimate with a small amount of noise");
    { //test with one
    val l1 = List(Rational(2,3))
    val e1 = new CRTRangeEstimator(l1)
    val P = Util.lcm(l1).toDouble
    //println("P = " + P)
    //val r0 = Uniform.constructFromMinMax(0,P).noise
    val r0 = P/2.0
    val sigma = 0.001
    val Phi = new Gaussian(0,sigma*sigma).noise
    val y = l1.map( l => Util.fracpart(r0/l.toDouble + Phi))
    val rhat = e1.estimateRange(y)
    println(r0,rhat, (Util.fracpart((r0 - rhat)/P)*P).abs)
    assertTrue( (Util.fracpart((r0 - rhat)/P)*P).abs < 30*sigma  )
    }
    { //test with two
    val l1 = List(Rational(2,3), Rational(3,7))
    val e1 = new CRTRangeEstimator(l1)
    val P = Util.lcm(l1).toDouble
    //val r0 = Uniform.constructFromMinMax(0,P).noise
    val r0 = P/2.0
    val sigma = 0.001
    val Phi = new Gaussian(0,sigma*sigma).noise
    val y = l1.map( l => Util.fracpart(r0/l.toDouble + Phi))
    val rhat = e1.estimateRange(y)
    println(r0,rhat)
    assertTrue( (Util.fracpart((r0 - rhat)/P)*P).abs < 30*sigma  )
    }
    { //test with two reordered
    val l1 = List(Rational(3,7),Rational(2,3))
    val e1 = new CRTRangeEstimator(l1)
    val P = Util.lcm(l1).toDouble
    //val r0 = Uniform.constructFromMinMax(0,P).noise
    val r0 = P/2.0
    val sigma = 0.001
    val Phi = new Gaussian(0,sigma*sigma).noise
    val y = l1.map( l => Util.fracpart(r0/l.toDouble + Phi))
    val rhat = e1.estimateRange(y)
    println(r0,rhat)
    assertTrue( (Util.fracpart((r0 - rhat)/P)*P).abs < 30*sigma  )
    }
    {
    val l1 = List(Rational(1,5),Rational(1,3),Rational(1,1),Rational(1,2), Rational(2,3), Rational(3,7))
    val e1 = new CRTRangeEstimator(l1)
    val P = Util.lcm(l1).toDouble
    //val r0 = Uniform.constructFromMinMax(0,P).noise
    val r0 = P/2.0
    val sigma = 0.001
    val Phi = new Gaussian(0,sigma*sigma).noise
    val y = l1.map( l => Util.fracpart(r0/l.toDouble + Phi))
    val rhat = e1.estimateRange(y)
    println(r0,rhat)
    assertTrue( (Util.fracpart((r0 - rhat)/P)*P).abs < 30*sigma  )
    }
    { //test reordered
    val l1 = List( Rational(3,7),Rational(1,5),Rational(1,1),Rational(1,2), Rational(2,3),Rational(1,3))
    val e1 = new CRTRangeEstimator(l1)
    val P = Util.lcm(l1).toDouble
    //val r0 = Uniform.constructFromMinMax(0,P).noise
    val r0 = P/2.0
    val sigma = 0.001
    val Phi = new Gaussian(0,sigma*sigma).noise
    val y = l1.map( l => Util.fracpart(r0/l.toDouble + Phi))
    val rhat = e1.estimateRange(y)
    println(r0,rhat)
    assertTrue( (Util.fracpart((r0 - rhat)/P)*P).abs < 30*sigma  )
    }
    {
    val l1 = List(Rational(2),Rational(3),Rational(5),Rational(7)) //list of wavelenghts
    val e1 = new CRTRangeEstimator(l1)
    val P = Util.lcm(l1).toDouble
    //val r0 = Uniform.constructFromMinMax(0,P).noise
    val r0 = P/2.0
    val sigma = 0.001
    val Phi = new Gaussian(0,sigma*sigma).noise
    val y = l1.map( l => Util.fracpart(r0/l.toDouble + Phi))
    val rhat = e1.estimateRange(y)
    println(r0,rhat)
    assertTrue( (Util.fracpart((r0 - rhat)/P)*P).abs < 30*sigma  )
    }
    {
    val l1 = List(Rational(2),Rational(42,13),Rational(210,37),Rational(7)) //list of wavelenghts
    val e1 = new CRTRangeEstimator(l1)
    val P = Util.lcm(l1).toDouble
    //val r0 = Uniform.constructFromMinMax(0,P).noise
    val r0 = P/2.0
    val sigma = 0.001
    val Phi = new Gaussian(0,sigma*sigma).noise
    val y = l1.map( l => Util.fracpart(r0/l.toDouble + Phi))
    val rhat = e1.estimateRange(y)
    println(r0,rhat)
    assertTrue( (Util.fracpart((r0 - rhat)/P)*P).abs < 30*sigma  )
    }
    {
    val l1 = List(Rational(210,79),Rational(210,61),Rational(210,41),Rational(210,31)) //list of wavelenghts
    val e1 = new CRTRangeEstimator(l1)
    val P = Util.lcm(l1).toDouble
    //val r0 = Uniform.constructFromMinMax(0,P).noise
    val r0 = P/2.0
    val sigma = 0.001
    val Phi = new Gaussian(0,sigma*sigma).noise
    val y = l1.map( l => Util.fracpart(r0/l.toDouble + Phi))
    val rhat = e1.estimateRange(y)
    println(r0,rhat)
    assertTrue( (Util.fracpart((r0 - rhat)/P)*P).abs < 30*sigma  )
    }
  }
  
  
  @Test
  def TwoStageCRTFirstStageRangeEstimator = {
    println("Test first stage of two stage CRT estimator");
    {
      val lambda = List(Rational(1,5),Rational(1,3),Rational(1,1),Rational(1,2), Rational(2,3), Rational(3,7))
      val subsets = List(List(0,1,2), List(3,4,5)) 
      val e = new TwoStageCRTRangeEstimator(lambda, subsets)
      val e1 = new CRTRangeEstimator(List(Rational(1,5),Rational(1,3),Rational(1,1)))
      val e2 = new CRTRangeEstimator(List(Rational(1,2), Rational(2,3), Rational(3,7)))
      val P = Util.lcm(lambda).toDouble
      //val r0 = Uniform.constructFromMinMax(0,P).noise
      val r0 = P/2.0
      val sigma = 0.001
      val Phi = new Gaussian(0,sigma*sigma).noise
      val y = lambda.map( l => Util.fracpart(r0/l.toDouble + Phi))
      val rhats = e.pseudoranges(y)
      val rhat1 = e1.estimateRange(y.slice(0,3))
      val rhat2 = e2.estimateRange(y.slice(3,6))
      assertTrue( (rhats(0) - rhat1).abs < 1e-7 )
      assertTrue( (rhats(1) - rhat2).abs < 1e-7 )
    } 
    {
      val lambda = List(Rational(1,5),Rational(1,3),Rational(1,1),Rational(1,2), Rational(2,3), Rational(3,7))
      val subsets = List(List(0,3), List(1,2,4,5)) 
      val e = new TwoStageCRTRangeEstimator(lambda, subsets)
      val e1 = new CRTRangeEstimator(List(lambda(0),lambda(3)))
      val e2 = new CRTRangeEstimator(List(lambda(1),lambda(2),lambda(4),lambda(5)))
      val P = Util.lcm(lambda).toDouble
      //val r0 = Uniform.constructFromMinMax(0,P).noise
      val r0 = P/2.0
      val sigma = 0.001
      val Phi = new Gaussian(0,sigma*sigma).noise
      val y = lambda.map( l => Util.fracpart(r0/l.toDouble + Phi))
      val rhats = e.pseudoranges(y)
      val rhat1 = e1.estimateRange( List(0,3).map( i => y(i)) )
      val rhat2 = e2.estimateRange( List(1,2,4,5).map( i => y(i)) )
      assertTrue( (rhats(0) - rhat1).abs < 1e-7 )
      assertTrue( (rhats(1) - rhat2).abs < 1e-7 )
    } 
    {
      val lambda = List(Rational(2310,877),Rational(2310,523),Rational(2310,277),Rational(2310,221),Rational(2310,211))
      val s1 = List(0,1)
      val s2 = List(2,3,4)
      val subsets = List(s1, s2) 
      val e = new TwoStageCRTRangeEstimator(lambda, subsets)
      val e1 = new CRTRangeEstimator( s1.map( i => lambda(i) ) )
      val e2 = new CRTRangeEstimator( s2.map( i => lambda(i) ) )
      val P = Util.lcm(lambda).toDouble
      //val r0 = Uniform.constructFromMinMax(0,P).noise
      val r0 = P/2.0
      val sigma = 0.001
      val Phi = new Gaussian(0,sigma*sigma).noise
      val y = lambda.map( l => Util.fracpart(r0/l.toDouble + Phi))
      val rhats = e.pseudoranges(y)
      val rhat1 = e1.estimateRange( s1.map( i => y(i)) )
      val rhat2 = e2.estimateRange( s2.map( i => y(i)) )
      assertTrue( (rhats(0) - rhat1).abs < 1e-7 )
      assertTrue( (rhats(1) - rhat2).abs < 1e-7 )
    } 
    {
      val lambda = List(Rational(2310,877),Rational(2310,523),Rational(2310,277),Rational(2310,221),Rational(2310,211))
      val s1 = List(0,1,3)
      val s2 = List(2,4)
      val subsets = List(s1, s2) 
      val e = new TwoStageCRTRangeEstimator(lambda, subsets)
      val e1 = new CRTRangeEstimator( s1.map( i => lambda(i) ) )
      val e2 = new CRTRangeEstimator( s2.map( i => lambda(i) ) )
      val P = Util.lcm(lambda).toDouble
      //val r0 = Uniform.constructFromMinMax(0,P).noise
      val r0 = P/2.0
      val sigma = 0.001
      val Phi = new Gaussian(0,sigma*sigma).noise
      val y = lambda.map( l => Util.fracpart(r0/l.toDouble + Phi))
      val rhats = e.pseudoranges(y)
      val rhat1 = e1.estimateRange( s1.map( i => y(i)) )
      val rhat2 = e2.estimateRange( s2.map( i => y(i)) )
      assertTrue( (rhats(0) - rhat1).abs < 1e-7 )
      assertTrue( (rhats(1) - rhat2).abs < 1e-7 )
    } 
  }

  @Test
  def testTwoStageCRTNoiseLess = {
    println("Test two stage CRT estimate with no noise");
    {
      val lambda = List(Rational(1,5),Rational(1,3),Rational(1,1),Rational(1,2), Rational(2,3), Rational(3,7))
      val subsets = List(List(0,1,2), List(3,4,5)) 
      val e = new TwoStageCRTRangeEstimator(lambda, subsets)
      val P = Util.lcm(lambda).toDouble
      val r0 = P/2.0
      val y = lambda.map( l => Util.fracpart(r0/l.toDouble) )
      println(y)
      val rhat = e.estimateRange(y)
      println(r0,rhat)
      assertTrue( (Util.fracpart((r0 - rhat)/P)*P).abs < tol  )
    }
    {
      val lambda = List(Rational(2),Rational(3),Rational(5),Rational(7))
      val subsets = List(List(0,1), List(2,3)) 
      val e = new TwoStageCRTRangeEstimator(lambda, subsets)
      val P = Util.lcm(lambda).toDouble
      val r0 = 20*scala.math.Pi
      val y = lambda.map( l => Util.fracpart(r0/l.toDouble) )
      val rhat = e.estimateRange(y)
      println(r0,rhat)
      assertTrue( (Util.fracpart((r0 - rhat)/P)*P).abs < tol  )
    }
    {
      val lambda = List(Rational(2),Rational(42,13),Rational(210,37),Rational(7))
      val subsets = List(List(0,1), List(2,3)) 
      val e = new TwoStageCRTRangeEstimator(lambda, subsets)
      val P = Util.lcm(lambda).toDouble
      val r0 = 20*scala.math.Pi
      val y = lambda.map( l => Util.fracpart(r0/l.toDouble) )
      val rhat = e.estimateRange(y)
      println(r0,rhat)
      assertTrue( (Util.fracpart((r0 - rhat)/P)*P).abs < tol  )
    }
    {
      val lambda = List(Rational(210,79),Rational(210,61),Rational(210,41),Rational(210,31))
      val subsets = List(List(0,1), List(2,3)) 
      val e = new TwoStageCRTRangeEstimator(lambda, subsets)
      val P = Util.lcm(lambda).toDouble
      val r0 = 20*scala.math.Pi
      val y = lambda.map( l => Util.fracpart(r0/l.toDouble) )
      val rhat = e.estimateRange(y)
      println(r0,rhat)
      assertTrue( (Util.fracpart((r0 - rhat)/P)*P).abs < tol  )
    }
    {
      val lambda = List(Rational(210,79),Rational(210,61),Rational(210,41),Rational(210,31))
      val subsets = List(List(0,1), List(2,3)) 
      val e = new TwoStageCRTRangeEstimator(lambda, subsets)
      val P = Util.lcm(lambda).toDouble
      val r0 = 20*scala.math.Pi
      val y = lambda.map( l => Util.fracpart(r0/l.toDouble) )
      val rhat = e.estimateRange(y)
      println(r0,rhat)
      assertTrue( (Util.fracpart((r0 - rhat)/P)*P).abs < tol  )
    }
  }
  
  @Test
  def testTwoStageCRTisOrderInvariant = {
    println("Test order invariance of multistage CRT");
      val lambda1 = List(Rational(2310,877),Rational(2310,523),Rational(2310,277),Rational(2310,221),Rational(2310,211))
      val subsets1 = List(List(0,1,2), List(3,4)) 
      val e1 = new TwoStageCRTRangeEstimator(lambda1, subsets1)
      
      val lambda2 = List(Rational(2310,523),Rational(2310,277),Rational(2310,221),Rational(2310,877),Rational(2310,211))
      val subsets2 = List(List(3,0,1), List(2,4)) 
      val e2 = new TwoStageCRTRangeEstimator(lambda2, subsets2)
      
      val r0 = 20*scala.math.Pi
      
      val y1 = lambda1.map( l => Util.fracpart(r0/l.toDouble) )
      val rhat1 = e1.estimateRange(y1)
      val y2 = lambda2.map( l => Util.fracpart(r0/l.toDouble) )
      val rhat2 = e2.estimateRange(y2)
    
      println(r0,rhat1,rhat2)
      assertTrue( (rhat1 - rhat2).abs < tol  )
  }
  
  
}
