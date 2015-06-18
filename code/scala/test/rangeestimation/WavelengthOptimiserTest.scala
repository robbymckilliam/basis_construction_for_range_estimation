package rangeestimation

import org.junit._
import Assert._
import numbers.Rational

class WavelengthOptimiserTest {
  
  @Test
  def testLandL1equivalence{
    {
      val w = List[Rational](Rational(1,2),Rational(1,3),Rational(3,5))
      val N = w.size+1
      val lambdamax = Rational(4,5)
      val lambdamin = Rational(1,5)
      val wopt = new WavelengthOptimiser(N,1,lambdamin,lambdamax,1)
      assertTrue(wopt.L(lambdamax :: w) == wopt.L1(w))
    }
    {
      val w = List[Rational](Rational(1,2),Rational(1,3),Rational(4,5),Rational(2,7),Rational(2,15))
      val N = w.size+1
      val lambdamax = Rational(5,4)
      val lambdamin = Rational(1,5)
      val wopt = new WavelengthOptimiser(N,1,lambdamin,lambdamax,1)
      assertTrue(wopt.L(lambdamax :: w) == wopt.L1(w))
    }
  }
  
  @Test
  def testL1andL2equivalence{
    {
      val w = List[Rational](Rational(1,2),Rational(1,3),Rational(3,5))
      val N = w.size+1
      val lambdamax = Rational(4,5)
      val lambdamin = Rational(1,5)
      val p = (0 until N-1).map( n => (w(n)/lambdamax).n )
      val q = (0 until N-1).map( n => (w(n)/lambdamax).d )
      val wopt = new WavelengthOptimiser(N,1,lambdamin,lambdamax,1)
      for( n <- 0 until N-1 ) assertTrue(w(n) == (lambdamax*p(n))/q(n))    
      assertTrue(wopt.L1(w) == wopt.L2(p,q))
    }

    {
      val w = List[Rational](Rational(1,2),Rational(1,3),Rational(4,5),Rational(2,7),Rational(2,15))
      val N = w.size+1
      val lambdamax = Rational(5,4)
      val lambdamin = Rational(1,5)
      val p = (0 until N-1).map( n => (w(n)/lambdamax).n )
      val q = (0 until N-1).map( n => (w(n)/lambdamax).d )
      val wopt = new WavelengthOptimiser(N,1,lambdamin,lambdamax,1)
      for( n <- 0 until N-1 ) assertTrue(w(n) == (lambdamax*p(n))/q(n))    
      assertTrue(wopt.L1(w) == wopt.L2(p,q))
    }
    
  }
  
  @Test
  def testPandQequivalence{
    {
      val w = List[Rational](Rational(1),Rational(1,3),Rational(3,5))
      val N = w.size
      val p = (1 until N).map( n => w(n).n )
      assertTrue(Util.lcm(w).isInteger)
      assertTrue(Util.lcm(w) == Rational(numbers.EuclideanDomain.lcm(p)))
    }
    {
      val w = List[Rational](Rational(1,2),Rational(1,3),Rational(4,5),Rational(2,7),Rational(2,15))
      val N = w.size
      val p = (1 until N).map( n => w(n).n )
      assertTrue(Util.lcm(w).isInteger)
      assertTrue(Util.lcm(w) == Rational(numbers.EuclideanDomain.lcm(p)))
    }
  }
  
  @Test 
  def testd {
    val w = List[Rational](Rational(1,2),Rational(1,3),Rational(3,5))
    val N = w.size+1
    val lambdamax = Rational(4,5)
    val lambdamin = Rational(1,5)
    val rmin = 5
    val wopt = new WavelengthOptimiser(N,rmin,lambdamin,lambdamax,1)
    val d = wopt.d
    assertTrue(lambdamax*d/(d+1) <= lambdamax)
    assertTrue(lambdamax*d/(d+1) >= lambdamin)
    assertTrue( Util.lcm(lambdamax, lambdamax*d/(d+1)) >= rmin)
  }
  
  @Test 
  def testOptimise {
    val N = 2
    val lambdamax = Rational(7)
    val lambdamin = Rational(2)
    val rmin = 210
    val beta = 10000
    val wopt = new WavelengthOptimiser(N,rmin,lambdamin,lambdamax,beta)
    assertTrue(wopt.wavelengths(0) == lambdamax)
    println(wopt.wavelengths)
  }

}
