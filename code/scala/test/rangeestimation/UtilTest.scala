package rangeestimation

import org.junit._
import Assert._
import numbers.Rational
import numbers.Integer
import util.Random

class UtilTest {

  @Test
  def testRationalLCM{
    println("Test rational least common multiple")
    assertTrue( Util.lcm(Rational(3,1), Rational(4,1))==Rational(12,1) )
    assertTrue( Util.lcm(Rational(3,2), Rational(2,1))==Rational(6,1) )
    assertTrue( Util.lcm(Rational(1,5), Rational(2,5))==Rational(2,5) )
    val l1 = List(Rational(1,5), Rational(2,5), Rational(3,10))
    assertTrue( Util.lcm(l1)==Rational(6,5) )
    
    //check that pure integers result from dividing the lcm
    val l = Util.lcm(l1)
    for( v <- l1 ) assertTrue( (l/v).d == Integer(1) )
    
    //check that pure integers result from dividing the lcm for a list of 20 random Rationals
    val N = 20
    val l2 = (1 to N).map( i => Rational( Random.nextInt, Random.nextInt ) )
    val p = Util.lcm(l2)
    for( v <- l2 ) assertTrue( (p/v).d == Integer(1) ) 
  }
  
  @Test
  def testInverseMod{
    println("Test modulo inverse");
    {
      val a = Integer(3); val b = Integer(4);
      val c = Util.inverse_mod(a,b)
      //println(a, b, c, a*c, (a*c) mod b)
      assertTrue(((a*c) mod b) == Integer.one)
    }
    {
      val a = Integer(32); val b = Integer(23);
      val c = Util.inverse_mod(a,b)
      //println(a, b, c, a*c, (a*c) mod b)
      assertTrue(((a*c) mod b) == Integer.one)
    }
    {
      val a = -Integer(41); val b = Integer(446);
      val c = Util.inverse_mod(a,b)
      //println(a, b, c, a*c, (a*c) mod b)
      assertTrue(((a*c) mod b) == Integer.one)
    }
    {
      val a = Integer(24605); val b = Integer(51);
      val c = Util.inverse_mod(a,b)
      //println(a, b, c, a*c, (a*c) mod b)
      assertTrue(((a*c) mod b) == Integer.one)
    }
  }
  
  @Test
  def testChineseRemainderTheorem{
    println("Test Chinese Remainder theorem");
    {
      //testing with pairwise relatively prime moduls m
      val a = List(Integer(1),Integer(2),Integer(3))
      val m = List(Integer(3),Integer(7),Integer(5))
      val x = Util.chinese_remainder(a,m)
      for( i <- a.indices ) assertTrue( (x mod m(i)) == (a(i) mod m(i)) )
    }
  }
}
