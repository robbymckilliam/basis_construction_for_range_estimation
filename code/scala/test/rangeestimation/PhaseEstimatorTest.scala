/*
 * @author Robby McKilliam.
 */

package rangeestimation

import org.junit._
import Assert._
import scala.math.sin
import scala.math.Pi
import numbers.finite.Complex
import numbers.finite.RectComplex
import numbers.finite.PolarComplex
import Util.fracpart

class PhaseEstimatorTest {
  
  val tol = 1e-6
  
  val f  = 1000.0 //frequency (Hertz)
  val phi = 0.2 //phase
  val d = 11.0 //distance (meters)
  val c = 343.0 //meters per second (something like the speed of sound)
  val T = 1.0/44100.0 //CD quality sample period
  val L : Int = 300 //number of samples
  val lambda = c/f //wavelength in meters
  val alpha = 0.4 //amplitude of recieved signal
  val theta = phi - d/lambda //recieved phase
  def x(t : Double) = sin(2*Pi*f*t + 2*Pi*phi) //transmitted signal
  def y(t: Double) = alpha*x(t - d/c) //recieved signal
  
  @Test
  def testAWithoutNoise = {
    println("Testing A without noise")
    
    val ys = (0 to L-1).map(i => y(i*T)) //construct sampled signal (no noise in this case)
    val est = new FourierCoefficientEstimator(f,L,T) //construct the phase estimator
    val a = est.calculate_a(ys) //the complex number a used in constructing a phase estimate
    val A = est.A
    val ea = PolarComplex(alpha,2*Pi*theta) - PolarComplex(alpha,-2*Pi*theta)*A
    //println( A + ", " + a + ", " + ea )
    assertTrue( (ea - a).magnitude < tol )
    
  }

  @Test
  def testAWithNoise = {
    println("Testing A with noise")
    
    val w = (0 to L-1).map(i => scala.math.random)
    val ys = (0 to L-1).map(i => y(i*T) + w(i)) //construct noisy sampled signal
    val est = new FourierCoefficientEstimator(f,L,T) //construct the phase estimator
    val a = est.calculate_a(ys) //the complex number a used in constructing a phase estimate
    val A = est.A
    val W = (0 to L-1).foldLeft(Complex.zero)( (acc, i) => acc + PolarComplex(1.0, -2*Pi*f*i*T)*w(i) ) * RectComplex(0.0,2.0/L) 
    val ea = PolarComplex(alpha,2*Pi*theta) - PolarComplex(alpha,-2*Pi*theta)*A + W
    //println( A + ", " + W + ", " + a + ", " + ea )
    assertTrue( (ea - a).magnitude < tol )
    
  }
  
}
