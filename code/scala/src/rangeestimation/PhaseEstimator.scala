/*
 * @author Robby McKilliam
 */

package rangeestimation

import numbers.finite.Complex
import numbers.finite.RectComplex
import numbers.finite.PolarComplex
import scala.math.Pi
import scala.math.sin

/** Base trait for estimators of phase from sampled sinusoids */
trait PhaseEstimator {

  def estimatePhase(y : Seq[Double]) : Double
  
}

/** Estimates the phase a sinusoidal signal of frequency f from L samples at rate T */
class FourierCoefficientEstimator(val f : Double, val L : Int, val T : Double) extends PhaseEstimator {
  
  /** The constant A resulting from a particular geometric progression */
  val A = PolarComplex(1.0,-2*Pi*f*T*(L-1))*sin(2*Pi*f*L*T)/sin(2*Pi*f*T)/L
  
  def estimatePhase(y : Seq[Double]) : Double = {
    if(y.length != L) throw new ArrayIndexOutOfBoundsException("Number of samples must be " + L)
    val a = calculate_a(y)
    val b = a  - a.conjugate * A
    return b.angle/2/Pi
  }
  
  def calculate_a(y : Seq[Double]) : Complex = {
    var sum = Complex.zero
    for( ell <- 0 to L-1) sum = sum + PolarComplex(1.0,-2*Pi*f*ell*T) * y(ell)
    return sum * RectComplex(0.0,2.0/L)
  }
  
}
