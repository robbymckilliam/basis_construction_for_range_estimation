/**
  * Run simulations of various range estimators.
  */
import numbers.Rational
import org.mckilliam.distributions.Gaussian
import org.mckilliam.distributions.Uniform
import rangeestimation.Util.fracpart
import rangeestimation.Util.lcm
import rangeestimation.RangeEstimator
import rangeestimation.LeastSquaresRangeEstimator
import rangeestimation.BabaiRangeEstimator
import rangeestimation.CRTRangeEstimator
import rangeestimation.TwoStageCRTRangeEstimator

val iters = 1e7.toInt
val r0 = 20*scala.math.Pi

//generates a geometric sequence min, min*step, min*step^2, ... up to first elemet less than max
def geomseq(min : Double, max : Double, step : Double) : List[Double] = {
  def f(s : Double, list : List[Double]) : List[Double] = if(s < max) f(s*step,s :: list) else list
  return f(min,List[Double]())
}

def timenow = (new java.util.Date).getTime 
val starttime = timenow

{
  val varseq = geomseq(5e-7,2e-2,1.25) //sequence of noise variances

  val A = List(Rational(2),Rational(3),Rational(5),Rational(7)) //list of wavelenghts
  runsim("LeastSquaresA", () => new LeastSquaresRangeEstimator(A), A, varseq)
  //runsim("BabaiA", () => new BabaiRangeEstimator(A), A, varseq)
  runsim("CRTA", () => new CRTRangeEstimator(A), A, varseq)

  val B = List(Rational(210,79),Rational(210,61),Rational(210,41),Rational(210,31)) //list of modified wavelenghts
  runsim("LeastSquaresB", () => new LeastSquaresRangeEstimator(B), B, varseq)
  //runsim("BabaiB", () => new BabaiRangeEstimator(B), B, varseq)
  runsim("CRTB", () => new CRTRangeEstimator(B), B, varseq)
  runsim("TwoStageCRTB", () => new TwoStageCRTRangeEstimator(B, List(List(0,2), List(1,3))), B, varseq)
}

{
  val varseq = geomseq(7e-9,2e-2,1.35) //sequence of noise variances

  val C = List(Rational(2),Rational(3),Rational(5),Rational(7),Rational(11)) //list of wavelenghts
  runsim("LeastSquaresC", () => new LeastSquaresRangeEstimator(C), C, varseq)
  //runsim("BabaiC", () => new BabaiRangeEstimator(C), C, varseq)
  runsim("CRTC", () => new CRTRangeEstimator(C), C, varseq)

  val D = List(Rational(2310,877),Rational(2310,523),Rational(2310,277),Rational(2310,221),Rational(2310,211)) //list of modified wavelenghts
  runsim("LeastSquaresD", () => new LeastSquaresRangeEstimator(D), D, varseq)
  //runsim("BabaiD", () => new BabaiRangeEstimator(D), D, varseq)
  runsim("CRTD", () => new CRTRangeEstimator(D), D, varseq)
  runsim("TwoStageCRTD", () => new TwoStageCRTRangeEstimator(D, List(List(0,1,2),List(3,4))), D, varseq)
}


val runtime = timenow - starttime
println("Simulation finshed in " + (runtime/1000.0) + " seconds.\n")

/** Runs a simulation with given parameters and stores output in a file */
def runsim(name : String, estimator_generator : () => RangeEstimator, lambda : List[Rational], vars : List[Double]) {

  print("Running " + name + " ")
  val eststarttime = timenow
  
  val mses = vars.par.map { variance => 
    val est = estimator_generator()
    val Phi = new Gaussian(0,variance)
    Phi.setSeed((variance / 6e-5 * 20).toLong)
    val mse = (1 to iters).foldLeft(0.0) { (sum, i) => 
      val y = lambda.map( l => fracpart(r0/l.toDouble + Phi.noise))
      val rhat = est.estimateRange(y)
      sum + (r0 - rhat)*(r0 - rhat)/iters
    }
    print(".")
    mse
  }.toList

  val file = new java.io.FileWriter("data/" + name) //write output to file
  (vars, mses).zipped.foreach((v,m) => file.write(v + "\t" + m + "\n"))
  file.close

  val estruntime = timenow - eststarttime
  println(" finished in " + (estruntime/1000.0) + " seconds.")

}



