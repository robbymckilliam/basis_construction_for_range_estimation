/**
* Run benchmarks.
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

val Ns = 2 to 70
val variance = 0.00001
val r0 = 20*scala.math.Pi


//time used for benchmark warmup in nanoseconds.  Benchmark runtime is approximately twice this.
val MIN_BENCH_DURATION : Long = 1000000000L; // (2 secs)  
val NUM_ITERATIONS = 200

runbench(Ns, (N : Int) => new LeastSquaresRangeEstimator((1 to N).map(Rational(_))), "LeastSquares")
runbench(Ns, (N : Int) => new CRTRangeEstimator((1 to N).map(Rational(_))), "CRT")

/** Runs a simulation with given parameters and stores output in a file */
def runbench(Ns: Seq[Int], estimator_generator : Int => RangeEstimator, name : String) {

  println("Running " + name)

  val runtimeslist = Ns.map{ N =>
    print(" N = " + N)
			    
    val Phi = new Gaussian(0,variance)
    val est = estimator_generator(N)
    val lambda = (1 to N).map(Rational(_))

    print(" ... Warming up ... ")
    var numiters = 0
    val warmupstarttime = System.nanoTime
    while(System.nanoTime - warmupstarttime < MIN_BENCH_DURATION){
      val y = lambda.map( l => fracpart(r0/l.toDouble + Phi.noise) )
      val rhat = est.estimateRange(y)
      numiters = numiters+1
    }
   
    //warmup override (comment out to fix benchmark time instead of iterations
   numiters = NUM_ITERATIONS
   
   print("Benchmarking ... ")
   val benchstarttime = System.nanoTime
   for( i <- 1 to numiters) {
     val y = lambda.map( l => fracpart(r0/l.toDouble + Phi.noise) )
     val rhat = est.estimateRange(y)
   }
   val tNano = (System.nanoTime - benchstarttime +(numiters)/2) / numiters //copied from Alan Eliasen java BigInteger benchmarks
   val tMilli = tNano/1000000.0 //time in milliseconds per iteration
   println(tMilli + " ms per iteration")
   tMilli //return time into list
  
  }.toList

  println
      
  //write output to files
  val filet = new java.io.FileWriter("data/" + name + "Benchmark")
  (runtimeslist, Ns).zipped.foreach
  { (time, N) =>
    filet.write(N.toString.replace('E', 'e') + "\t" + time.toString.replace('E', 'e')  + "\n") 
  }
  filet.close
 
}
