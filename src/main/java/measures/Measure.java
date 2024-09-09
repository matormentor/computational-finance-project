package measures;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;
import statespaces.StateSpace;

/**
 * Interface for different measures as object-oriented design.
 *
 * @author Group 05
 */

public interface Measure {
	
	/**
	 * Return the numeraire at a given time index (unadjusted without control variate).
	 *
	 * @param process The discretization process generating this model. The process provides call backs for TimeDiscretization and allows calls to getProcessValue for timeIndices less or equal the given one.
	 * @param time The time t for which the numeraire N(t) should be returned.
	 * @return The numeraire at the specified time as type RandomVariable
	 */
	RandomVariable getNumeraireUnadjusted(MonteCarloProcess process, double time) throws CalculationException;

	/**
	 * Return the drift, i.e. the coefficient vector \mu =  (\mu_{1}, ..., \mu_{n})
	 * such that X = f(Y) and
	 * dY_j = \mu_j dt + \lambda_{1,j} dW_1 + ... + \lambda_{n,j} dW_n;
	 * in an m-factor model. Here j denotes the index of the component of the resulting process.
	 * <p>
	 * Since the model is provided only on a time discretization, the method may also return the drift
	 * as \frac{1}{t_{i+1}-t_{i}} \int_{t_{i}}^{t_{i+1}} \mu(s) \mathrm{d}s
	 *
	 * @param process The discretization process generating this model. The process provides call backs for TimeDiscretization and allows calls to getProcessValue for timeIndices less or equal the given one.
	 * @param timeIndex The time index (related to the model times discretization).
	 * @param realizationAtTimeIndex The given realization at timeIndex
	 * @param realizationPredictor The given realization at <code>timeIndex+1</code> or null if no predictor is available.
	 * @return The drift or average drift from timeIndex i to i+1.
	 */
	RandomVariable[] getDrift(MonteCarloProcess process, int timeIndex, RandomVariable[] realizationAtTimeIndex, 
			RandomVariable[] realizationPredictor, StateSpace stateSpace);


}
