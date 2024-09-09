package statespaces;

import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;

/**
 * Interface for different state-spaces as object-oriented design.
 *
 * @author Group 05
 */

public interface StateSpace {

	/**
	 * Applies the state space transform f_{i} to the given state random variable
	 * such that Y_{i} \mapsto f_{i}(Y_{i}) := X_{i}.
	 *
	 * @param process The discretization process generating this model. The process provides call backs for TimeDiscretization and allows calls to getProcessValue for timeIndices less or equal the given one.
	 * @param timeIndex The time index (related to the model times discretization).
	 * @param componentIndex The component index i.
	 * @param randomVariable The state random variable Y_{i}.
	 * @return New random variable holding the result of the state space transformation.
	 */
	RandomVariable applyStateSpaceTransform(MonteCarloProcess process, int timeIndex, int componentIndex, RandomVariable randomVariable);

	/**
	 * Applies the inverse state space transform f_{i}^{−1} to the given random variable
	 * such that f_{i}^{-1} (\xi) := \xi.
	 *
	 * @param process The discretization process generating this model. The process provides call backs for TimeDiscretization and allows calls to getProcessValue for timeIndices less or equal the given one.
	 * @param timeIndex The time index (related to the model times discretization).
	 * @param componentIndex The component index <i>i</i>.
	 * @param randomVariable The state random variable <i>X<sub>i</sub></i>.
	 * @return New random variable holding the result of the state space transformation.
	 */
	default RandomVariable applyStateSpaceTransformInverse(MonteCarloProcess process, int timeIndex, int componentIndex, RandomVariable randomVariable) {
		throw new UnsupportedOperationException("Inverse of state space transform not set");
	}

	/**
	 * return the drift adjustment resulting from the state space transform of the underlying process.
	 *
	 * @param process The discretization process generating this model. The process provides call backs for TimeDiscretization and allows calls to getProcessValue for timeIndices less or equal the given one.
	 * @param timeIndex The time index (related to the model times discretization).
	 * @param realizationAtTimeIndex The given realization at timeIndex
	 * @param realizationPredictor The given realization at timeIndex+1 or null if no predictor is available.
	 * @return The drift or average drift adjustment (approximation) from timeIndex to timeIndex+1
	 */
	RandomVariable[] getDriftAdjustment(MonteCarloProcess process, int timeIndex, RandomVariable[] realizationAtTimeIndex, RandomVariable[] realizationPredictor);

	/**
	 * return the one-step measure transform adjustment resulting from the state space transform of the underlying process.
	 *
	 * @param process					The discretization process generating this model.
	 * @param timeIndex					The time index (related to the model times discretization).
	 * @param componentIndex			The LIBOR Index to which the adjustment to be computed
	 * @param realizationAtTimeIndex	The given realization of the discretized LIBOR process at timeIndex
	 * @return							The one-step measure transform adjustment
	 */
	RandomVariable getOneStepMeasureTransformAdjustment(MonteCarloProcess process, int timeIndex, int componentIndex, RandomVariable[] realizationAtTimeIndex);

}

