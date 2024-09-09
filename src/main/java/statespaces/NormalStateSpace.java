package statespaces;

import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;

/**
 * Apply the normal state space for the LIBOR market model.
 *
 * @author Group 05
 */
public class NormalStateSpace implements StateSpace {

	@Override
	public RandomVariable applyStateSpaceTransform(
			final MonteCarloProcess process, 
			final int timeIndex, 
			final int componentIndex,
			final RandomVariable randomVariable) {
		return randomVariable;
	}

	@Override
	public RandomVariable applyStateSpaceTransformInverse(
			final MonteCarloProcess process, 
			final int timeIndex, 
			final int componentIndex, 
			final RandomVariable randomVariable) {
		return randomVariable;
	}

	@Override
	public RandomVariable[] getDriftAdjustment(
			final MonteCarloProcess process, 
			final int timeIndex,
			final RandomVariable[] realizationAtTimeIndex, 
			final RandomVariable[] realizationPredictor) {

		// initialize the drift adjustment vector
		final RandomVariable[]	driftAdjustment = new RandomVariable[process.getNumberOfComponents()];
		for(int componentIndex = 0; componentIndex < process.getNumberOfComponents(); componentIndex++) {
			driftAdjustment[componentIndex] = Scalar.of(0.0);
		}
		return driftAdjustment; // no adjustment here for normal
	}

	@Override
	public RandomVariable getOneStepMeasureTransformAdjustment(
			final  MonteCarloProcess process, 
			final int timeIndex,
			final int componentIndex, 
			final RandomVariable[] realizationAtTimeIndex) {
		return Scalar.of(1.0);
	}
}

