package statespaces;

import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;

/**
 * Apply the log-normal state space for the LIBOR market model.
 *
 * @author Group 05
 */

public class LognormalStateSpace implements StateSpace {
	
	
	@Override
	public RandomVariable applyStateSpaceTransform(final MonteCarloProcess process, final int timeIndex, final int componentIndex,
			final RandomVariable randomVariable) {
		return randomVariable.exp();// exponential of the state-space
	}
	
	
	@Override
	public RandomVariable applyStateSpaceTransformInverse(final MonteCarloProcess process, final int timeIndex, final int componentIndex, final RandomVariable randomVariable) {
		return randomVariable.log();// logarithm of the state-space
	}

	@Override
	public RandomVariable[] getDriftAdjustment(final MonteCarloProcess process, final int timeIndex,
			RandomVariable[] realizationAtTimeIndex, RandomVariable[] realizationPredictor) {

		final LIBORMarketModel model = (LIBORMarketModel) process.getModel();
		final LIBORCovarianceModel covarianceModel = model.getCovarianceModel();

		final double	time					= process.getTime(timeIndex);				// t - current simulation time
		int				firstForwardRateIndex	= model.getLiborPeriodIndex(time) + 1;		// m(t)+1 - the end of the current period
		if(firstForwardRateIndex < 0) {
			firstForwardRateIndex = -firstForwardRateIndex-1 + 1;
		}

		// initialize the drift adjustment vector
		final RandomVariable[]	driftAdjustment = new RandomVariable[process.getNumberOfComponents()];
		for(int componentIndex = 0; componentIndex < process.getNumberOfComponents(); componentIndex++) {
			driftAdjustment[componentIndex] = Scalar.of(0.0);
		}

		// Drift adjustment for log-coordinate in each component
		for(int componentIndex = firstForwardRateIndex; componentIndex < process.getNumberOfComponents(); componentIndex++) {
			final RandomVariable variance = covarianceModel.getCovariance(time, componentIndex, componentIndex, realizationAtTimeIndex);
			driftAdjustment[componentIndex] = driftAdjustment[componentIndex].addProduct(variance, -0.5); // mu - sigma^2/2
		}

		return driftAdjustment;
	}

	@Override
	public RandomVariable getOneStepMeasureTransformAdjustment(final MonteCarloProcess process, final int timeIndex,
			final int componentIndex, final RandomVariable[] realizationAtTimeIndex) {
        return realizationAtTimeIndex[componentIndex];
	}
}
