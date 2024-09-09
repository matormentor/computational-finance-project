package statespaces;

import java.util.Arrays;

import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;

/**
 * Apply the displaced log-normal state space for the LIBOR market model with the displacements
 *
 * @author Group 05
 */

public class DisplacedLognormalStateSpace extends ParametricStateSpace {

	private final double[] displacements; // d_0, ..., d_{n-1}
	private final boolean isDisplacementConstant;

	public DisplacedLognormalStateSpace(final double[] displacements, final boolean isDisplacementConstant) {
		
		this.displacements = displacements;
		this.isDisplacementConstant = isDisplacementConstant;
		
	}

	public DisplacedLognormalStateSpace(final double[] displacements) {
		this(displacements, false);
	}

	@Override
	public double[] getParameterAsDouble() {
		if (isDisplacementConstant) {
			return new double[] {displacements[0]}; //If constant displacement, just use the first element
		} else {
			double[] parameters = new double[displacements.length];
			System.arraycopy(displacements, 0, parameters, 0, parameters.length); // Copy displacements to parameters
			return parameters;
		}
	}

	@Override
	public ParametricStateSpace getCloneWithModifiedParameters(final double[] parameters) {
		if(parameters == null || parameters.length == 0) {
			return this; 
		}

		final double [] displacementParameters = new double[displacements.length];
		if (!isDisplacementConstant) {
			System.arraycopy(parameters, 0, displacementParameters, 0, displacementParameters.length); //pass parameters to displacement parameters
		} else {
			Arrays.fill(displacementParameters, parameters[0]); 
		}

		return new DisplacedLognormalStateSpace(displacementParameters, isDisplacementConstant);
	}

	@Override
	public RandomVariable applyStateSpaceTransform(
			final MonteCarloProcess process, 
			final int timeIndex, 
			final int componentIndex,
			final RandomVariable randomVariable) {
		return randomVariable.exp().sub(displacements[componentIndex]); // Y = exp(X)-d
	}

	@Override
	public RandomVariable applyStateSpaceTransformInverse(
			final MonteCarloProcess process, 
			final int timeIndex, 
			final int componentIndex, 
			final RandomVariable randomVariable) {
		return randomVariable.add(displacements[componentIndex]).log();// log(Y+d)
	}

	@Override
	public RandomVariable[] getDriftAdjustment(
			final MonteCarloProcess process,
			final int timeIndex,
			final RandomVariable[] realizationAtTimeIndex, 
			final RandomVariable[] realizationPredictor) {
		final LIBORMarketModel model = (LIBORMarketModel) process.getModel();
		final LIBORCovarianceModel covarianceModel = model.getCovarianceModel();

		final double	time					= process.getTime(timeIndex);				// t - current simulation time
		int				firstForwardRateIndex	= model.getLiborPeriodIndex(time)+1;		// m(t)+1 - the end of the current period
		if(firstForwardRateIndex < 0) {
			firstForwardRateIndex = -firstForwardRateIndex-1 + 1;
		}

		// Initialize the drift adjustment vector
		final RandomVariable[]	driftAdjustment = new RandomVariable[process.getNumberOfComponents()]; // initialize drift adjustment
		for(int componentIndex = 0; componentIndex< process.getNumberOfComponents(); componentIndex++) {
			driftAdjustment[componentIndex] = Scalar.of(0.0);
		}

		// Drift adjustment for log-coordinate in each component
		for(int componentIndex = firstForwardRateIndex; componentIndex < process.getNumberOfComponents(); componentIndex++) {
			final RandomVariable		variance		= covarianceModel.getCovariance(time, componentIndex, componentIndex, realizationAtTimeIndex);
			driftAdjustment[componentIndex] = driftAdjustment[componentIndex].addProduct(variance, -0.5); // mu - sigma^2/2
		}

		return driftAdjustment;
	}

	@Override
	public RandomVariable getOneStepMeasureTransformAdjustment(MonteCarloProcess process, int timeIndex,
			int componentIndex, RandomVariable[] realizationAtTimeIndex) {
        // L[i]+d
        return realizationAtTimeIndex[componentIndex].add(displacements[componentIndex]);
	}

	public double[] getDisplacements() {
		return displacements;
	}
}
