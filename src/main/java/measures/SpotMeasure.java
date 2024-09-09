package measures;


import java.util.Arrays;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;
import net.finmath.time.TimeDiscretization;
import statespaces.StateSpace;

/**
 * Implements the spot measure with the single period rolling bond N as numeraire.
 *
 * @author Group 05
 */
public class SpotMeasure implements Measure {

	@Override
	public RandomVariable getNumeraireUnadjusted(final MonteCarloProcess process, final double time) throws CalculationException {
		final LIBORMarketModel model = (LIBORMarketModel) process.getModel();

		final int liborTimeIndex = model.getLiborPeriodIndex(time);
        model.getRandomVariableForConstant(1.0);
        RandomVariable numeraireUnadjusted;	// Initialize to 1.0

		if (liborTimeIndex < 0) {
			final int upperIndex = -liborTimeIndex - 1;
			final int lowerIndex = upperIndex - 1;
			if (lowerIndex < 0) {
				throw new IllegalArgumentException("Numeraire requested for time " + time + ". Unsupported");
			}

			// Depending on the measure, we calculate the corresponding (unadjusted) numeraire.
			numeraireUnadjusted = getNumeraireUnadjusted(process, model.getLiborPeriod(upperIndex));

			//Multiply with short period bond
			numeraireUnadjusted = numeraireUnadjusted.discount(
					model.getForwardRate(process, time, time, model.getLiborPeriod(upperIndex)), model.getLiborPeriod(upperIndex) - time);

			return numeraireUnadjusted;
		}
		else {
			/*
			 * Calculate the numeraire, when time is part of liborPeriodDiscretization
			 */
			return getNumeraireUnadjustedAtLIBORIndex(process, liborTimeIndex);
		}
	}

	private RandomVariable getNumeraireUnadjustedAtLIBORIndex(final MonteCarloProcess process, final int liborTimeIndex) throws CalculationException {
		final LIBORMarketModel model = (LIBORMarketModel) process.getModel();
		final TimeDiscretization liborPeriodDiscretization = model.getLiborPeriodDiscretization();

		/*
		 * Depending on the measure, we calculate the corresponding (unadjusted) numeraire.
		 * If numeraire is not N(0), multiply (1 + L_{i-1}(T_{i-1})*(T_i-T_{i-1}) on N(T_{i-1})
		 */
		
		RandomVariable numeraireUnadjusted = model.getRandomVariableForConstant(1.0);
		if (liborTimeIndex != 0) {
			int timeIndex = process.getTimeIndex(liborPeriodDiscretization.getTime(liborTimeIndex-1));
			if(timeIndex < 0) {
				timeIndex = -timeIndex - 1;
			}

			final double periodLength = liborPeriodDiscretization.getTimeStep(liborTimeIndex - 1);
			final RandomVariable libor = process.getProcessValue(timeIndex, liborTimeIndex - 1);
			numeraireUnadjusted = getNumeraireUnadjustedAtLIBORIndex(process, liborTimeIndex - 1).accrue(libor, periodLength);
		}

		return numeraireUnadjusted;
	}

	@Override
	public RandomVariable[] getDrift(final MonteCarloProcess process, final int timeIndex, final RandomVariable[] realizationAtTimeIndex,
			final RandomVariable[] realizationPredictor, final StateSpace stateSpace) {
		final LIBORMarketModel model = (LIBORMarketModel) process.getModel();

		final double	time					= process.getTime(timeIndex);				// t - current simulation time
		int				firstForwardRateIndex	= model.getLiborPeriodIndex(time) + 1;		// m(t)+1 - the end of the current period
		if(firstForwardRateIndex < 0) {
			firstForwardRateIndex = -firstForwardRateIndex-1 + 1;
		}

		final RandomVariable		zero	= Scalar.of(0.0);

		// Initialize the drift vector 
		final RandomVariable[]	drift = new RandomVariable[process.getNumberOfComponents()];
		for(int componentIndex = 0; componentIndex < process.getNumberOfComponents(); componentIndex++) {
			drift[componentIndex] = zero;
		}

		// Allocate array (for each k) for the sums of delta_{i}/(1+L_{i} \delta_i) f_{i,k} (+ for spot measure, - for terminal measure)
		final RandomVariable[]	factorLoadingsSums	= new RandomVariable[process.getNumberOfFactors()];		// alpha_k
		Arrays.fill(factorLoadingsSums, zero);
		

		// Calculate drift for the component componentIndex (starting at firstForwardRateIndex, others are zero)
		for(int componentIndex = firstForwardRateIndex; componentIndex < process.getNumberOfComponents(); componentIndex++) { // forward from m(t)+1, ..., n-1
			final RandomVariable	forwardRate		= realizationAtTimeIndex[componentIndex];
			final double			periodLength	= model.getLiborPeriodDiscretization().getTimeStep(componentIndex);
			RandomVariable			oneStepMeasureTransform = Scalar.of(periodLength).discount(forwardRate, periodLength);

			oneStepMeasureTransform = oneStepMeasureTransform.mult(
					stateSpace.getOneStepMeasureTransformAdjustment(process, timeIndex, componentIndex, realizationAtTimeIndex));

			final RandomVariable[]	factorLoading   	= model.getFactorLoading(process, timeIndex, componentIndex, realizationAtTimeIndex);
			for(int factorIndex = 0; factorIndex < factorLoading.length; factorIndex++) {
				factorLoadingsSums[factorIndex] = factorLoadingsSums[factorIndex].addProduct(oneStepMeasureTransform, factorLoading[factorIndex]);
			}
			drift[componentIndex] = drift[componentIndex].addSumProduct(factorLoadingsSums, factorLoading);
		}

		// Drift adjustment in each component depending on the state space
		final RandomVariable[] driftAdjustments = stateSpace.getDriftAdjustment(process, timeIndex, realizationAtTimeIndex, realizationPredictor);

		for(int componentIndex = firstForwardRateIndex; componentIndex < process.getNumberOfComponents(); componentIndex++) {
			drift[componentIndex] = drift[componentIndex].add(driftAdjustments[componentIndex]);
		}

		return drift;
	}

	@Override
	public String toString() {
		return "Spot Measure";
	}
}
