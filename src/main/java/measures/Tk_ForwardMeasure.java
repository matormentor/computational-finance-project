package measures;


import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Arrays;
import java.util.Locale;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;
import net.finmath.time.TimeDiscretization;
import statespaces.StateSpace;

/**
 * Implements the Tk-forward measure, where the numéraire is the zero coupon bond with maturity Tk and after maturity
 * a rolling bond into the next single period.
 *
 * @author Group 05
 */

public class Tk_ForwardMeasure implements Measure {

	private final int maturityIndex_k;	// The maturity index of the maturity T_k

	/**
	 * @param maturityIndex_k The maturity index of the maturity T_k
	 */
	public Tk_ForwardMeasure(int maturityIndex_k) {
		this.maturityIndex_k = maturityIndex_k;
	}

	@Override
	public RandomVariable getNumeraireUnadjusted(MonteCarloProcess process, double time) throws CalculationException {
		final LIBORMarketModel model = (LIBORMarketModel) process.getModel();
		final TimeDiscretization liborPeriodDiscretization = model.getLiborPeriodDiscretization();

		/*
		 * Check if numéraire is on LIBOR time grid
		 */
		
		final int liborTimeIndex = model.getLiborPeriodIndex(time);
		RandomVariable numeraireUnadjusted = model.getRandomVariableForConstant(1.0);	// Initialize to 1.0;

		if (liborTimeIndex < 0) {
			final int upperIndex = -liborTimeIndex - 1;
			final int lowerIndex = upperIndex - 1;
			if (lowerIndex < 0) {
				throw new IllegalArgumentException("Numeraire requested for time " + time + "is unsupported");
			}

			/*
			 * Depending on the measure, we calculate the corresponding numéraire.
			 */
			
			if (time <= liborPeriodDiscretization.getTime(maturityIndex_k)) {
				for (int liborIndex = upperIndex; liborIndex < maturityIndex_k; liborIndex++) {
					final RandomVariable libor = process.getProcessValue(
							process.getTimeIndex(Math.min(time, liborPeriodDiscretization.getTime(liborIndex))), liborIndex); // Get LIBOR
					final double periodLength = liborPeriodDiscretization.getTimeStep(liborIndex);
					numeraireUnadjusted = numeraireUnadjusted.discount(libor, periodLength);
				}
			} else {
				numeraireUnadjusted = getNumeraireUnadjusted(process, model.getLiborPeriod(upperIndex));
			}

			/*
			 * Multiply with short period bond
			 */
			numeraireUnadjusted = numeraireUnadjusted.discount(
					model.getForwardRate(process, time, time, model.getLiborPeriod(upperIndex)), model.getLiborPeriod(upperIndex) - time);

			return numeraireUnadjusted;
		}
		else {
			/*
			 * Calculate the numéraire, when time is part of liborPeriodDiscretization
			 */
			return getNumeraireUnadjustedAtLIBORIndex(process, liborTimeIndex);
		}
	}

	private RandomVariable getNumeraireUnadjustedAtLIBORIndex(MonteCarloProcess process, int liborTimeIndex) throws CalculationException {
		final LIBORMarketModel model = (LIBORMarketModel) process.getModel();
		final TimeDiscretization liborPeriodDiscretization = model.getLiborPeriodDiscretization();

		RandomVariable numeraireUnadjusted = model.getRandomVariableForConstant(1.0);

		/*
		 * Depending on the measure, we calculate the corresponding (unadjusted) numéraire.
		 */
		if (liborTimeIndex < maturityIndex_k) {
			int timeIndex = process.getTimeIndex(liborPeriodDiscretization.getTime(liborTimeIndex));
			if(timeIndex < 0) {
				timeIndex = -timeIndex - 1;
			}

			for (int liborIndex = liborTimeIndex; liborIndex < maturityIndex_k; liborIndex++) {
				final RandomVariable libor = process.getProcessValue(timeIndex, liborIndex);
				final double periodLength = liborPeriodDiscretization.getTimeStep(liborIndex);
				numeraireUnadjusted = numeraireUnadjusted.discount(libor, periodLength);
			}
		} else {
			// If numeraire is not N(T_k), multiply (1 + L_{i-1}(T_{i-1})*(T_i-T_{i-1}) on N(T_{i-1})
			if (liborTimeIndex != maturityIndex_k) {
				int timeIndex = process.getTimeIndex(liborPeriodDiscretization.getTime(liborTimeIndex-1));
				if(timeIndex < 0) {
					timeIndex = -timeIndex - 1;
				}

				final RandomVariable libor = process.getProcessValue(timeIndex, liborTimeIndex - 1);
				final double periodLength = liborPeriodDiscretization.getTimeStep(liborTimeIndex - 1);
				numeraireUnadjusted = getNumeraireUnadjustedAtLIBORIndex(process, liborTimeIndex - 1).accrue(libor, periodLength);
			}
		}

		return numeraireUnadjusted;
	}

	@Override
	public RandomVariable[] getDrift(final MonteCarloProcess process, final int timeIndex, final RandomVariable[] realizationAtTimeIndex,
			final RandomVariable[] realizationPredictor, final StateSpace stateSpace) {
		final LIBORMarketModel model = (LIBORMarketModel) process.getModel();

		final double	time					= process.getTime(timeIndex);			// t - current simulation time
		int				firstForwardRateIndex	= model.getLiborPeriodIndex(time) + 1;		// m(t)+1 - the end of the current period
		if(firstForwardRateIndex < 0) {
			firstForwardRateIndex = -firstForwardRateIndex-1 + 1;
		}

		final RandomVariable		zero	= Scalar.of(0.0);
		
		
		// Step 1: Allocate and initialize the drift vector 
		final RandomVariable[]	drift = new RandomVariable[process.getNumberOfComponents()];
		for(int componentIndex = 0; componentIndex < process.getNumberOfComponents(); componentIndex++) {
			drift[componentIndex] = zero;
		}

		// Step 2:  Allocate array (for each k) for the sums of delta_{i}/(1+L_{i} \delta_i) f_{i,k}
		final RandomVariable[]	factorLoadingsSums	= new RandomVariable[process.getNumberOfFactors()];
		Arrays.fill(factorLoadingsSums, zero);

		// Step 3: Calculate drift for the component componentIndex
		if (time <= model.getLiborPeriod(maturityIndex_k)) { // Case 1: for t smaller than or equal T_k

			// first, backwards starting at k-1 till firstForwardRateIndex
			for(int componentIndex = maturityIndex_k - 1; componentIndex >= firstForwardRateIndex; componentIndex--) {
				final RandomVariable	forwardRate		= realizationAtTimeIndex[componentIndex];
				final double			periodLength	= model.getLiborPeriodDiscretization().getTimeStep(componentIndex);
				RandomVariable			oneStepMeasureTransform = Scalar.of(-periodLength).discount(forwardRate, periodLength);

				oneStepMeasureTransform = oneStepMeasureTransform.mult(stateSpace.getOneStepMeasureTransformAdjustment(process, timeIndex, componentIndex, realizationAtTimeIndex));

				final RandomVariable[]	factorLoading   	= model.getFactorLoading(process, timeIndex, componentIndex, realizationAtTimeIndex);
				drift[componentIndex] = drift[componentIndex].addSumProduct(factorLoadingsSums, factorLoading);
				for(int factorIndex = 0; factorIndex < factorLoading.length; factorIndex++) {
					factorLoadingsSums[factorIndex] = factorLoadingsSums[factorIndex].addProduct(oneStepMeasureTransform, factorLoading[factorIndex]);
				}
			}

			// next forwards starting at maturityIndex_k, others are zero
			Arrays.fill(factorLoadingsSums, zero);
			for(int componentIndex = maturityIndex_k; componentIndex < process.getNumberOfComponents(); componentIndex++) { // forward from k, ..., n-1
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
		} else { // Case 2: for t greater than T_k
			for(int componentIndex = firstForwardRateIndex; componentIndex < process.getNumberOfComponents(); componentIndex++) { // forward from m(t) + 1, ..., n-1
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
		}
		// Step 4:  Drift adjustment in each component depending on the state space
		final RandomVariable[] driftAdjustments = stateSpace.getDriftAdjustment(process, timeIndex, realizationAtTimeIndex, realizationPredictor);

		for(int componentIndex = firstForwardRateIndex; componentIndex < process.getNumberOfComponents(); componentIndex++) {
			drift[componentIndex] = drift[componentIndex].add(driftAdjustments[componentIndex]);
		}

		return drift;
	}

	@Override
	public String toString() {
		final DecimalFormat formatterValue	= new DecimalFormat("00", new DecimalFormatSymbols(Locale.ENGLISH));

		return "T_k-Forward Measure; [k = " + formatterValue.format(maturityIndex_k) + "]";
	}

	/**
	 * @return the maturityIndexTk
	 */
	public int getMaturityIndexTk() {
		return maturityIndex_k;
	}
}

