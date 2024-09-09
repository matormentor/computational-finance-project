package models;

import java.time.LocalDateTime;

import measures.Measure;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.AnalyticModel;
import net.finmath.marketdata.model.curves.DiscountCurve;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;
import statespaces.StateSpace;

/**
 * A base class and interface description for the LIBOR market model
 * with specification of the state space and of the measure.
 *
 * @author Group 05
 */
public abstract class AbstractExtendedLIBORMarketModel implements LIBORMarketModel {

	private final TimeDiscretization		liborPeriodDiscretization;

	private final AnalyticModel			curveModel;
	private final ForwardCurve			forwardRateCurve;
	private final DiscountCurve			discountCurve;

	private final RandomVariableFactory	randomVariableFactory;
	
	private final Measure 	measure;
	private final StateSpace  stateSpace;

	/**
	 * @param liborPeriodDiscretization	The discretization of the interest rate curve into forward rates (tenor structure).
	 * @param curveModel				The associated analytic model of this model (containing the associated market data objects like curve).
	 * @param forwardRateCurve			The initial values for the forward rates.
	 * @param discountCurve				The discount curve to use. This will create an LMM model with a deterministic zero-spread discounting adjustment.
	 * @param randomVariableFactory		The random variable factory used to create the initial values of the model.
	 */
	public AbstractExtendedLIBORMarketModel(final TimeDiscretization liborPeriodDiscretization, final AnalyticModel curveModel,
			final ForwardCurve forwardRateCurve, final DiscountCurve discountCurve, final RandomVariableFactory randomVariableFactory,final Measure measure,final StateSpace stateSpace ) {
		super();
		
		this.liborPeriodDiscretization = liborPeriodDiscretization;
		this.curveModel = curveModel;
		this.forwardRateCurve = forwardRateCurve;
		this.discountCurve = discountCurve;
		this.randomVariableFactory = randomVariableFactory;
		this.measure = measure;
		this.stateSpace = stateSpace;
	}

	@Override
	public LocalDateTime getReferenceDate() {
		return forwardRateCurve.getReferenceDate() != null ? forwardRateCurve.getReferenceDate().atStartOfDay() : null;
	}

	@Override
	public RandomVariable getLIBOR(final MonteCarloProcess process, final int timeIndex, final int liborIndex)
			throws CalculationException {
		return process.getProcessValue(timeIndex, liborIndex);
	}

	@Override
	public TimeDiscretization getLiborPeriodDiscretization() {
		return liborPeriodDiscretization;
	}

	@Override
	public int getNumberOfLibors() {
		return liborPeriodDiscretization.getNumberOfTimeSteps();
	}

	@Override
	public double getLiborPeriod(final int timeIndex) {
		if(timeIndex >= liborPeriodDiscretization.getNumberOfTimes() || timeIndex < 0) {
			throw new ArrayIndexOutOfBoundsException("Index for LIBOR period discretization out of bounds: " + timeIndex + ".");
		}
		return liborPeriodDiscretization.getTime(timeIndex);
	}

	@Override
	public int getLiborPeriodIndex(final double time) {
		return liborPeriodDiscretization.getTimeIndex(time);
	}

	@Override
	public AnalyticModel getAnalyticModel() {
		return curveModel;
	}

	@Override
	public DiscountCurve getDiscountCurve() {
		return discountCurve;
	}

	@Override
	public ForwardCurve getForwardRateCurve() {
		return forwardRateCurve;
	}

	public RandomVariableFactory getRandomVariableFactory() {
		return randomVariableFactory;
	}

	@Override
	public int getNumberOfComponents() {
		return getNumberOfLibors();
	}

	@Override
	public RandomVariable[] getInitialState(final MonteCarloProcess process) {
		final RandomVariable[] liborInitialStates = new RandomVariable [liborPeriodDiscretization.getNumberOfTimeSteps()];
		for(int componentIndex = 0; componentIndex < liborPeriodDiscretization.getNumberOfTimeSteps(); componentIndex++) {
			final double rate = forwardRateCurve.getForward(
					curveModel, liborPeriodDiscretization.getTime(componentIndex), liborPeriodDiscretization.getTimeStep(componentIndex));
			liborInitialStates[componentIndex] = applyStateSpaceTransformInverse(process, 0, componentIndex, getRandomVariableForConstant(rate));
		}

		return liborInitialStates;
	}

	/**
	 * Returns the initial value of the model.
	 *
	 * @param process The discretization process generating this model. The process provides call backs for TimeDiscretization and allows calls to getProcessValue for timeIndices less or equal the given one.
	 * @return The initial value of the model.
	 */
	
	public RandomVariable[] getInitialValue(final MonteCarloProcess process) {
		final RandomVariable[] initialState = getInitialState(process);

		final RandomVariable[] value = new RandomVariable[initialState.length];
		for(int i= 0; i < value.length; i++) {
			value[i] = applyStateSpaceTransform(process, 0, i, initialState[i]);
		}

		return value;
	}

	@Override
	public RandomVariable getRandomVariableForConstant(final double value) {
		return randomVariableFactory.createRandomVariable(value);
	}

	public Measure getMeasure() {
		return measure;
	}

	public StateSpace getStateSpace() {
		return stateSpace;
	}
	
}
