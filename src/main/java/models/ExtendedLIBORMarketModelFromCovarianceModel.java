package models;

import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

import measures.Measure;
import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.AnalyticModel;
import net.finmath.marketdata.model.curves.DiscountCurve;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.interestrate.CalibrationProduct;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.models.covariance.AbstractLIBORCovarianceModelParametric;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.TimeDiscretization;
import statespaces.ParametricStateSpace;
import statespaces.StateSpace;

/**
 * This class duplicates the class LIBORMarketModelFromCovarianceModel from the finmath lib to implement a
 * generalized LIBOR market model using different state spaces and the universal pricing measure
 * under which it is to be simulated.
 *
 * @author Group 05
 */
public class ExtendedLIBORMarketModelFromCovarianceModel extends AbstractExtendedLIBORMarketModel {


	public enum Driftapproximation	{ EULER, LINE_INTEGRAL }
	public enum InterpolationMethod	{ LINEAR, LOG_LINEAR_UNCORRECTED, LOG_LINEAR_CORRECTED }
	public enum SimulationTimeInterpolationMethod { ROUND_DOWN, ROUND_NEAREST }

	private final Measure 	measure;
	private final StateSpace 	stateSpace;

	private LIBORCovarianceModel 			covarianceModel;

	private SimulationTimeInterpolationMethod	simulationTimeInterpolationMethod		= SimulationTimeInterpolationMethod.ROUND_NEAREST;
	private InterpolationMethod					interpolationMethod					= InterpolationMethod.LOG_LINEAR_UNCORRECTED;

	private double				liborCap					= 1E5;

	// This is a cache of the integrated covariance.
	private double[][][]		integratedLIBORCovariance;

	// Cache for the numéraires, needs to be invalidated if process changes - move out of the object
	private final transient ConcurrentHashMap<Integer, RandomVariable>	numeraires = new ConcurrentHashMap<>();

	/**
	 * Creates a  LIBOR Market Model for given covariance model.
	 * 
	 * @param measure 		The equivalent martingale measure under which the simulation is performed
	 * @param stateSpace 	The state space transform for the discretization scheme
	 * @param liborPeriodDiscretization The discretization of the interest rate curve into forward rates (tenor structure).
	 * @param curveModel The associated analytic model of this model (containing the associated market data objects like curve).
	 * @param forwardRateCurve 	The initial values for the forward rates.
	 * @param discountCurve 	The discount curve to use. This will create an LMM model with a deterministic zero-spread discounting adjustment.
	 * @param randomVariableFactory The random variable factory used to create the initial values of the model.
	 * @param covarianceModel 		The covariance model to use.
	 * @param calibrationProducts 	The vector of calibration items (a union of a product, target value and weight) for the objective function sum weight(i) * (modelValue(i)-targetValue(i)).
	 * @param properties 	Key value map specifying other model properties and calibration (hyper)parameters.
	 * @throws net.finmath.exception.CalculationException Thrown if the valuation fails, specific cause may be available via the <code>cause()</code> method.
	 */
	public ExtendedLIBORMarketModelFromCovarianceModel(
			final Measure measure,
			final StateSpace stateSpace,
			final TimeDiscretization	liborPeriodDiscretization,
			final AnalyticModel			curveModel,
			final ForwardCurve			forwardRateCurve,
			final DiscountCurve			discountCurve,
			final RandomVariableFactory	randomVariableFactory,
			final LIBORCovarianceModel	covarianceModel,
			final CalibrationProduct[]	calibrationProducts,
			final Map<String, ?>		properties
			) throws CalculationException {
		super(liborPeriodDiscretization, curveModel, forwardRateCurve, discountCurve, randomVariableFactory,measure,stateSpace);

		this.measure = measure;
		this.stateSpace = stateSpace;

		// Set some properties
		if(properties != null && properties.containsKey("interpolationMethod")) {
			interpolationMethod	= InterpolationMethod.valueOf(((String)properties.get("interpolationMethod")).toUpperCase());
		}
		if(properties != null && properties.containsKey("simulationTimeInterpolationMethod")) {
			simulationTimeInterpolationMethod	= SimulationTimeInterpolationMethod.valueOf(((String)properties.get("simulationTimeInterpolationMethod")).toUpperCase());
		}
		if(properties != null && properties.containsKey("liborCap")) {
			liborCap	= (Double)properties.get("liborCap");
		}
		this.covarianceModel	= covarianceModel;
		
	}

	/**
	 * Creates a LIBOR Market Model for given covariance.
	 *
	 * @param measure 		The equivalent martingale measure under which the simulation is performed
	 * @param stateSpace 	The state space transform for the discretization scheme
	 * @param liborPeriodDiscretization The discretization of the interest rate curve into forward rates (tenor structure).
	 * @param analyticModel 	The associated analytic model of this model (containing the associated market data objects like curve).
	 * @param forwardRateCurve 	The initial values for the forward rates.
	 * @param discountCurve 	The discount curve to use. This will create an LMM model with a deterministic zero-spread discounting adjustment.
	 * @param randomVariableFactory The random variable factory used to create the initial values of the model.
	 * @param covarianceModel 	The covariance model to use.
	 * @param properties Key value map specifying other model properties and calibration (hyper)parameters.
	 * @throws net.finmath.exception.CalculationException Thrown if the valuation fails, specific cause may be available via the <code>cause()</code> method.
	 */
	public ExtendedLIBORMarketModelFromCovarianceModel(
			final Measure measure,
			final StateSpace stateSpace,
			final TimeDiscretization	liborPeriodDiscretization,
			final AnalyticModel			analyticModel,
			final ForwardCurve			forwardRateCurve,
			final DiscountCurve			discountCurve,
			final RandomVariableFactory	randomVariableFactory,
			final LIBORCovarianceModel	covarianceModel,
			final Map<String, ?>		properties
			) throws CalculationException {
		this(measure, stateSpace, liborPeriodDiscretization, analyticModel, forwardRateCurve, discountCurve, randomVariableFactory, covarianceModel, null, properties);
	}

	/**
	 * Creates a LIBOR Market Model for given covariance with a calibration (if calibration items are given).
	 * 
	 * @param measure 		The equivalent martingale measure under which the simulation is performed
	 * @param stateSpace 	The state space transform for the discretization scheme
	 * @param liborPeriodDiscretization The discretization of the interest rate curve into forward rates (tenor structure).
	 * @param analyticModel 	The associated analytic model of this model (containing the associated market data objects like curve).
	 * @param forwardRateCurve 	The initial values for the forward rates.
	 * @param discountCurve 	The discount curve to use. This will create an LMM model with a deterministic zero-spread discounting adjustment.
	 * @param randomVariableFactory The random variable factory used to create the initial values of the model.
	 * @param covarianceModel 		The covariance model to use.
	 * @param calibrationProducts 	The vector of calibration items (a union of a product, target value and weight) for the objective function sum weight(i) * (modelValue(i)-targetValue(i).
	 * @param properties Key value map specifying other model properties and calibration (hyper)parameters.
	 * @return A new instance of LIBORMarketModelFromCovarianceModel, possibly calibrated.
	 * @throws net.finmath.exception.CalculationException Thrown if the valuation fails, specific cause may be available via the <code>cause()</code> method.
	 */
	public static ExtendedLIBORMarketModelFromCovarianceModel of(
			final Measure measure,
			final StateSpace stateSpace,
			final TimeDiscretization	liborPeriodDiscretization,
			final AnalyticModel			analyticModel,
			final ForwardCurve			forwardRateCurve,
			final DiscountCurve			discountCurve,
			final RandomVariableFactory	randomVariableFactory,
			final LIBORCovarianceModel	covarianceModel,
			final CalibrationProduct[]	calibrationProducts,
			final Map<String, ?>		properties
			) throws CalculationException {

        return new ExtendedLIBORMarketModelFromCovarianceModel(
                measure, stateSpace, liborPeriodDiscretization, analyticModel, forwardRateCurve, discountCurve, randomVariableFactory, covarianceModel, properties);
		
	}

	public double[] getParameterAsDouble() {
		final double[] covarianceParameters = ((AbstractLIBORCovarianceModelParametric) covarianceModel).getParameterAsDouble();
		
		/*
		 * Check if it is displaced lognomal
		 */
		if(!(stateSpace instanceof ParametricStateSpace) || ((ParametricStateSpace) stateSpace).getParameterAsDouble() == null) {
			return covarianceParameters;
		}
		final double[] displacementParameters = ((ParametricStateSpace) stateSpace).getParameterAsDouble();
		if (covarianceParameters == null) {
			return displacementParameters;
		}

		// Append displacement to the end of covarianceParameters
		final double[] jointParameters = new double[covarianceParameters.length + displacementParameters.length];
		System.arraycopy(covarianceParameters, 0, jointParameters, 0, covarianceParameters.length);
		System.arraycopy(displacementParameters, 0, jointParameters, covarianceParameters.length, displacementParameters.length);
		
		return jointParameters;
	}

	@Override
	public LIBORMarketModel getCloneWithModifiedCovarianceModel(final LIBORCovarianceModel calibrationCovarianceModel) {
		final ExtendedLIBORMarketModelFromCovarianceModel model = (ExtendedLIBORMarketModelFromCovarianceModel) this.clone();
		model.covarianceModel = calibrationCovarianceModel;
		return model;
	}

	@Override
	public RandomVariable getForwardRate(final MonteCarloProcess process, double time, final double periodStart, final double periodEnd)
			throws CalculationException {
		final int periodStartIndex    = getLiborPeriodIndex(periodStart);
		final int periodEndIndex      = getLiborPeriodIndex(periodEnd);

		// If time is beyond fixing, use the fixing time.
		time = Math.min(time, periodStart);
		int timeIndex           = process.getTimeIndex(time);

		// Simulation time interpolation: If time is not part of the discretization, use the nearest available point.
		if(timeIndex < 0) {
			timeIndex = -timeIndex - 2;
			if(simulationTimeInterpolationMethod == SimulationTimeInterpolationMethod.ROUND_NEAREST && time-process.getTime(timeIndex) > process.getTime(timeIndex+1)-time) {
				timeIndex++;
			}
		}

		// The forward rates are provided on fractional tenor discretization points using linear interpolation

		// Interpolation on tenor using interpolationMethod
		// Interpolation on the period end of the forward rate
		
		if(periodEndIndex < 0) {
			final int		previousEndIndex	= (-periodEndIndex - 1) - 1;
			final double	nextEndTime			= getLiborPeriod(previousEndIndex + 1);

			// Interpolate libor from periodStart to periodEnd on periodEnd
			final RandomVariable onePlusLongLIBORdt         = getForwardRate(process, time, periodStart, nextEndTime).mult(nextEndTime - periodStart).add(1.0);
			final RandomVariable onePlusInterpolatedLIBORDt = getOnePlusInterpolatedLIBORDt(process, timeIndex, periodEnd, previousEndIndex);
			return onePlusLongLIBORdt.div(onePlusInterpolatedLIBORDt).sub(1.0).div(periodEnd - periodStart);
		}

		// Interpolation on the period start of the forward rate
		if(periodStartIndex < 0) {
			final int	previousStartIndex   = (-periodStartIndex - 1) - 1;
			final double nextStartTime	 = getLiborPeriod(previousStartIndex + 1);
			if(nextStartTime > periodEnd) {
				throw new AssertionError("Interpolation is not possible.");
			}
			if(nextStartTime == periodEnd) {
				return getOnePlusInterpolatedLIBORDt(process, timeIndex, periodStart, previousStartIndex).sub(1.0).div(periodEnd - periodStart);
			}
			final RandomVariable onePlusLongLIBORdt         = getForwardRate(process, time, nextStartTime, periodEnd).mult(periodEnd - nextStartTime).add(1.0);
			final RandomVariable onePlusInterpolatedLIBORDt = getOnePlusInterpolatedLIBORDt(process, timeIndex, periodStart, previousStartIndex);
			return onePlusLongLIBORdt.mult(onePlusInterpolatedLIBORDt).sub(1.0).div(periodEnd - periodStart);
		}

		if(periodStartIndex < 0 || periodEndIndex < 0) {
			throw new AssertionError("LIBOR requested outside libor discretization points and interpolation was not performed.");
		}

		// If this is a model primitive then return it
		if(periodStartIndex + 1 == periodEndIndex) {
			return getLIBOR(process, timeIndex, periodStartIndex);
		}

		RandomVariable accrualAccount = null;

		// Calculate the value of the forward bond
		for(int periodIndex = periodStartIndex; periodIndex < periodEndIndex; periodIndex++)
		{
			final double subPeriodLength = getLiborPeriod(periodIndex + 1) - getLiborPeriod(periodIndex);
			final RandomVariable liborOverSubPeriod = getLIBOR(process, timeIndex, periodIndex);

			accrualAccount = accrualAccount == null ? liborOverSubPeriod.mult(subPeriodLength).add(1.0) : accrualAccount.accrue(liborOverSubPeriod, subPeriodLength);
		}

		final RandomVariable libor = accrualAccount.sub(1.0).div(periodEnd - periodStart);

		return libor;
	}

	/**
	 * Implement the interpolation of the forward rate in tenor time.
	 * The method provides the forward rate \( F(t_{i}, S, T_{j+1}) \) where \( S \in [T_{j}, T_{j+1}] \).
	 *
	 * @param timeIndex 		The time index associated with the simulation time. The index i in \( t_{i} \).
	 * @param periodStartTime 	The period start time S (on which we interpolate).
	 * @param liborPeriodIndex 	The period index j for which \( S \in [T_{j}, T_{j+1}] \) (to avoid another lookup).
	 * @return The interpolated forward rate.
	 * @throws CalculationException Thrown if valuation failed.
//	 */
	private RandomVariable getOnePlusInterpolatedLIBORDt(final MonteCarloProcess process, int timeIndex, final double periodStartTime, final int liborPeriodIndex) throws CalculationException
	{
		final double tenorPeriodStartTime       = getLiborPeriod(liborPeriodIndex);
		final double tenorPeriodEndTime         = getLiborPeriod(liborPeriodIndex + 1);
		final double tenorDt                    = tenorPeriodEndTime - tenorPeriodStartTime;
		if(tenorPeriodStartTime < process.getTime(timeIndex)) {
			// Fixed at Long LIBOR period Start.
			timeIndex  = Math.min(timeIndex, process.getTimeIndex(tenorPeriodStartTime));
			if(timeIndex < 0) {
				//				timeIndex = -timeIndex-2;			// mapping to last known fixing.
				throw new IllegalArgumentException("Tenor discretization not part of time discretization.");
			}
		}
		final RandomVariable onePlusLongLIBORDt = getLIBOR(process, timeIndex , liborPeriodIndex).mult(tenorDt).add(1.0);

		final double smallDt                    = tenorPeriodEndTime - periodStartTime;
		final double alpha                      = smallDt / tenorDt;

		RandomVariable onePlusInterpolatedLIBORDt;
		switch(interpolationMethod)
		{
		case LINEAR:
			onePlusInterpolatedLIBORDt = onePlusLongLIBORDt.mult(alpha).add(1 - alpha);
			break;
		case LOG_LINEAR_UNCORRECTED:
		case LOG_LINEAR_CORRECTED:
			onePlusInterpolatedLIBORDt = onePlusLongLIBORDt.log().mult(alpha).exp();
			break;
		default: throw new IllegalArgumentException("Method for enum " + interpolationMethod.name() + " not implemented!");
		}

	
		final double analyticOnePlusLongLIBORDt   = 1 + getForwardRateCurve().getForward(getAnalyticModel(), tenorPeriodStartTime, tenorDt) * tenorDt;
		final double analyticOnePlusShortLIBORDt	= 1 + getForwardRateCurve().getForward(getAnalyticModel(), periodStartTime, smallDt) * smallDt;

		double analyticOnePlusInterpolatedLIBORDt;
		switch(interpolationMethod)
		{
		case LINEAR:
			analyticOnePlusInterpolatedLIBORDt = analyticOnePlusLongLIBORDt * alpha + (1-alpha);
			break;
		case LOG_LINEAR_UNCORRECTED:
		case LOG_LINEAR_CORRECTED:
			analyticOnePlusInterpolatedLIBORDt = Math.exp(Math.log(analyticOnePlusLongLIBORDt) * alpha);
			break;
		default: throw new IllegalArgumentException("Method for enum " + interpolationMethod.name() + " not implemented!");
		}
		onePlusInterpolatedLIBORDt = onePlusInterpolatedLIBORDt.mult(analyticOnePlusShortLIBORDt / analyticOnePlusInterpolatedLIBORDt);

		return onePlusInterpolatedLIBORDt;
	}

	@Override
	public RandomVariable applyStateSpaceTransform(final MonteCarloProcess process, final int timeIndex, final int componentIndex,
			final RandomVariable randomVariable) {
		RandomVariable value = stateSpace.applyStateSpaceTransform(process, timeIndex, componentIndex, randomVariable);

		return value;
	}

	@Override
	public RandomVariable applyStateSpaceTransformInverse(final MonteCarloProcess process, final int timeIndex, final int componentIndex, final RandomVariable randomVariable) {
		RandomVariable value = stateSpace.applyStateSpaceTransformInverse(process, timeIndex, componentIndex, randomVariable);

		return value;
	}

	/**
	 * Return the numeraire at a given time.
	 * 
	 * @param time Time time <i>t</i> for which the numeraire should be returned <i>N(t)</i>.
	 * @return The numeraire at the specified time as <code>RandomVariable</code>
	 * @throws net.finmath.exception.CalculationException Thrown if the valuation fails, specific cause may be available via the <code>cause()</code> method.
	 */
	@Override
	public RandomVariable getNumeraire(final MonteCarloProcess process, final double time) throws CalculationException {
		if(time < 0) {
			return getRandomVariableFactory().createRandomVariable(getDiscountCurve().getDiscountFactor(getAnalyticModel(), time));
		}

		RandomVariable numeraire = getNumeraireUnAdjusted(process, time);	// N*(T)
		return numeraire;
	}

	private RandomVariable getNumeraireUnAdjusted(final MonteCarloProcess process, final double time) throws CalculationException {
		return measure.getNumeraireUnadjusted(process, time);
	}

	protected RandomVariable getNumerairetUnAdjustedAtLIBORIndex(final MonteCarloProcess process, final int liborTimeIndex) throws CalculationException {
		/*
		 * synchronize lazy init cache
		 */
		synchronized(numeraires) {
			
			RandomVariable numeraireUnadjusted = numeraires.get(liborTimeIndex);

			if (numeraireUnadjusted == null) {
				numeraireUnadjusted = measure.getNumeraireUnadjusted(process, getLiborPeriodDiscretization().getTime(liborTimeIndex));
				numeraires.put(liborTimeIndex, numeraireUnadjusted);
			}
			return numeraireUnadjusted;
		}
	}

	/**
	 * Return the complete vector of the drift for the time index timeIndex, given that current state is realizationAtTimeIndex.
	 * The drift is zero for rates being already fixed.
	 *
	 * The method currently provides the drift for any equivalent martingale measure - depending how the model object was constructed.
	 *
	 * @param process    The discretization process generating this model. The process provides call backs for TimeDiscretization and allows calls to getProcessValue for timeIndices less or equal the given one.
	 * @param timeIndex  Time index i for which the drift should be returned.
	 * @param realizationAtTimeIndex Time current forward rate vector at time index i which should be used in the calculation.
	 * @return The drift vector μ(t_{i}) as type RandomVariableFromDoubleArray
	 */
	@Override
	public RandomVariable[] getDrift(final MonteCarloProcess process, final int timeIndex, final RandomVariable[] realizationAtTimeIndex,
			final RandomVariable[] realizationPredictor) {
		return measure.getDrift(process, timeIndex, realizationAtTimeIndex, realizationPredictor, stateSpace);
	}

	@Override
	public RandomVariable[] getFactorLoading(final MonteCarloProcess process, final int timeIndex, final int componentIndex,
			RandomVariable[] realizationAtTimeIndex) {
		return covarianceModel.getFactorLoading(process.getTime(timeIndex), getLiborPeriod(componentIndex), realizationAtTimeIndex);
	}

	@Override
	public double[][][] getIntegratedLIBORCovariance(final TimeDiscretization simulationTimeDiscretization) {
			if(integratedLIBORCovariance == null) {
				final TimeDiscretization liborPeriodDiscretization = getLiborPeriodDiscretization();

				integratedLIBORCovariance = new double[simulationTimeDiscretization.getNumberOfTimeSteps()][liborPeriodDiscretization.getNumberOfTimeSteps()][liborPeriodDiscretization.getNumberOfTimeSteps()];
				for(int timeIndex = 0; timeIndex < simulationTimeDiscretization.getNumberOfTimeSteps(); timeIndex++) {
					final double dt = simulationTimeDiscretization.getTimeStep(timeIndex);
					final RandomVariable[][] factorLoadings = new RandomVariable[liborPeriodDiscretization.getNumberOfTimeSteps()][];
					
					for(int componentIndex = 0; componentIndex < liborPeriodDiscretization.getNumberOfTimeSteps(); componentIndex++) {
						factorLoadings[componentIndex] = covarianceModel.getFactorLoading(
								simulationTimeDiscretization.getTime(timeIndex), liborPeriodDiscretization.getTime(componentIndex), null);
					}

					for(int componentIndex1 = 0; componentIndex1 < liborPeriodDiscretization.getNumberOfTimeSteps(); componentIndex1++) {
						final RandomVariable[] factorLoadingOfComponent1 = factorLoadings[componentIndex1];
						// Sum the libor cross terms
						for(int componentIndex2 = componentIndex1; componentIndex2 < liborPeriodDiscretization.getNumberOfTimeSteps(); componentIndex2++) {
							double integratedLIBORCovarianceValue = 0.0;
							if(getLiborPeriod(componentIndex1) > simulationTimeDiscretization.getTime(timeIndex)) {
								final RandomVariable[] factorLoadingOfComponent2 = factorLoadings[componentIndex2];
								for(int factorIndex = 0; factorIndex < factorLoadingOfComponent2.length; factorIndex++) {
									integratedLIBORCovarianceValue += factorLoadingOfComponent1[factorIndex].get(0) * factorLoadingOfComponent2[factorIndex].get(0) * dt;
								}
							}
							integratedLIBORCovariance[timeIndex][componentIndex1][componentIndex2] = integratedLIBORCovarianceValue;
							integratedLIBORCovariance[timeIndex][componentIndex2][componentIndex1] = integratedLIBORCovarianceValue; // use symmetry
						}
					}
				}

				// Integrate over time (i.e. sum up).
				for(int timeIndex = 1; timeIndex < simulationTimeDiscretization.getNumberOfTimeSteps(); timeIndex++) {
					final double[][] prevIntegratedLIBORCovariance = integratedLIBORCovariance[timeIndex - 1];
					final double[][] thisIntegratedLIBORCovariance = integratedLIBORCovariance[timeIndex];
					for(int componentIndex1 = 0; componentIndex1 < liborPeriodDiscretization.getNumberOfTimeSteps(); componentIndex1++) {
						for(int componentIndex2 = componentIndex1; componentIndex2 < liborPeriodDiscretization.getNumberOfTimeSteps(); componentIndex2++) {
							thisIntegratedLIBORCovariance[componentIndex1][componentIndex2] = prevIntegratedLIBORCovariance[componentIndex1][componentIndex2] + thisIntegratedLIBORCovariance[componentIndex1][componentIndex2];
							thisIntegratedLIBORCovariance[componentIndex2][componentIndex1] = thisIntegratedLIBORCovariance[componentIndex1][componentIndex2];
						}
					}
				}
			}

		return integratedLIBORCovariance;
	}

	@Override
	public ExtendedLIBORMarketModelFromCovarianceModel getCloneWithModifiedData(final Map<String, Object> dataModified) throws CalculationException {
		Measure measure										= this.measure;
		StateSpace stateSpace 								= this.stateSpace;
		RandomVariableFactory 	randomVariableFactory		= getRandomVariableFactory();
		TimeDiscretization		liborPeriodDiscretization	= getLiborPeriodDiscretization();
		AnalyticModel			analyticModel				= getAnalyticModel();
		ForwardCurve			forwardRateCurve			= getForwardRateCurve();
		DiscountCurve			discountCurve				= getDiscountCurve();
		LIBORCovarianceModel	covarianceModel				= getCovarianceModel();

		final Map<String, Object> properties = new HashMap<>();

		properties.put("interpolationMethod", interpolationMethod.name());
		properties.put("liborCap", liborCap);

		if(dataModified != null) {
			
			measure = (Measure) dataModified.getOrDefault("measure", measure);
			stateSpace = (StateSpace) dataModified.getOrDefault("stateSpace", stateSpace);
			randomVariableFactory = (RandomVariableFactory)dataModified.getOrDefault("randomVariableFactory", randomVariableFactory);
			liborPeriodDiscretization = (TimeDiscretization)dataModified.getOrDefault("liborPeriodDiscretization", liborPeriodDiscretization);
			analyticModel = (AnalyticModel)dataModified.getOrDefault("analyticModel", analyticModel);
			forwardRateCurve = (ForwardCurve)dataModified.getOrDefault("forwardRateCurve", forwardRateCurve);
			discountCurve = (DiscountCurve)dataModified.getOrDefault("discountCurve", discountCurve);
			covarianceModel = (LIBORCovarianceModel)dataModified.getOrDefault("covarianceModel", covarianceModel);
		}

		final ExtendedLIBORMarketModelFromCovarianceModel newModel = ExtendedLIBORMarketModelFromCovarianceModel.of(measure, stateSpace, liborPeriodDiscretization, analyticModel, forwardRateCurve, discountCurve, randomVariableFactory, covarianceModel, null, properties);
		return newModel;
	}

	@Override
	public Object clone() {
		try {
			final Map<String, Object>				properties					= new HashMap<>();
			properties.put("interpolationMethod", interpolationMethod.name());
			properties.put("liborCap", liborCap);
			return ExtendedLIBORMarketModelFromCovarianceModel.of(
					getMeasure(), getStateSpace(), getLiborPeriodDiscretization(), getAnalyticModel(), getForwardRateCurve(), getDiscountCurve(),
					getRandomVariableFactory(), covarianceModel, null, properties);
		} catch (final CalculationException e) {
			return null;
		}
	}

	@Override
	public Measure getMeasure() {
		return measure;
	}

	@Override
	public StateSpace getStateSpace() {
		return stateSpace;
	}

	@Override
	public LIBORCovarianceModel getCovarianceModel() {
		return covarianceModel;
	}

	@Override
	public int getNumberOfFactors() {
		return covarianceModel.getNumberOfFactors();
	}
}

