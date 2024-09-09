package models;

import java.util.Map;
import java.util.stream.DoubleStream;

import measures.Measure;
import statespaces.DisplacedLognormalStateSpace;
import statespaces.LognormalStateSpace;
import statespaces.NormalStateSpace;
import statespaces.StateSpace;

import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.curves.DiscountCurve;
import net.finmath.marketdata.model.curves.DiscountCurveFromForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurveInterpolation;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.BrownianMotionFromMersenneRandomNumbers;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.interestrate.CalibrationProduct;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.LIBORMonteCarloSimulationFromLIBORModel;
import net.finmath.montecarlo.interestrate.TermStructureMonteCarloSimulationModel;
import net.finmath.montecarlo.interestrate.models.covariance.AbstractLIBORCovarianceModelParametric;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCorrelationModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCorrelationModelExponentialDecay;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModelFromVolatilityAndCorrelation;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModelTwoParameterExponentialForm;
import net.finmath.montecarlo.process.EulerSchemeFromProcessModel;
import net.finmath.montecarlo.process.MonteCarloProcess;
import net.finmath.stochastic.Scalar;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;

/**
 * This class creates the Euler-scheme simulation of the generalized LIBOR market model,
 * where we can specify the pricing measure(spot, terminal and or Tk-forward) and the state space dynamics under which the simulation is performed.
 *
 * @author Group 05
 */
public class GeneralLIBORMarketModelFactory {

	public static TermStructureMonteCarloSimulationModel createLIBORMarketModel(
			final RandomVariableFactory randomVariableFactory,
			final Measure measure,
			final StateSpace stateSpace,
			final double simulationTimeStep,
			final String simulationTimeInterpolationMethod,
			final double liborPeriodLength, // T_i-T_{i-1}, here assume to be fixed
			final double liborRateTimeHorizon, //T_n
			final double[] fixingForGivenForwards,
			final double[] givenForwards,
			final boolean ifUseDiscountCurve,
			final double[] covarianceParameters,	// the volatility parameter pair (\sigma, c), \sigma_{i}(t_{j}) = \sigma * exp(-c (T_{i}-t_{j})) and the decay of the correlation between LIBOR rates
			final int numberOfFactors,
			final int numberOfPaths,
			final int seed,
			final boolean ifUseLinearVolatilityModel
			) throws CalculationException {

		// Step 1: create the time discretization for the simulation of the processes
		final TimeDiscretization timeDiscretization = new TimeDiscretizationFromArray(
				0.0, (int) (liborRateTimeHorizon / simulationTimeStep), simulationTimeStep);

		// Step 2: create the time discretization for the tenure structure (i.e., the dates T_1,..,T_n)
		final int numberOfLIBOR = (int) (liborRateTimeHorizon / liborPeriodLength);
		final TimeDiscretization LIBORPeriodDiscretization = new
				TimeDiscretizationFromArray(0.0, numberOfLIBOR, liborPeriodLength);

		// Step 3: create the forward curve (initial values for the LIBOR market model)
		
		final ForwardCurve forwardCurve = ForwardCurveInterpolation.createForwardCurveFromForwards(
				"forwardCurve", // name of the curve
				fixingForGivenForwards, // fixings of the forward
				givenForwards, // the forwards we have
				liborPeriodLength	// tenor period length
				);

		// Step 4: create the discount curve from the corresponding forward curve; may be used for numeraire adjustment
		final DiscountCurve discountCurve = ifUseDiscountCurve? new DiscountCurveFromForwardCurve(forwardCurve) : null;

		// Step 5: create the two-parameter LIBOR volatility model, \sigma_{i}(t_{j}) = \sigma * exp(-c (T_{i}-t_{j}))
		final LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelTwoParameterExponentialForm(timeDiscretization, LIBORPeriodDiscretization,
				covarianceParameters[0], covarianceParameters[1]);

		// Step 6: create the LIBOR correlation model rho_{i,j} = exp(−alpha ∗ |T_i −T_j|)
		final LIBORCorrelationModel correlationModel = new LIBORCorrelationModelExponentialDecay(timeDiscretization, LIBORPeriodDiscretization,
				numberOfFactors, covarianceParameters[2]);

		// Step 7: create the LIBOR covariance model by combining volatility model and correlation model
		final LIBORCovarianceModel covarianceModel = new LIBORCovarianceModelFromVolatilityAndCorrelation(
				timeDiscretization, LIBORPeriodDiscretization, volatilityModel, correlationModel);

		//Step 8: create the model (i.e., the object of type LiborMarketModel)
		final LIBORMarketModel LIBORMarketModel;

		/*
		 * Set model properties
		 */
		final Map<String, String> properties = Map.of(
				"interpolationMethod", "linear",														// Interpolation of the tenor
				"simulationTimeInterpolationMethod", simulationTimeInterpolationMethod					// Interpolation of the simulation time
				);

		/*
		 *  Empty array of calibration items, to be given to the constructor of ExtendedLIBORMarketModelFromCovarianceModel
		 *  for model calibration.
		 */
		final CalibrationProduct[] calibrationItems = new CalibrationProduct[0];
		/*
		 * if the argument ifUSeLinearVolatilityModel is true, substitute the covariance model with a new one
		 * built on top of it depending on the state space specified. The new model is given by
		 * (L_i(t) + d_i) F_i where d_i is the displacement parameter for the i-th component of L,
		 * and F is the factor loading from the given covariance model.
		 *
		 * Otherwise, simulate the model as with suggested state space transformation without another
		 * covariance model on top of the underlying covariance model.
		 */
		if (!(stateSpace instanceof NormalStateSpace) && ifUseLinearVolatilityModel) {
			final LIBORCovarianceModel linearCovarianceModel;

			if (stateSpace instanceof DisplacedLognormalStateSpace) {
				linearCovarianceModel = new MultiCurveDisplacedLocalVolatilityModel(	// use our class implement displaced log-normal model allowing for different displacements
						(AbstractLIBORCovarianceModelParametric) covarianceModel,
						Scalar.arrayOf(((DisplacedLognormalStateSpace) stateSpace).getDisplacements()));
			} else {
				linearCovarianceModel = new MultiCurveDisplacedLocalVolatilityModel(
						(AbstractLIBORCovarianceModelParametric) covarianceModel, 0);	// zero displacement in case of log-normal state space
			}
			LIBORMarketModel = ExtendedLIBORMarketModelFromCovarianceModel.of(
					measure,
					new NormalStateSpace(),
					LIBORPeriodDiscretization,
					null /* analyticModel */, forwardCurve, discountCurve,
					randomVariableFactory, linearCovarianceModel, calibrationItems, properties);
		} else {
			LIBORMarketModel = ExtendedLIBORMarketModelFromCovarianceModel.of(
					measure,
					stateSpace,
					LIBORPeriodDiscretization,
					null /* analyticModel */, forwardCurve, discountCurve,
					randomVariableFactory, covarianceModel, calibrationItems, properties);
		}

		// Step 9: create an Euler scheme of the LIBOR model defined above, possibly with factor reduction
		final BrownianMotion brownianMotion = new BrownianMotionFromMersenneRandomNumbers(
				timeDiscretization,
				numberOfFactors, // may perform factor reduction here
				numberOfPaths,
				seed
				);

		final MonteCarloProcess process = new EulerSchemeFromProcessModel(
				LIBORMarketModel, brownianMotion, EulerSchemeFromProcessModel.Scheme.EULER);

		// Step 10: perform the simulation by giving the Euler scheme to the constructor of LIBORMonteCarloSimulationFromLIBORModel
		return new LIBORMonteCarloSimulationFromLIBORModel(process);
	}

	/*
	 * Some version with default arguments
	 *
	 * method to create the lognormal LIBOR market model
	 */
	public static TermStructureMonteCarloSimulationModel createLognormalLIBORMarketModel(
			final RandomVariableFactory randomVariableFactory,
			final Measure measure,
			final double simulationTimeStep,
			final double liborPeriodLength, //T_i-T_{i-1}, we suppose it to be fixed
			final double liborRateTimeHorizon, //T_n
			final double[] fixingForGivenForwards,
			final double[] givenForwards,
			final boolean ifUseDiscountCurve,
			final double[] covarianceParameters,	// the volatility parameter pair (\sigma, c), \sigma_{i}(t_{j}) = \sigma * exp(-c (T_{i}-t_{j})) and the decay of the correlation between LIBOR rates
			final int numberOfFactors,
			final int numberOfPaths,
			final int seed,
			final boolean ifUseLinearVolatilityModel
			) throws CalculationException {

		return createLIBORMarketModel(randomVariableFactory, measure, new LognormalStateSpace(),
				simulationTimeStep, "round_down" /* simulationTimeInterpolationMethod */,
				liborPeriodLength, liborRateTimeHorizon, fixingForGivenForwards, givenForwards, ifUseDiscountCurve,
				covarianceParameters,
				numberOfFactors, numberOfPaths, seed, ifUseLinearVolatilityModel);
	}

	/*
	 * method to create the displaced lognormal LIBOR market model
	 */
	public static TermStructureMonteCarloSimulationModel createDisplacedLognormalLIBORMarketModel(
			final RandomVariableFactory randomVariableFactory,
			final Measure measure,
			final double simulationTimeStep,
			final double liborPeriodLength, //T_i-T_{i-1}, we suppose it to be fixed
			final double liborRateTimeHorizon, //T_n
			final double[] fixingForGivenForwards,
			final double[] givenForwards,
			final boolean ifUseDiscountCurve,
			final double[] covarianceParameters,	// the volatility parameter pair (\sigma, c), \sigma_{i}(t_{j}) = \sigma * exp(-c (T_{i}-t_{j})) and the decay of the correlation between LIBOR rates
			final double[] displacements,
			final int numberOfFactors,
			final int numberOfPaths,
			final int seed,
			boolean ifUseLinearVolatilityModel
			) throws CalculationException {

		return createLIBORMarketModel(randomVariableFactory, measure, new DisplacedLognormalStateSpace(displacements),
				simulationTimeStep, "round_down" /* simulationTimeInterpolationMethod */,
				liborPeriodLength, liborRateTimeHorizon, fixingForGivenForwards, givenForwards, ifUseDiscountCurve,
				covarianceParameters,
				numberOfFactors, numberOfPaths, seed, ifUseLinearVolatilityModel);
	}

	/*
	 * method to create the displaced lognormal LIBOR market model with the same displacement for each forward rate
	 */
	public static TermStructureMonteCarloSimulationModel createBlendedLognormalLIBORMarketModel(
			RandomVariableFactory randomVariableFactory,
			Measure measure,
			double simulationTimeStep,
			double liborPeriodLength, //T_i-T_{i-1}, we suppose it to be fixed
			double liborRateTimeHorizon, //T_n
			double[] fixingForGivenForwards,
			double[] givenForwards,
			final boolean ifUseDiscountCurve,
			final double[] covarianceParameters,	// the volatility parameter pair (\sigma, c), \sigma_{i}(t_{j}) = \sigma * exp(-c (T_{i}-t_{j})) and the decay of the correlation between LIBOR rates
			double displacement,
			int numberOfFactors,
			int numberOfPaths,
			int seed,
			boolean ifUseLinearVolatilityModel
			) throws CalculationException {

		final int numberOfLIBOR = (int) (liborRateTimeHorizon / liborPeriodLength);
		final double[] displacements = DoubleStream.generate(()-> displacement).limit(numberOfLIBOR).toArray();

		return createLIBORMarketModel(randomVariableFactory, measure, new DisplacedLognormalStateSpace(displacements),
				simulationTimeStep, "round_down" /* simulationTimeInterpolationMethod */,
				liborPeriodLength, liborRateTimeHorizon, fixingForGivenForwards, givenForwards, ifUseDiscountCurve,
				covarianceParameters,
				numberOfFactors, numberOfPaths, seed, ifUseLinearVolatilityModel);
	}

	/*
	 * method to create the LIBOR market model with an initial forward rate curve flat
	 */
	public static TermStructureMonteCarloSimulationModel createLIBORMarketModelWithFlatInitialForward (
			RandomVariableFactory randomVariableFactory,
			Measure measure,
			StateSpace stateSpace,
			double simulationTimeStep,
			double liborPeriodLength, //T_i-T_{i-1}, we suppose it to be fixed
			double liborRateTimeHorizon, //T_n
			double forwardRate,
			final boolean ifUseDiscountCurve,
			final double[] covarianceParameters,	// the volatility parameter pair (\sigma, c), \sigma_{i}(t_{j}) = \sigma * exp(-c (T_{i}-t_{j})) and the decay of the correlation between LIBOR rates
			int numberOfFactors,
			int numberOfPaths,
			int seed,
			boolean ifUseLinearVolatilityModel
			) throws CalculationException {

		/*
		 * to be used to create the interest rate curve (here a flat (constant) value of forwardRate for the curve of forward rates
		 */
		double[] fixingForForwards = new double[] { 0.5, 1.0, 2.0, 5.0, 40.0};	/* fixings of the forward */
		double[] forwards = new double[] {forwardRate, forwardRate, forwardRate, forwardRate, forwardRate};	/* forwards */

		return createLIBORMarketModel(randomVariableFactory, measure, stateSpace, simulationTimeStep, "round_down" /* simulationTimeInterpolationMethod */,
				liborPeriodLength, liborRateTimeHorizon, fixingForForwards, forwards, ifUseDiscountCurve,
				covarianceParameters,
				numberOfFactors, numberOfPaths, seed, ifUseLinearVolatilityModel);
	}

}
