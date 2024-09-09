package info.quantlab.computationalfinance.assignments.montecarlo.interestrates.simplehybrid.product;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.stream.DoubleStream;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import info.quantlab.computationfinance.lecture.montecarlo.interestrates.project.Floater;
import measures.Measure;
import measures.SpotMeasure;
import measures.Tk_ForwardMeasure;
import measures.TerminalMeasure;
import models.GeneralLIBORMarketModelFactory;

import net.finmath.functions.AnalyticFormulas;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.RandomVariableFromArrayFactory;
import net.finmath.montecarlo.interestrate.TermStructureMonteCarloSimulationModel;
import net.finmath.montecarlo.interestrate.products.Bond;
import net.finmath.montecarlo.interestrate.products.Caplet;
import net.finmath.montecarlo.interestrate.products.TermStructureMonteCarloProduct;

/**
 * This class checks our implementation of extended LIBOR market model,
 * It tests the accuracy of the Forward and the Caplet under different measures using our extendedLIBOR market model.
 *
 * @author Group 05
 */
public class ValuationTest {

	private static final DecimalFormat formatterDeviation	= new DecimalFormat(" 0.000E00;-0.000E00", new DecimalFormatSymbols(Locale.ENGLISH));

	private static final RandomVariableFactory randomVariableFactory = new RandomVariableFromArrayFactory();

	//parameters for the two time discretizations
	private static final double simulationTimeStep = 0.5;
	private static final double periodLength = 0.5;				// semi-annual tenor discretization

	private static final double LIBORRateTimeHorizon = 20;
	private static final int numberOfLIBOR = (int) (LIBORRateTimeHorizon / periodLength);

	//fixing times for the forwards: the forwards corresponding to other fixing times will be interpolated
	private static final double forwardRate = 0.05;			// constant forward rate
	private static final double[] fixingForGivenForwards = { 0.5, 1.0, 2.0, 3.0};
	private static final double[] forwardsForCurve = { forwardRate, forwardRate, forwardRate, forwardRate};		// flat forward curve

	private static final double volatility = 0.30;					// constant volatility
	private static final double volatilityDecayParam = 0.0;			// no volatility exponential decay
	private static final double correlationDecayParam = 0.0;		// one factor, correlation of all drivers is 1
	private static final double[] covarianceParameters = new double[] {volatility, volatilityDecayParam, correlationDecayParam};
	private static final double displacement = -0.02;				// displaced log-normal model

	private static final int	numberOfPaths	= 50000;
	private static final int	numberOfFactors	= 1;
	private static final int	seed			= 3141;

	
	@Test
	public void testBondUnderMeasureAndStateSpace() throws Exception {
		final boolean useDiscountCurve = false;

		List<Measure> measures = new ArrayList<>();

		measures.add(new SpotMeasure());	// The spot measure is the same as the T_0 terminal measure
		
		for (int k = 0; k <= numberOfLIBOR; k += 5) {
			measures.add(new Tk_ForwardMeasure(k));
			
		}
		measures.add(new TerminalMeasure());	// The terminal measure is the same as the T_k forward measure with k = numberOfLIBOR

		for (final boolean ifUseLocalVolModel : new boolean[] { false, true } ) {
			for(final Measure measure : measures) {
				double errorSum = 0.0;

				final TermStructureMonteCarloSimulationModel simulationModel = GeneralLIBORMarketModelFactory.createDisplacedLognormalLIBORMarketModel(
						randomVariableFactory,
						measure,
						simulationTimeStep,
						periodLength,
						LIBORRateTimeHorizon,
						fixingForGivenForwards, forwardsForCurve, useDiscountCurve,
						covarianceParameters,
						DoubleStream.generate(()-> displacement).limit(numberOfLIBOR).toArray(),
						numberOfFactors,
						numberOfPaths,
						seed,
						ifUseLocalVolModel);

				for(double maturity = periodLength; maturity < LIBORRateTimeHorizon; maturity += periodLength) {
					final TermStructureMonteCarloProduct product = new Bond(maturity);
					final double value = product.getValue(simulationModel);
					final double yieldMonteCarlo = - Math.log(value) / maturity;

					final double valueAnalytic = 1.0 / Math.pow((1 + forwardRate * periodLength), maturity / periodLength);
					final double yieldAnalytic = - Math.log(valueAnalytic) / maturity;
					final double error		= yieldMonteCarlo - yieldAnalytic;
					errorSum += Math.abs(error);		// sum of error
				}
				final double averageError = errorSum / numberOfLIBOR;

				System.out.println("Zero bond average error using " + measure +
						(ifUseLocalVolModel? " under normal simulation with a linear displaced local volatility      ":
								" under displaced log-normal simulation without a linear local volatility")
						+ ": " + formatterDeviation.format(averageError));

				Assertions.assertEquals(0, averageError, 3E-4, "Testing implementation failed.");
			}
			System.out.println("____________________________________________________________________________________________________________________________________________________________________________________");
		}
	}

	@Test
	public void testForwardUnderMeasureAndStateSpace() throws Exception {
		final boolean useDiscountCurve = false;

		List<Measure> measures = new ArrayList<>();

		measures.add(new SpotMeasure());	
		for (int k = 0; k <= numberOfLIBOR; k += 5) {
			measures.add(new Tk_ForwardMeasure(k));
		}
		measures.add(new TerminalMeasure());	
		for (final boolean ifUseLocalVolModel : new boolean[] { false, true } ) {
			for(final Measure measure : measures) {

				double errorSum = 0.0;

				final TermStructureMonteCarloSimulationModel simulationModel = GeneralLIBORMarketModelFactory.createDisplacedLognormalLIBORMarketModel(
						randomVariableFactory,
						measure,
						simulationTimeStep,
						periodLength,
						LIBORRateTimeHorizon,
						fixingForGivenForwards, forwardsForCurve, useDiscountCurve,
						covarianceParameters,
						DoubleStream.generate(()-> displacement).limit(numberOfLIBOR).toArray(),
						numberOfFactors,
						numberOfPaths,
						seed,
						ifUseLocalVolModel);

				for(double fixing = periodLength; fixing < LIBORRateTimeHorizon; fixing += periodLength) {
					final TermStructureMonteCarloProduct floater = new Floater(fixing, fixing, fixing + periodLength, fixing + periodLength, periodLength);

					final double valueBondAnalytic = 1.0 / Math.pow((1 + forwardRate * periodLength), (fixing + periodLength) / periodLength);

					/*
					 * By discounting the floater with analytic bond value, it will reduce the numerical error for the spot measure
					 * where the error is accumulating with fixing but increase the error for the terminal measure where the error
					 * is kind of resulting from random noise.
					 */
					
					final double value = floater.getValue(simulationModel) / valueBondAnalytic;

					final double valueAnalytic = forwardRate * periodLength;
					final double error = (value - valueAnalytic);

					errorSum += Math.abs(error);			
				}

				final double averageError = errorSum / numberOfLIBOR;

				System.out.println("Forward rate average errors using " + measure +
						(ifUseLocalVolModel? " under normal simulation with a linear displaced local volatility      ":
								" under displaced log-normal simulation without a linear local volatility")
						+ ": " + formatterDeviation.format(averageError));

				/*
				 * The closer to the terminal error, the larger the forward rate error. Also, the use of a linear
				 * local volatility model without state space transformation increases the error.
				 */
				if (!((measure instanceof TerminalMeasure) || (measure instanceof Tk_ForwardMeasure &&
						(((Tk_ForwardMeasure) measure).getMaturityIndexTk() == numberOfLIBOR) ))){
					if (!(measure instanceof Tk_ForwardMeasure && (((Tk_ForwardMeasure) measure).getMaturityIndexTk() == numberOfLIBOR - 5))) {
						Assertions.assertEquals(0, averageError, 1E-4, "Testing implementation failed.");
					}
				} else {
					Assertions.assertEquals(0, averageError, 3E-4, "Testing implementation failed.");
				}
			}
			System.out.println("____________________________________________________________________________________________________________________________________________________________________________________");
		}
	}

	@Test
	public void testCapletSmileUnderMeasure() throws Exception {
		final boolean useDiscountCurve = false;

		List<Measure> measures = new ArrayList<>();

		measures.add(new SpotMeasure());	

		for (int k = 0; k <= numberOfLIBOR; k += 5) {
			measures.add(new Tk_ForwardMeasure(k));
		}
		measures.add(new TerminalMeasure());	

		for (final boolean ifUseLocalVolModel : new boolean[] { false, true } ) {
			for(final Measure measure : measures) {

				final TermStructureMonteCarloSimulationModel simulationModel = GeneralLIBORMarketModelFactory.createDisplacedLognormalLIBORMarketModel(
						randomVariableFactory,
						measure,
						simulationTimeStep,
						periodLength,
						LIBORRateTimeHorizon,
						fixingForGivenForwards, forwardsForCurve, useDiscountCurve,
						covarianceParameters,
						DoubleStream.generate(()-> displacement).limit(numberOfLIBOR).toArray(),
						numberOfFactors,
						numberOfPaths,
						seed,
						ifUseLocalVolModel);

				/*
				 * Value different products
				 */
				final double maturity = 5.0;

				double errorSum = 0;
				for(double strike = 0.025; strike <= 0.10; strike += 0.0025) {

					final TermStructureMonteCarloProduct product = new Caplet(maturity, periodLength, strike);
					final double value = product.getValue(simulationModel);

					/*
					 * Conversion to implied volatility
					 */
					
					// Determine the zero bond at payment numerically
					final TermStructureMonteCarloProduct bondAtPayment = new Bond(maturity + periodLength);
					final double discountFactor = bondAtPayment.getValue(simulationModel);

					// Determine the forward rate at fixing numerically
					final TermStructureMonteCarloProduct floaterAtMaturity = new Floater(maturity, periodLength, periodLength);
					final double forward = floaterAtMaturity.getValue(simulationModel) / discountFactor / periodLength;

                    final double payoffUnit = discountFactor * periodLength;

                    final double impliedVolatility = AnalyticFormulas.blackScholesOptionImpliedVolatility(forward + displacement, maturity, strike + displacement, payoffUnit, value);

					final double error = impliedVolatility - volatility;
					errorSum += Math.abs(error);			
				}

				final double averageError = errorSum / numberOfLIBOR;

				System.out.println("Caplet average implied volatility error using " + measure + "(ATM strike = " + forwardRate + ")" +
						(ifUseLocalVolModel? " under normal simulation with a linear displaced local volatility      ":
								" under displaced log-normal simulation without a linear local volatility")
						+ ": " + formatterDeviation.format(averageError));

				/*
				 * The closer to the terminal error, the larger the error. Also, the use of a linear
				 * local volatility model without state space transformation increases the error.
				 */
				
				if (!((measure instanceof TerminalMeasure) || (measure instanceof Tk_ForwardMeasure &&
						(((Tk_ForwardMeasure) measure).getMaturityIndexTk() == numberOfLIBOR) ))
						&& !ifUseLocalVolModel) {
					if (!(measure instanceof Tk_ForwardMeasure && (((Tk_ForwardMeasure) measure).getMaturityIndexTk() == numberOfLIBOR - 5))) {
						Assertions.assertEquals(0, averageError, 1E-3, "Testing implementation failed.");
					}
				} else {
					Assertions.assertEquals(0, averageError, 1E-2, "Testing implementation failed.");
				}
			}
			System.out.println("____________________________________________________________________________________________________________________________________________________________________________________");
		}
	}

	@Test
	public void testCapletATMImpliedVolUnderMeasure() throws Exception {

		final boolean useDiscountCurve = false;

		List<Measure> measures = new ArrayList<>();

		measures.add(new SpotMeasure());

		for (int k = 0; k <= numberOfLIBOR; k += 5) {
			measures.add(new Tk_ForwardMeasure(k));
		}
		measures.add(new TerminalMeasure());	
		
		for (final boolean ifUseLocalVolModel : new boolean[] { false, true } ) {
			for(final Measure measure : measures) {

				final TermStructureMonteCarloSimulationModel simulationModel = GeneralLIBORMarketModelFactory.createBlendedLognormalLIBORMarketModel(
						randomVariableFactory,
						measure,
						simulationTimeStep,
						periodLength,
						LIBORRateTimeHorizon,
						fixingForGivenForwards, forwardsForCurve, useDiscountCurve,
						covarianceParameters,
						displacement,
						numberOfFactors,
						numberOfPaths,
						seed,
						ifUseLocalVolModel);

				/*
				 * Value different products
				 */
				
				final double strike = forwardRate;
				double errorSum = 0.0;

				for(double maturity = periodLength; maturity < LIBORRateTimeHorizon; maturity += periodLength) {

					final TermStructureMonteCarloProduct product = new Caplet(maturity, periodLength, strike);
					final double value = product.getValue(simulationModel);

					/*
					 * Conversion to implied volatility
					 */
					
					// Determine the zero bond at payment (numerically)
					final TermStructureMonteCarloProduct bondAtPayment = new Bond(maturity + periodLength);
					final double discountFactor = bondAtPayment.getValue(simulationModel);

					// Determine the forward rate at fixing (numerically)
					final TermStructureMonteCarloProduct forwardRateProduct = new Floater(maturity, maturity, maturity + periodLength, maturity + periodLength, periodLength);
					final double forward = forwardRateProduct.getValue(simulationModel) / discountFactor / periodLength;

                    final double payoffUnit = discountFactor * periodLength;

                    final double impliedVolatility = AnalyticFormulas.blackScholesOptionImpliedVolatility(forward + displacement, maturity, strike + displacement, payoffUnit, value);

					final double error = (impliedVolatility - volatility);
					errorSum += Math.abs(error);			
				}

				final double averageError = errorSum / numberOfLIBOR;

				System.out.println("ATM Caplet average implied volatility error using " + measure +
						(ifUseLocalVolModel? " under normal simulation with a linear displaced local volatility      ":
								" under displaced log-normal simulation without a linear local volatility")
						+ ": " + formatterDeviation.format(averageError));

				Assertions.assertEquals(0, averageError, 8E-3, "Testing implementation failed.");
			}
			System.out.println("____________________________________________________________________________________________________________________________________________________________________________________");
		}
	}
}
