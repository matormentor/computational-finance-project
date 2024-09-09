package models;

import java.io.Serial;
import java.util.Map;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.models.covariance.AbstractLIBORCovarianceModelParametric;
import net.finmath.stochastic.RandomVariable;
import java.util.stream.DoubleStream;

import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.stochastic.Scalar;

/**
 * implements the multi-curve displaced local volatility model with general displacements for each forward rate component
 * <p>
 * The model constructed for the i-th factor loading is (L_i(t) + d_i) F_i(t)
 * where
 * d_i is the displacement for the i-th component of the stochastic process and
 * F_i(t) is the factor loading from the given covariance model.
 *
 * 
 * @author Group 05
 */
public class MultiCurveDisplacedLocalVolatilityModel extends AbstractLIBORCovarianceModelParametric {

	@Serial
	private static final long serialVersionUID = 1L;

	private final AbstractLIBORCovarianceModelParametric covarianceModel;
	private final RandomVariable[] displacements;

	/**
	 * Displaced local volatility model build on top of a standard covariance model.
	 * @param covarianceModel 	The given covariance model specifying the factor loadings.
	 * @param displacements 	The displacement d_0, ..., d_{n-1}.
	 */
	public MultiCurveDisplacedLocalVolatilityModel(final AbstractLIBORCovarianceModelParametric covarianceModel, final RandomVariable[] displacements) {
		super(covarianceModel.getTimeDiscretization(), covarianceModel.getLiborPeriodDiscretization(), covarianceModel.getNumberOfFactors());
		this.covarianceModel	= covarianceModel;
		this.displacements		= displacements;
	}

	
	public MultiCurveDisplacedLocalVolatilityModel(final AbstractLIBORCovarianceModelParametric covarianceModel, final double[] displacements) {
		super(covarianceModel.getTimeDiscretization(), covarianceModel.getLiborPeriodDiscretization(), covarianceModel.getNumberOfFactors());
		this.covarianceModel	= covarianceModel;
		this.displacements		= Scalar.arrayOf(displacements);
	}

	
	public MultiCurveDisplacedLocalVolatilityModel(final AbstractLIBORCovarianceModelParametric covarianceModel, final double displacement) {
		super(covarianceModel.getTimeDiscretization(), covarianceModel.getLiborPeriodDiscretization(), covarianceModel.getNumberOfFactors());
		this.covarianceModel	= covarianceModel;
		final int numberOfLibor = covarianceModel.getLiborPeriodDiscretization().getNumberOfTimeSteps();
		this.displacements		= Scalar.arrayOf(DoubleStream.generate(() -> displacement).limit(numberOfLibor).toArray());
	}

	/**
	 * Returns the base covariance model, i.e., the model providing the factor loading F
	 * such that this model's i-th factor loading is (L_i(t) + d_i) F_i(t)
	 * @return The covariance model.
	 */
	public AbstractLIBORCovarianceModelParametric getBaseCovarianceModel() {
		return covarianceModel;
	}

	@Override
	public RandomVariable[] getParameter() {
		return covarianceModel.getParameter();
	}

	@Override
	public double[] getParameterAsDouble() {
		final RandomVariable[] parameters = getParameter();
		final double[] parametersAsDouble = new double[parameters.length];
		for(int i = 0; i < parameters.length; i++) {
			parametersAsDouble[i] = parameters[i].doubleValue();
		}
		return parametersAsDouble;
	}

	@Override
	public Object clone() {
		return new MultiCurveDisplacedLocalVolatilityModel((AbstractLIBORCovarianceModelParametric) covarianceModel.clone(), displacements);
	}

	@Override
	public AbstractLIBORCovarianceModelParametric getCloneWithModifiedParameters(final RandomVariable[] parameters) {
		if(parameters == null || parameters.length == 0) {
			return this;
		}


		return new MultiCurveDisplacedLocalVolatilityModel(covarianceModel.getCloneWithModifiedParameters(parameters), displacements);
	}

	@Override
	public AbstractLIBORCovarianceModelParametric getCloneWithModifiedParameters(final double[] parameters) {
		return getCloneWithModifiedParameters(Scalar.arrayOf(parameters));
	}

	@Override
	public RandomVariable[] getFactorLoading(final int timeIndex, final int component, final RandomVariable[] realizationAtTimeIndex) {
		final RandomVariable[] factorLoading = covarianceModel.getFactorLoading(timeIndex, component, realizationAtTimeIndex);

		if(realizationAtTimeIndex != null && realizationAtTimeIndex[component] != null) {
			final RandomVariable localVolatilityFactor = realizationAtTimeIndex[component].add(displacements[component]);
			for (int factorIndex = 0; factorIndex < factorLoading.length; factorIndex++) {
				factorLoading[factorIndex] = factorLoading[factorIndex].mult(localVolatilityFactor);
			}
		}

		return factorLoading;
	}

	@Override
	public RandomVariable getFactorLoadingPseudoInverse(final int timeIndex, final int component, final int factor, final RandomVariable[] realizationAtTimeIndex) {
		throw new UnsupportedOperationException();
	}

	@Override
	public AbstractLIBORCovarianceModelParametric getCloneWithModifiedData(Map<String, Object> dataModified)
			throws CalculationException {
		RandomVariable[] displacements = this.displacements;
		AbstractLIBORCovarianceModelParametric covarianceModel = this.covarianceModel;
		RandomVariableFactory randomVariableFactory = null;

		if(dataModified != null) {
			if(dataModified.containsKey("randomVariableFactory")) {
				randomVariableFactory = (RandomVariableFactory)dataModified.get("randomVariableFactory");
				for (int i = 0; i < displacements.length; i++) {
					displacements[i] = randomVariableFactory.createRandomVariable((displacements[i]).doubleValue());
				}
			}
			if (!dataModified.containsKey("covarianceModel")) {
				covarianceModel = covarianceModel.getCloneWithModifiedData(dataModified);
			}
			covarianceModel = (AbstractLIBORCovarianceModelParametric)dataModified.getOrDefault("covarianceModel", covarianceModel);
			
			if (dataModified.getOrDefault("displacements", displacements) instanceof RandomVariable[]) {
				displacements = (RandomVariable[]) dataModified.getOrDefault("displacement", displacements);
			} else if (randomVariableFactory == null) {
				displacements = Scalar.arrayOf((double[]) dataModified.get("displacements"));
			} else {
				displacements = randomVariableFactory.createRandomVariableArray((double[]) dataModified.get("displacements"));
			}
		}

        return new MultiCurveDisplacedLocalVolatilityModel(covarianceModel, displacements);
	}

	public RandomVariable[] getDisplacements() {
		return displacements;
	}
}
