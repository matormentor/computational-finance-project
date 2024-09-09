package statespaces;

import net.finmath.stochastic.RandomVariable;
import net.finmath.stochastic.Scalar;

/**
 * Apply an abstract class to implement a state space with parameters
 *
 * @author Group 05
 */
public abstract class ParametricStateSpace implements StateSpace {

	/**
	 * Get the parameters of determining the parametric LIBOR market model.
	 *
	 * @return Parameter vector.
	 */
	public RandomVariable[]	getParameter() {
		final double[] parameterAsDouble = this.getParameterAsDouble();

		if (parameterAsDouble == null) {
			return null;
		}
		final RandomVariable[] parameter = new RandomVariable[parameterAsDouble.length];
		for(int i = 0; i < parameter.length; i++) {
			parameter[i] = new Scalar(parameterAsDouble[i]);
		}
		return parameter;
	}

	/**
	 * Get the parameters of determining this parametric LIBOR market model.
	 *
	 * @return Parameter vector.
	 */
	public abstract double[]	getParameterAsDouble();
	
	/**
	 * Return an instance of this state space using a new set of parameters.
	 *
	 * @param parameters The new set of parameters.
	 * @return An instance of ParametricStateSpace with modified parameters.
	 */
	
	public ParametricStateSpace getCloneWithModifiedParameters(final RandomVariable[] parameters) {
		final double[] parameterAsDouble = new double[parameters.length];

		for(int i = 0; i < parameterAsDouble.length; i++) {
			parameterAsDouble[i] = parameters[i].doubleValue();
		}
		return getCloneWithModifiedParameters(parameterAsDouble);
	}

	/**
	 * Return an instance of this model using a new set of parameters.
	 *
	 * @param parameters The new set of parameters.
	 * @return An instance of ParametricStateSpace with modified parameters.
	 */
	public abstract ParametricStateSpace getCloneWithModifiedParameters(double[] parameters);
}
