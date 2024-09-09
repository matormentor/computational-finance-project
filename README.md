# computational-finance-2022-2-project

**2023**

The project consists of the improvements of a Monte-Carlo simulation of an Euler-scheme discretization of a
Discrete Term Structure Model (LIBOR Market Model) and the investigation of properties of interest rate products
using this model.


## Project description
Found in the project description [pdf](Lecture%202022.2%20Project%20-%20Discrete%20Term%20Structure%20Simulation%20and%20IR%20Product%20Valuation%20V0_7.pdf)

Link to the Powerpoint [Presentation](ProjectPresentation.pptx)


## Importing as Maven project

The project comes with a pre-configured Maven pom.xml file referencing finmath-lib.

- Clone the project using git (`git clone`), then import the project in your favourite IDE as Maven project.


## Notes

### Unit Tests

You may consider adding unit tests to your solution. You find a small stub / sample test in `src/main/test`. You can run unit test from Eclipse or via maven on the command line

```
mvn clean test
```

(run from the project directory).


### Code Style

If you like to improve your code-style you may run ``checkstyle`` via the maven command

```
mvn checkstyle:check 
```

(run from the project directory).

Checkstyle will report style issues of you code. If you like to clean up the formatting, you may use *Source -> Clean up...* in Eclipse.
