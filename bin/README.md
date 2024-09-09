# computational-finance-2022-2-project

**January 2023**

The project consists of the improvements of a Monte-Carlo simulation of an Euler-scheme discretization of a
Discrete Term Structure Model (LIBOR Market Model) and the investigation of properties of interest rate products
using this model.


## Project description

The current version of the project description can be found in the lecture's Moodle's announcements section.


## Importing as Maven project

The project comes with a pre-configured Maven pom.xml file referencing finmath-lib.

- Clone the project using git (`git clone`), then import the project in your favourite IDE as Maven project.


## Notes

### Data

You find the data for the calibration exercise in the class `LectureProjectData`. It is part of this project.

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


## Working in Eclipse

Import this git repository into Eclipse and start working.

- Click on the link to your repository (the link starts with qntlb/computational-finance… )
- Click on “Clone or download” and copy the URL to your clipboard.
- Go to Eclipse and select File -> Import -> Git -> Projects from Git.
- Select “Clone URI” and paste the GitHub URL from step 2.
- Select “master”, then Next -> Next
- In the Wizard for Project Import select “Import existing Eclipse projects”, then Next -> Finish
