To make sure we would not miss any pattern, we have tried all statistical models that are relevant to the continuous dependent variable and are available in R through the caret package. 

In total, we attempted to use 88 models, from which we could not perform 7 models for various technical reasons.

From this set, we fitted each model to each dependent variable in our dataset. For models that require tunning parameters we used adaptive tuning procedure using 10-fold cross-validation. This procedure was used for each model 10 times to get best error estimates of the models’ quality.
For the quality of fit, we used RMSE - Root Mean Standard Error. We always calculated this statistic on the testing set using cross-validation.

In total, we attempted to fit

* 26 Linear Regressions
* 7 Support Vector Machines
* 2 Relevance Vector Machines (counted as Support Vector Machines later on)
* 3 Gaussian Processes
* 7 Random Forests
* 8 Tree-Based Models (counted together with Random Forests later on)
* 4 Multivariate Adaptive Regression Splines
* 11 Neural Networks

After fitting all the models, we kept all the models which have the quality of fit - RMSE - within interquartile range of the best performing model.

Attached picture shows RMSE for the models.

As can be seen, there are no clear winners - there are many completely different models that fit the data just as good as the best one.


##  Variable importances

In lieu of coefficient significances of ordinary least squares regression, we computed permutation-based feature (variable) importances, as described in https://arxiv.org/abs/1801.01489 using the _DALEX_ library for R. In short, variable's importance is a difference between the objective function (here: RMSA) of the model on the whole datasaet and the objective function calculated on the modified dataset, where the values of the variable are randomly permuted. If the variable has no explanatory power, replacing it with its random permutation should not reduce the power of the model. Distribution of variable importances, if calculated on a sufficiently large number of bootstrapped samples can be used as a basis of non-parametric test for regression coeffient significance.

For each variable in the model we gathered 1000 samples of the permutation importances. 

Boxplots of importances with the LASSO regression coefficients superimposed can be seen on the charts <chart references>. Importances on the chart are scaled in such a way that 100% accounts for the whole model predictive power measured by the RMSA. They are illustrated as boxplots sorted by their mean, and annotated with the LASSO regression coeffients. Whiskers of the boxes span an interquartile range of the distribution of the importances for a given variable. 
​

## Which variables are relevant to analysis? Do the models treat use the variables in simlar way?

After gathering all the models, we calculated the variables’ importances using library _DALEX_, which is an implementation of the idea described in: https://arxiv.org/abs/1801.01489 . In short, for a given predictor, we see what is the expected loss of the prediction if we replace the predictor with its random permutation.

We get estimates of the variable relevance for each model and each predictor. To visualize if there are groups of models in respect to the variable importances, we performed a multidimensional scaling.

A’priori we expected that models from the same family would be grouped together. It turned out, that in general all the models are in agreement, which reinforced our belief, that in this problem choice of the model is not critical. Because there is a difficulty in visualizing 62 models in one chart aesthetically, we reduced the number of models to include only regressions. Still, the charts are barely readable.

Size of the circle (indicated as weight_cnt) measures amount of variables that are relevant in the given model.

The variable importances allow us to investigate if independent variables can be grouped with respect to the pattern in which they are relevant in various models.
It turned out, that amount of variables is too large to produce fully readable charts. Fortunately one can still read the outstanding variables, which gives us information which variables are the most controversial, in the sense that their relevance is the most non-uniform across the models. 



## LASSO regression

The final model we have chosen was _glmnet_ - an efficient implementation of the LASSO regression (https://arxiv.org/abs/1311.6529). LASSO linear regression (https://sci-hub.se/10.2307/2346178) is a linear regression with a single modification that adds to the objective function, (which is the residual sum of squares) a sum of absolute value of all regression coefficients, multiplied by the tunable parameter lambda. The modification acts as a form of penalty that discourages large models, i.e. models with a lot of non-zero and large regression coefficients. (Adding a penalty for large models is also called *regularization*).

The LASSO regression is defined as finding a vector of coefficients $\mathbf{a}$ that defines a regression function $$y = \mathbf{a} \mathbf{x}^T + x_0$$ such as it minimizes a LASSO cost function that is a modified square error of prediction, i.e. $$\mathbf{a} = \operatorname{argmin} \sum^N{\frac{\left( y - \hat{y} \right)^2}{N}} + \lambda \sum{ \left| \mathbf{a} \right|}$$ where $\lambda > 0$. Only the last term is specific to LASSO regression and adds additional term to the minimized expression. Parameter $\lambda$ is fitted to minimize the square error of prediction on validation data (this time without the LASSO term).

We choose the LASSO regression for its property, that minimization of its objective function automatically leads to variable selection. It is also fast and very extensively used in the machine learning and, on much more complicated models, in deep learning.
On the other hand the LASSO regression does not allow for close-form significance measures of regression coefficients. This lack of significance measures is not a big problem though, because we use statistical methods (described below) that bypass the problem and have an added benefit that they are known to be insensitive to regression’s restrictive assumptions, especially to the assumption of normally-distributed residuals (maybe citation?).
cross-validation

LASSO regression has one tuning parameter, $\lambda$, that is used to influence parsimony of the model. It cannot be calculated from the training sample, but like any other tuning parameter, it can be found by optimizing the RMSE (not the LASSO’s objective function) on the independent dataset (usually called “test dataset”).
Splitting the dataset into a “train” and “test” dataset is both wasteful (because we underuse the information in the test dataset) and arbitrary. To minimize the impact of those two problems we perform a procedure called 10-fold cross-validation. To minimize the bias further, we performed each 10-fold cross-validation ten times and averaged the results.


