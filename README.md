# adapt-discr-efficient-code
Data from orientation adaptation 2AFC experiment and efficient coding model.

Mao, J., Rothkopf, C. A., & Stocker, A. A. Adaptation optimizes sensory encoding of future stimuli. 

## Dependencies
* Philipp Berens (2022). Circular Statistics Toolbox (Directional Statistics) (https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics), MATLAB Central File Exchange. 
* We have included the function needed (circ_vmpdf.m), but would like to explicitly acknowledge the use of the package.

## data
Data from adaptation experiment using a 2AFC task measuring discrimination threshold of visual orientation after adapting to an oblique (45 or 22.5 deg) adaptor or a control adaptor (all orientations).

* test: test orientation relative to adaptor.
* dtheta: difference between comparison and test stimuli.
* nTot: number of trials completed for each dtheta and each test.
* nRight: number of trials where subject answered that comparison was more clockwise than test stimulus, for each dtheta and each test.

## analyze
* bootstrap\_psychometric.m: fit psychometric curve to original data, bootstrap, and fit psychometric curve to bootstraped data.

## model
Model simulations and fits use parallel computing.

* fit\_ctrl\_par: higher coding accuracy at cardinal orientations; fit to 45 & 22.5 control conditions jointly.
* fit\_2peak\_par: (overall best model) increased coding accuracy at adaptor & 90deg off; fit to 45 & 22.5 oblique conditions jointly.
* fit\_1peak\_par: increased coding accuracy at adaptor only; fit to 45 & 22.5 oblique conditions jointly.
* fit\_2peakF\_par: increased coding accuracy at adaptor & 90deg off, allow different total Fisher from control condition; fit to 45 & 22.5 oblique conditions jointly.
* fit\_2peakK\_par: increased coding accuracy at adaptor & 90deg off, allow different kernels for 45 & 22.5 adaptors; fit to 45 & 22.5 oblique conditions separatly.

## plot
* threshold\_data\_fit\_avesub.m - Fig.3: discrimination threshold data and model fits averaged across subjects.
* threshold\_data\_fit\_indivsub.m - Fig.4: discrimination threshold data and model fits for individual subjects.
* kernel\_fit\_indivsub.m - Fig.5a: adaptation kernels for individual subjects.
* totFisher\_BIC\_indivsub.m - Fig.5b&c: total Fisher in control v.s. oblique condition; BIC of different models.
