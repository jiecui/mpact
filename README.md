Matching Pursuit Based Adaptive Chirplet Transform (MPACT)
==========================================================

**MPACT** is an open-source MATLAB® toolbox that implements multi-component
chirplet signal analysis

$$s\left( t\right) = {\left( \frac{\alpha }{\pi }\right)} ^{1/4}\exp \left[
-\frac{\alpha }{2}{\left( t-t_{0}\right)} ^{2}\right] \times \exp \left[
j\frac{\beta }{2}{\left( t-t_{0}\right)} ^{2}\right] \exp \left[ j\omega
_{0}\left( t-t_{0}\right) \right]$$

using matching pursuit algorithm.  The code repository for **MPACT** is
hosted on GitHub at https://github.com/jiecui/mpact.

Installation
------------

1. Copy the files to the directory of your choice, e.g. ~/mpact/
1. Within MATLAB, go the directory where you have copied the files
   e.g. >> cd(‘~/mpact/’);
1. Add the directory and its subdirectories into MATLAB search path.
1. Type MATLAB command, i.e. >> doc, to bring up MATLAB help browser.
1. Find “Adaptive Chirplet Transform with Matching Pursuit Toolbox” in
   the section of Supplemental Software.
1. Follow the instructions on the screen.

References
----------
* __Cui, J.__ and D. Wang (2017). "Biosignal Analysis with Matching-Pursuit
  Based Adaptive Chirplet." arXive (pre-print
  https://arxiv.org/abs/1709.08328)

* __Cui, J.__ and W. Wong (2006). "The adaptive chirplet transform and
  visual evoked potentials." IEEE Transactions on Biomedical Engineering
  53(7): 1378-1384.

* __Cui, J.__, Wong, W., & Mann, S. (2005). "Time-frequency analysis of
  visual evoked potentials using chirplet transform." Electronics Letters,
  41(4), 217-218. (PDF
  http://individual.utoronto.ca/jiecui/_private/Time_frequency_chirplet_vep.pdf)

License
-------
**MPACT** is distributed by the GPL v3 Open Source License.
