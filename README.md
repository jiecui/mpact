Matching Pursuit based Adatpive Chirplet Transform - MPACT
==========================================================
**MPACT** is an open-source MatLab(R) toolbox that implements multi-component chirplet signal analysis 

<!-- $$
s\left( t\right)
	 = {\left( \frac{\alpha }{\pi }\right)} ^{1/4}\exp \left[ -\frac{\alpha }{2}{\left( t-t_{0}\right)} ^{2}\right]
	 \times \exp \left[ j\frac{\beta }{2}{\left( t-t_{0}\right)} ^{2}\right] \exp \left[ j\omega _{0}\left( t-t_{0}\right) \right]
$$ --> 

<div align="center"><img style="background: white;" src="https://render.githubusercontent.com/render/math?math=s%5Cleft(%20t%5Cright)%0A%09%20%3D%20%7B%5Cleft(%20%5Cfrac%7B%5Calpha%20%7D%7B%5Cpi%20%7D%5Cright)%7D%20%5E%7B1%2F4%7D%5Cexp%20%5Cleft%5B%20-%5Cfrac%7B%5Calpha%20%7D%7B2%7D%7B%5Cleft(%20t-t_%7B0%7D%5Cright)%7D%20%5E%7B2%7D%5Cright%5D%0A%09%20%5Ctimes%20%5Cexp%20%5Cleft%5B%20j%5Cfrac%7B%5Cbeta%20%7D%7B2%7D%7B%5Cleft(%20t-t_%7B0%7D%5Cright)%7D%20%5E%7B2%7D%5Cright%5D%20%5Cexp%20%5Cleft%5B%20j%5Comega%20_%7B0%7D%5Cleft(%20t-t_%7B0%7D%5Cright)%20%5Cright%5D"></div>

using matching pursuit algorithm.

The code repository for **MPACT** is hosted on GitHub at https://github.com/jiecui/mpact.

Installation
------------
1. Copy the files to the directory of your choice, e.g. ~/mpact/
1. Within Matlab, go the directory where you have copied the files
   e.g. >> cd(‘~/mpact/’);
1. Add the directory and its subdirectories into MatLab search path.
1. Type MatLab command, i.e. >> doc, to bring up MatLab help browser.
1. Find “Adaptive Chirplet Transform with Matching Pursuit Toolbox” in
   the section of Supplemental Software.
1. Follow the instructions on the screen.

References
----------
* Cui, J. and D. Wang (2017). "Biosignal Analysis with Matching-Pursuit Based Adaptive Chirplet." arXive (pre-print https://arxiv.org/abs/1709.08328)

* Cui, J. and W. Wong (2006). "The adaptive chirplet transform and visual evoked potentials." IEEE Transactions on Biomedical Engineering 53(7): 1378-1384.

* Cui, J., Wong, W., & Mann, S. (2005). "Time-frequency analysis of visual evoked potentials using chirplet transform." Electronics Letters, 41(4), 217-218. (PDF http://individual.utoronto.ca/jiecui/_private/Time_frequency_chirplet_vep.pdf)

License
-------
**MPACT** is distributed by the GPL v3 Open Source License.
