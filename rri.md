# Meetings

## Check Possibilities


-   Representative papers
    -   **Master:** [PRL](file:///home/paulo/Documents/Mendeley Desktop/Galuzio, Lopes, Viana - 2010 - Two-State On-Off Intermittency and the Onset of Turbulence in a Spatiotemporally Chaotic System.pdf)
    -   **PhD:** [PRERC](file:///home/paulo/Documents/Mendeley Desktop/Galuzio, Viana, Lopes - 2014 - Control of extreme events in the bubbling onset of wave turbulence.pdf)
    -   **PostDoc:** [COBEM](file:///home/paulo/Documents/PostDoc/COBEM2019/Template_COBEM2019_Latex_Final_Paper/cobem2019.pdf)
    -   **RP:** [Mammography](file:///home/paulo/Documents/Mendeley Desktop/Prado et al. - 2011 - Spatial recurrence analysis a sensitive and fast detection tool in digital mammography.pdf)
    -   **Hamiltonian:** [Stickiness](file:///home/paulo/Documents/Mendeley Desktop/Krüger et al. - 2015 - Mechanism for stickiness suppression during extreme events in Hamiltonian systems.pdf)
-   CV [File](file:///home/paulo/Documents/Curriculum/paulo/cv.pdf)

## New Meeting w/ Alhaji


-   gonna work with the preprocessing idea we first talked about

### DONE One Drive

-   enviar as informações do one drive para ele

## First meeting w/ Peter


-   Honorariums USD 5k for first deadline + 2 X USD 3K

# Research

## PreProcessing

:GeneralInfo:
-   Term in statistics for gap filling: **imputation**

:END:

### Ideas

-   Attractor reconstruction


    1.  Use GP (or spline) of first variable to reconstruct attractor.
    2.  Optimize delays and other parameters to best fit the other variables
    3.  Output is spline or GP of reconstructed series

-   TODO Recurrence plot

    -   RQA is implemented through pyRQA [Usage](https://pypi.org/project/PyRQA/#usage)
    -   [ ] Maximize determinism in cross recurrence plot between series

-   TODO Singular spectrum analysis (SSA)

-   Gaussian Processes

    -   [ ] Kernel choices
        -   Comnination of fourier decomposition
        -   Some sort of broadband kernel, where one can provide the exponent for fourier spctrum exponential fit [Possible Reference](https://aaltodoc.aalto.fi/handle/123456789/37143)

-   mrDMD


    -   Use all the points to go through a DMD slow moving modes
    -   Use data with better resolution to Iterate through faster modes

-   TODO Topics to go through

    -   [-] Use of Neural Networks
        -   [X] [PRL Iten 2020](file:///home/paulo/Documents/RRI/References/papers/time_series/neural networks/PhysRevLett.124.010508.pdf) and [Suplementary information](file:///home/paulo/Documents/RRI/References/papers/time_series/neural networks/Supplementary_information.pdf)
                 Very good paper, but without an interest in modeling the equations but rather to only use the neural network to build a low dimensional representation of the system.
        -   [X] [Nature Daniels 2015](file:///home/paulo/Documents/RRI/References/papers/time_series/neural networks/daniels2015.pdf)   [Code](file:///home/paulo/Documents/RRI/References/Codes/SirIsaac)
            -   They created a code that does model selection, although there are not much details of how they do it, they do provide their code
            -   [ ] More details are provided in the [PlosOne paper](file:///home/paulo/Documents/RRI/References/papers/time_series/neural networks/daniels2015plos.pdf)
        -   [ ] [PNAS brunton 2016](file:///home/paulo/Documents/RRI/References/papers/time_series/neural networks/brunton2016.pdf)
        -   [ ] [arXiv Lusch 2018](file:///home/paulo/Documents/RRI/References/papers/time_series/neural networks/1712.09707.pdf)
        -   [ ] [arXiv Takeishi 2018](file:///home/paulo/Documents/RRI/References/papers/time_series/neural networks/1710.04340.pdf)
        -   [ ] [arXiv Otto 2019](file:///home/paulo/Documents/RRI/References/papers/time_series/neural networks/1712.01378.pdf)
        -   [ ] [arXiv Raissi 2018](file:///home/paulo/Documents/RRI/References/papers/time_series/neural networks/1801.06637.pdf)
        -   [ ] [arXiv Zhang 2019](file:///home/paulo/Documents/RRI/References/papers/time_series/neural networks/1905.01205.pdf)
        -   [ ] [J Comp Phys Raissi 2018](file:///home/paulo/Documents/RRI/References/papers/time_series/neural networks/raissi2018.pdf)
        -   [ ] [J Comp Phys Raissi 2019](file:///home/paulo/Documents/RRI/References/papers/time_series/neural networks/raissi2019.pdf)
    -   [ ] singular spectrum analysis (SSA)
        -   **Notes:** PCA based technic, in fact: **SSA is a PCA applied to the embeded time series**
        -   **References**
            -   [X] [wiki](https://en.wikipedia.org/wiki/Singular_spectrum_analysis)
                -   **notes:** just took an overall look
            -   [ ] [Kondrashov2010](file:///home/paulo/Documents/RRI/References/papers/time_series/singular_spectrum_analysis/2010GL044138.pdf) **ON IT**
                -   **notes:** gap filling using SSA
            -   [X] [tese](file:///home/paulo/Documents/RRI/References/papers/time_series/singular_spectrum_analysis/24787_4.PDF)
                -   **Notes:** good reference for detailed explanantion of basic SSA, however maybe it will be necessary a more complex variation of the method to implement in our problem correctly
            -   [ ] [Vautard1992](file:///home/paulo/Documents/RRI/References/papers/time_series/singular_spectrum_analysis/vautard1992.pdf)
            -   [ ] [MPRA 2007](file:///home/paulo/Documents/RRI/References/papers/time_series/singular_spectrum_analysis/MPRA_paper_4991.pdf)
        -   **Ideas**
            -   `1st thoughts`
                1.  Make SSA on time instants were all three components are known
                2.  Write better time series in terms of combination of PCAS
                3.  Use same components to describe other series
                4.  maybe do some iteration to minimize errors in all known points
    -   [ ] imputation
        -   **references:** -   [ ] [wikipedia](https://en.wikipedia.org/wiki/Imputation_(statistics))
            -   [ ] [Barnard 1999](file:///home/paulo/Documents/RRI/References/papers/time_series/multiple_imputation/barnard1999.pdf)
            -   [ ] [Horton 2007](file:///home/paulo/Documents/RRI/References/papers/time_series/multiple_imputation/horton2007.pdf)
    -   [ ] Gap Filling in Time Series&#x2026; by Bulhões 2017 [File](file:///home/paulo/Documents/RRI/References/papers/time_series/bulhoes2017.pdf)
        -   spectral analysis on time series
        -   application of Hammerstein-Wiener models
    -   [ ] Kohonen self-organizing maps
        -   [ ] [Dergachev 2001 Cambridge](file:///home/paulo/Documents/RRI/References/papers/time_series/Self-Organizing maps/filling_of_gaps_in_geophysical_time_series_by_artificial_neural_networks.pdf)
        -   [ ] [Dergachev 2002](file:///home/paulo/Documents/RRI/References/papers/time_series/Self-Organizing maps/gai01384.pdf)
        -   [ ] [Dergachev 2001](file:///home/paulo/Documents/RRI/References/papers/time_series/Self-Organizing maps/geo_20_07.pdf)

    -   [ ] Sergey V Skakun and Ruslan M Basarab. Reconstruction of missing data in time-series of optical satellite images using self-organizing kohonen maps. Journal of Automation and Information Sciences, 46(12), 2014.

    -   [X] ROCHTU Yannick, Jean-Baptiste SAUBIN, and Jean-Luc BERTRAND-KRAJEWSKI. Filling gaps in time series in urban hydrology. Laboratoire de Gnie Civil et d’Ingnierie Environnementale (LGCIE), 2014. [File](file:///home/paulo/Documents/RRI/References/papers/time_series/Self-Organizing maps/RUG01-002153908_2014_0001_AC.pdf)
        -   The paper reviews the literature, and the Kohonen self-organizing maps appears in the review, maybe the method proposed is something else.

        -   General description of the method
            1.  K-means clustering

                **Poorly written English**, very hard to understand

    -   [ ] multiple imputation method
        -   Jeffrey C Wayman. Multiple imputation for missing data: What is it and how can i use it. In Annual Meeting of the American Educational Research Association, Chicago, IL, pages 2–16, 2003.
        -   John W Graham, Allison E Olchowski, and Tamika D Gilreath. How many imputations are really needed? some practical clarifications of multiple imputation theory. Prevention Science, 8(3):206–213, 2007.
    -   [ ] multiple regression analysis
        -   Marjolaine Métadier. Treatment and Analysis of Continuous Time Series of Turbidity for Formulation and Filetest of Model of Urban Discharges in Rainy Weather. Available in French. PhD thesis, PhD thesis, INSA-Lyon, France, 2011.
        -   Jon K Eischeid, C Bruce Baker, Thomas R Karl, and Henry F Diaz. The quality control of long-term climatological data using objective data analysis. Journal of Applied Meteorology, 34(12):2787–2795, 1995.
    -   [ ] kringing
    -   [ ] Hot Deck: a missing value was imputed from a randomly selected similar record
    -   [ ] Statistical tests:
        -   [ ] ljung-box
        -   [ ] Kendall test
    -   [-] **DMD**
        -   [ ] Koopman spectral decomposition
        -   [ ] koopman modes
        -   [-] Papers
            -   [X] [Kutz 2015](file:///home/paulo/Documents/RRI/References/papers/time_series/DMD/1506.00564.pdf)
                -   Linear approximation of the non-linear system: \(\frac{d \tilde x}{dt} = A \tilde x\)
                -   Eigenfunction expansion:

                    \[\tilde x(t) = \sum_{k=1}^{K} b_{k} \phi_{k} e^{\omega_{k t}} \]

                -   **\(N\):** number of spatial measurements per time snapshot: **NUMBER OF EQUATIONS**
                -   **\(M\):** number of snapshots taken in time
                -   Combination of Proper orthogonal decomposition (POD) with fourier transforms in time
                -   Data snapshots in \(N\times M\) matrix:

                    \[\mathbb{X} = \left[\vec x (t_{1})\; \vec x (t_{2})\; \vec x (t_{3}) \cdots \vec x (t_{M})\right]\]

                -   **Definition:** \(\mathbb{X}_{j}^{k} = [\vec x (t_{j}) \cdots \vec x (t_{k})]\)

                -   The evolution of the states can be approximated by:

                    \[\mathbb{X}_{2}^{M} \approx A \mathbb{X}_{1}^{M-1}\]

                    The operator \(A\) advances each snapshot column a single time step \(\Delta t\).

                -   **DMD algorithm**
                    1.  Take the SVD of \(\mathbb{X}_{1}^{M-1}\):

                        \[\mathbb{X}_{1}^{M-1}= \mathbb{U}\Sigma\mathbb{V}^{*}\]

                        \(*\) is conjugate transpose; \(\mathbb{U} \in C^{N\times K}\), \(\Sigma \in C^{K\times K}\), \(\mathbb{V} \in C^{M-1 \times K}\), \(K\) is the rank of the reduced SVD approximation.

                    2.  Compute \(\tilde A\), the \(K\times K\) projection of the full matrix \(A\) onto POD modes:

                        \[A = \mathbb{X}_{2}^{M} (\mathbb{X}_{1}^{M-1})^{-1} = \mathbb{X}_{2}^{M} \mathbb{V} \Sigma^{-1} \mathbb{U}^{*}\]

                        \[ \tilde A = \mathbb{U}^{*} A \mathbb{U} = \mathbb{U}^{*} \mathbb{X}_{2}^{M} \mathbb{V} \Sigma^{-1}\]

                    3.  Compute the eingendecomposition of \(\tilde A\):

                        \[\tilde A \mathbb{W} = \Lambda \mathbb{W}\]

                        \(\Lambda\) is a diagonal matrix with eigenvalues \(\lambda_{k}\)

                    4.  Reconstruct the eigendecomposition of \(A\), where \(\Psi\) are the eigenvalues of \(A\)

                        \[ \Psi = \mathbb{X}_{2}^{M} \mathbb{V} \Sigma^{-1} \mathbb{W}  \]

                    5.  The approximate solution at future times is:

                        \[\vec x_{DMD} (t) = \sum_{k=1}^{K} b_{k}(0) \psi_{k}(\vec \xi) e^{\omega_{k} t}\]

                        \(\omega_{k} = \ln (\lambda_{k})/\Delta t\); \(\vec \xi\) are the spatial coordinates (?); \(b_{k}(0)\) is the initial amplitude of each mode.

                    6.  \(\vec b = \Psi^{+} \vec x_{1}\)

                -   The value of \(M\) is chosen so that a full rank approximation can be acomplished
                -   The multi-resolution regards the separation of the DMD eigenvalue expansion is "slow" and "fast" modes, and then reapply DMD on the reconstructed series with fast modes only and proceed with the separation
            -   [ ] Equation free multi-scale modelling
                -   [ ] [kevrekidis 2003](file:///home/paulo/Documents/RRI/References/papers/time_series/DMD/euclid.cms.1119655353.pdf)
                -   [ ] [kevrekidis 2004](file:///home/paulo/Documents/RRI/References/papers/time_series/DMD/Kevrekidis-EqFree-AIChE04.pdf)
            -   [X] [Sparse spatial and temporal sampling (by Kutz 2014)](file:///home/paulo/Documents/RRI/References/papers/time_series/DMD/1401.7047.pdf)
                -   **Nyquist–Shannon sampling theorem:** It establishes a sufficient condition for a sample rate that permits a discrete sequence of samples to capture all the information from a continuous-time signal of finite bandwidth.

                    `If a function \(x ( t )\)  contains no frequencies higher than \(B\) hertz, it is completely determined by giving its ordinates at a series of points spaced \( 1 / ( 2 B )\) seconds apart. A sufficient sample-rate is therefore anything larger than \(2 B\)  samples per second.`

                -   **numerical method:** scalar signals
                    -   **signal:** \(f \in \mathbb{R}^{n}\), typically not sparse in the standard basis for \(\mathbb{R}^{n}\)

                    -   **compressibility:** \(f\) is compressible if there is a basis \(\psi\) such that the representation of \(f\) in \(\psi\) is approximately sparse. \(f\) is **k-sparse** in the basis \(\psi\) if:

                        \[f=\Psi \hat  f\]

                        where \(\psi \in \mathbb{R}^{n\times n}\) and \(\hat f \in \mathbb{R}^{n}\) with only \(k \ll n\).

                    -   all we know is an m-dimensional linear measurement:

                        \[g=\Phi^{T}f\]    (2)

                        \(\Phi\) is an \(n\times m\) matrix. when \(m\ll n\) then eq. (2) is under-determined and we cannot solve for \(f\) from a knowledge of \(g\); substituting \(f\):

                        \[g=\Phi^{T} \Psi \hat f\]     (3) &#x2013;> assuming \(k<m\ll n\)

                    -   In this work **\(\Psi\) is chosen to be the DFT basis**

                            import scipy as sp

                            def dftmtx(N):
                                return sp.fft(sp.eye(N))

                            Psi = dftmtx(N)

                    -   **Optimization problem:** to find \(\hat f\)

                        \[\min_{\hat f \in \mathbb{R}^{n}} ||\hat f ||_{1}\] subject to \[g=\Phi^{T} \Psi \hat f\]

                    -   \(g\) is composed by values of \(f\) sampled at particular instants of time.

                    -   If \(g\) is sampled at \(t=0\) and \(t=2\Delta t\), for example, then the first column of \(\Phi\) is \((1,0,0,\ldots)^{T}\) and the second column is \((0,0,1,0,\ldots)^{T}\). Each column of \(\Phi\) contains a subset of the standard basis.

                    -   **Vector valued signals:** \(F\in \mathbb{R}^{n\times p}\) where \(p\) is the number of measures in a given instant of time

                    -   Now the optimization problem becomes:

                        \[\min_{\hat F \in \mathbb{R}^{n\times p}} ||\hat F ||_{1,q}\] subject to \(G=\Phi^{T}\Psi \hat F\)

                        Where the mixed norm is defined as:

                        \[||M||_{1,q} = \sum_{i=0}^{n-1} \left|\left(\sum_{j=0}^{p-1}|M_{i,j}|^{q}\right)^{\frac{1}{q}} \right|\]

                            import numpy as np

                            Norm = np.linalg.norm([np.linalg.norm(y,q) for y in x],1)

                    -   **The choice of \(q\) penalizes nonzero elements within rows:** a way of seeing it is as if \(q\) somehow increased the coupling between the different \(p\) components of your vector. If \(q=1\) we have the 1-norm and the components are uncoupled.

            -   [X] [Irregularly spaced times](file:///home/paulo/Documents/RRI/References/papers/time_series/DMD/1312.0041.pdf) (**partial**)
                -   **Review:** I. Mezic. Analysis of fluid flows via spectral properties of the Koopman operator. Annu. Rev. Fluid Mech., 45:357{378, 2013

                -   **Alogrithm:** Sequential set of data vectors \(\{z_{0}, z_{1},\ldots,z_{m}\}\), \(z_{k}\in \mathbb{R}^{n}\)

                    \[z_{k+1}=Az_{k}\]

                    -   In the case of continuous system:  \(z_{k}=z(k\Delta t)\), and a fixed sampling time is assumed.

                    -   the algorithm is just the same as the standard DMD algorithm, except that now the data is presented by pairs of data \(\{(x_{k}, y_{k})\}\), and the data matrices are constructed as:

                        \[X=\left[ x_{1} \; x_{2}\; \cdots x_{m}\right]\]

                        \[Y=\left[ y_{1} \; y_{2}\; \cdots y_{m}\right]\]

                        and \(Y=AX\)

                    -   **Definition 2:** (linear consistency) two \(n\times m\) matrices are linearly consistent if, whenever \(Xc=0\), then \(Yc=0\) as well. Which is satisfied trivially in the case that the columns of \(X\) and \(Y\) are linearly independent.

                    -   **Theorem:** Define \(A=YX^{+}\) then \(Y=AX\) if and only if \(Y\) and \(X\) are linearly consistent.

                    -   **OBS:** the reference [pyrunner](http://www.pyrunner.com/weblog/2016/07/25/dmd-python/) implement the method in this reference
            -   [ ] [Review SIAM 2017](file:///home/paulo/Documents/RRI/References/papers/time_series/DMD/Hankel-DMD paper.pdf)

### Coding

-   DONE Data Generation


    -   [X] Implement Lorentz oscilator
        -   File: [lorentz.py](file:///home/paulo/Documents/RRI/Codes/Lorentz/lorentz.py)
        -   [Wiki](https://en.wikipedia.org/wiki/Lorenz_system)
        -   Equations:
            xprime = sigma \* (y - x)
            yprime = x \* (rho - z) - y
            zprime = x \* y - beta \* z
        -   **Parameters:** sigma = 10
            beta = 8/3
            rho = 28
        -   Output as pandas dataframe
            -   Define sampling for each variable, for every non-sampled time replace wit np.nan
    -   [X] Documentation
    -   [X] Verify sampling if it is done appropriately

-   Data Analysis



    -   General properties of every module

        -   <span class="underline">Inputs</span>
            -   **t:** array with times when to sample the series
            -   **df:** dataframe with the time series
        -   <span class="underline">Outputs</span>
            -   **dfRight:** dataframe with time series sampled at **t** moments
        -   [X] Checking input for analysis in [Analysis](file:///home/paulo/Documents/RRI/Codes/Analysis)

    -   DONE Spline


        -   [X] Implement spline

    -   Attractor reconstruction

    -   VERIFY Gaussian Process



        -   DONE Kernel

            -   Implement different kernel selection
                -   **RBF:** The RBF kernel is a stationary kernel. It is also known as the “squared exponential” kernel. It is parameterized by a length-scale parameter length<sub>scale</sub>>0, This kernel is infinitely differentiable. **parameters:** length<sub>scale</sub>, amplitude
            -   **ExpSineSquared:** The ExpSineSquared kernel allows modeling periodic functions. It is parameterized by a length-scale parameter length<sub>scale</sub>>0 and a periodicity parameter periodicity>0
            -   **parameters:** L: length scale, same for all kernels?
                RBF<sub>amp</sub>: amplitude for RBF kernel
                periodicity: (`how many?`)
                Sine<sub>amp</sub>: amplitude for periodic kernel (`how many?`)

        -   DONE Kernel sum

            -   [X] verify wether or not the amplitudes for summing different kernels has any impact on the final outcome
            -   [X] verify how parameter optimization take place in GP

        -   TODO Different kernel

            -   [ ] apply different kernel parameter per each variable

    -   PARTIAL DMD


        -   [X] Pip install of [pyDMD](https://mathlab.github.io/PyDMD)
        -   [ ] check implementation step by step side by side w/ theory
        -   [ ] verify why it doesn't work on lorenz system, maybe test on some of Kutz data

        -   References

            -   [ ] [pyrunner](http://www.pyrunner.com/weblog/2016/07/25/dmd-python/)
            -   [ ] <https://github.com/arbabiha/DMD-for-ergodic-systems>
            -   [ ] [Arbabi2017](file:///home/paulo/Documents/RRI/References/papers/time_series/DMD/1611.06664.pdf)
            -   [ ] <file:///home/paulo/Documents/RRI/References/papers/time_series/DMD/1206.3164.pdf>
            -   [X] **[Introdução](file:///home/paulo/Documents/RRI/References/papers/time_series/DMD/KoopmanIntro.pdf):** by the end there are several references of **variations of DMD**
            -   [ ] **[Mezic2005](file:///home/paulo/Documents/RRI/References/papers/time_series/DMD/art3a10.10072fs11071-005-2824-x.pdf):** Both papers also discussed the idea of applying Koopman methodology to capture
                the regular components of data in systems with combination of chaotic and regular
                behavior.
            -   [ ] [Rowley2009](file:///home/paulo/Documents/RRI/References/papers/time_series/DMD/10.1.1.707.7377.pdf) pp. 313 or 323 showed that KMD can be computed through a numerical decomposition technique known as Dynamic Mode Decomposition (DMD)

    -   TODO Sparse space and time analysis (Compressed DMD)


        -   General Algorithm
            1.  [X] Build G, the measurement vector
            2.  [X] Build \(\Phi\) with the sampling times of G
            3.  [X] Build \(\Psi\) as the DFT matrix
            4.  [X] Implement mixed norm
            5.  [ ] Implement POD to reduce system dimensionality
            6.  [X] Implement the constrained optimization problem :: the constraint was replaced by:
                \[\left|\left|G-\Phi^{T}\Psi\hat F\right|\right|_{1,2}=0\]
            7.  [ ] Approximate the signal by a GP using \(\hat F\) as an approximation to the Fourier modes
            8.  [ ] Resample the signal using few points from the GP approximation and repeat the process.

        -   topics to fix
            1.  [ ] Minimize the constraint, that is the best periodic approximation
            2.  [ ] if the minimization is successful, minimize the norm
            3.  [ ] verify values
            4.  [ ] write documentation
            5.  [ ] check all zero values
            6.  [ ] check how the inverse is being calculated

        -   **First Results:** shitty shit (all zero)

    -   POSTPONED mrDMD


        1.  Calculate

            \[A = X_{2}^{M}(X_{1}^{M-1})^{-1}\]
        2.  Eigendecomposition of \(A\):

            \[AW = \Lambda W\]
        3.  Initial condition

            \[\vec c = W^{-1}\vec x (t=0)\]
        4.  Build approximation

            \[\vec x^{(1)} (t) = \sum_{k=1}^{M_{3}} c_{k} e^{\omega_{k} t} \vec w_{k}\]    (1)

            where \(\omega_{k}=\frac{\ln(\lambda_{k})}{\Delta t}\)
        5.  Separate the series in \(M_{3}\) windows
        6.  Use (1) to approximate \(x_{3}\)
        7.  Repeat process in windows

        8.  Implement process in smaller scales using sub windows

    -   VERIFY GP and DMD


        1.  Just addapting from GP adn Spline

    -   VERIFY GP and Spline


        1.  Make a GP approximation and sample it at pre-selected points
        2.  Use the points from the GP to do spline
        3.  Use spline to predict the known points
        4.  calculate R2
        5.  optimeze the previous process for maximum R2
        6.  repeat the process with optimal values for the GP

        7.  [X] document the module
        8.  [X] [verify of selection of non null elements is correct](file:///home/paulo/Documents/RRI/Codes/Analysis/GPandSpline.py)
        9.  [X] [Implement number of frequencies as parameter](file:///home/paulo/Documents/RRI/Codes/Analysis/GPandSpline.py)
        10. [X] verify error in execution
        11. [X] Check argument passing
        12. [X] Implement normalized squared error
        13. [X] force frequencies to be different through inequality constraint abs(x[i]-x[j]) > 0 for all combinations of x[i] and x[j], [Correct minimal difference between frequencies](file:///home/paulo/Documents/RRI/Codes/Analysis/GPandSpline.py)
        14. [ ] Implement weights proportional to sampling [Here](file:///home/paulo/Documents/RRI/Codes/Analysis/GPandSpline.py)
        15. [X] correct constraints and make sure they are being satisfied during optimization [Reference](https://prog.world/scipy-conditions-optimization/)
        16. [X] Implement my Bayesian optimization
        17. [X] [check allowed values for nFreq](file:///home/paulo/Documents/RRI/Codes/Analysis/GPandSpline.py)

        18. [ ] Sometimes the code is breaking, need a little debugging

-   VERIFY Bayesian Optimization


    -   [ ] document the module
    -   [ ] check why it stops sometimes
    -   [X] Implement RBF kernel

-   TODO Tools <code>[3/4]</code>


    -   [X] Module to plot dataframes
    -   [X] Scikit learner mean<sub>squared</sub><sub>error</sub>
    -   [X] Increment function to plot dataframes together
    -   [ ] Metrics

## Changes track


:FileTree:
.
├── Analysis
│   ├── checkInput.py
│   ├── CompressedDMD.py
│   ├── dmd.py
│   ├── gp.py
│   ├── GPandDMD.py
│   ├── GPandspline.py
│   ├── mrdmd.py
│   ├── <span class="underline"><span class="underline">init</span></span>.py
│   └── spline.py
├── <span class="underline"><span class="underline">init</span></span>.py
├── Lorentz
│   ├── <span class="underline"><span class="underline">init</span></span>.py
│   └── lorentz.py
├── scripts
│   ├── pathmagic.py
│   ├── testCompressedDMD.py
│   ├── testDMD.py
│   ├── testFrequencies.py
│   ├── testmrDMD.py
│   ├── testSpline.py
│   ├── testeGPandDMD.py
│   ├── testeGPandSpline.py
│   └── testGP.py
└── Tools
    └── Plots.py
:END:
-   <span class="timestamp-wrapper"><span class="timestamp">&lt;2019-11-26 ter&gt;</span></span> Added small random perturbation to lorentz sampling
-   **<span class="timestamp-wrapper"><span class="timestamp">&lt;2019-11-26 ter&gt;</span></span> Implemented different kernels in GP:** kernel properties are now passed as arguments of the function
-   **<span class="timestamp-wrapper"><span class="timestamp">&lt;2019-11-27 qua&gt;</span></span> gp.py:** (i) change noise<sub>level</sub> and white kernel for alpha; (ii) Kernel amplitude through constant kernel; (iii) Include optimize boolean argument to module
-   **<span class="timestamp-wrapper"><span class="timestamp">&lt;2019-11-27 qua&gt;</span></span> dmd.py:** (i) created the file and started to implement the simplest possible implementation of it
-   **<span class="timestamp-wrapper"><span class="timestamp">&lt;2019-11-27 qua&gt;</span></span> GPandSPline.py:** (i) created the file and started implementation; (ii) started first tests and debug
-   **<span class="timestamp-wrapper"><span class="timestamp">&lt;2019-12-02 seg&gt;</span></span> GPandSpline.py:** (i) corrected error in execution (size of t array to spline module); (ii) implemented number of frequencies as argument of the module; (iii) implemented bounds as arguments; (iv) implemented a difference between frequencies possible; (v) save time series got from GP as well; (vi) implemented bayesian
-   **<span class="timestamp-wrapper"><span class="timestamp">&lt;2019-12-02 seg&gt;</span></span> Plots.py:** (i) Added optional title to plots; (ii) Plot true series with data
-   **<span class="timestamp-wrapper"><span class="timestamp">&lt;2019-12-03 ter&gt;</span></span> GPandspline.py:** (i) comments and squared error
-   **<span class="timestamp-wrapper"><span class="timestamp">&lt;2019-12-03 ter&gt;</span></span> pathmagic.py:** (i) comply with pep and include path in sys
-   **<span class="timestamp-wrapper"><span class="timestamp">&lt;2019-12-03 ter&gt;</span></span> dmd.py:** (i) implemented a simple version based on Pyrunner &#x2013; not working, need to go by detail
-   **<span class="timestamp-wrapper"><span class="timestamp">&lt;2019-12-03 ter&gt;</span></span> testDMD.py:** simple script to run dmd approximation
-   **<span class="timestamp-wrapper"><span class="timestamp">&lt;2019-12-16 seg&gt;</span></span> mrdmd.py:** created the file; (i) implemented simpler possible version of DMD without SVD using only coarse sampling of the time series
-   **<span class="timestamp-wrapper"><span class="timestamp">&lt;2019-12-16 seg&gt;</span></span> testmrDMD.py:** created the file as copy from testDMD.py
-   **<span class="timestamp-wrapper"><span class="timestamp">&lt;2019-12-16 seg&gt;</span></span> dmd.py:** reimplemented it based on the kutz algorithm [DMD.m](file:///home/paulo/Documents/RRI/References/Codes/KutzDMD/CH01_INTRO/DMD.m)
-   **<span class="timestamp-wrapper"><span class="timestamp">&lt;2019-12-19 qui&gt;</span></span> CompressedDMD.py:** created the file
-   **<span class="timestamp-wrapper"><span class="timestamp">&lt;2019-12-19 qui&gt;</span></span> CompressedDMD.py:** created a class with the problem (i) try to find the best periodic approximation, uses only n<sub>modes</sub> fourier modes, and try to minimize the constraint. If the constraint ever gets to zero, solve the minimization problem
-   **<span class="timestamp-wrapper"><span class="timestamp">&lt;2019-12-19 qui&gt;</span></span> testCompressedDMD.py:** created the script based on the DMD one

## Papers

-   [ ] good review [Rochtus 2014](file:///home/paulo/Documents/RRI/References/papers/time_series/Self-Organizing maps/RUG01-002153908_2014_0001_AC.pdf)

# Studying (20% of time)

-   [X] Singular value decomposition [Youtube](https://www.youtube.com/watch?v=EokL7E6o1AE)
-   [X] Principal component analysis [Youtube](https://www.youtube.com/watch?v=a9jdQGybYmE)
-   [ ] Python decorators
-   [ ] Tortoise SVN
-   [ ] Python module importing [PEP](https://www.python.org/dev/peps/pep-0420/)
-   [ ] Take a look at Orthogonal forward regression
-   [ ] Armoldi algorithm

## READ Papers

-   [ ] [Nature Communications 2017](file:///home/paulo/Documents/RRI/References/papers/TakeALook/s41467-017-00030-8.pdf)
-   [ ] [Estimating Equations Inference with Missing Data](file:///home/paulo/Documents/RRI/References/papers/TakeALook/27640153.pdf)
-   [ ] [Hutz on DMD](file:///home/paulo/Documents/RRI/References/Books/J. Nathan Kutz, Steven L. Brunton, Bingni W. Brunton, Joshua L. Proctor - Dynamic Mode Decomposition_ Data-Driven Modeling of Complex Systems-SIAM-Society for Industrial and Applied Mathematics (2016).pdf)

# Reports

## Update New Year


General conditions in all simulations:
Day = 0.01
Step = min(Day \* 0.1, 0.01)
xstep = 2 \* Day
ystep = 7 \* Day
zstep = 120 \* Day
T = 365 \* Day

### Algorithms Tried so far

-   **Compressed DMD:** poor result
-
