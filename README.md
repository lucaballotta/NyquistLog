# NyquistLog

NyquistLog plots the Nyquist diagram of continuous-time causal stable LTI systems, where the real and imaginary parts of the frequency response are in logarithmic scale, allowing for better visualization and interpretation of the system response. Nyquist stability criterion and visual methods - e.g., phase margin and vector margin - hold like in the standard polar plot.\

Also, this function is able to handle both poles at the origin and pairs of complex conjugates poles on the imaginary axis, by means of corresponding semi-circles "at infinity".
Further, it performs stability analysis via the Nyquist stability criterion and prints the number of unstable poles of the closed-loop system, together with the number of net encirclements of the Nyquist diagram around the -1 point.\

Any comments or bug reports will be highly appreciated.\

References:\
Trond Andresen (2021). Nyquist plot with logarithmic amplitudes (https://www.mathworks.com/matlabcentral/fileexchange/7444-nyquist-plot-with-logarithmic-amplitudes), MATLAB Central File Exchange. Retrieved June 3, 2021.\
Federica Grossi (2021). Closed Logarithmic Nyquist plot (https://www.mathworks.com/matlabcentral/fileexchange/43768-closed-logarithmic-nyquist-plot), MATLAB Central File Exchange. Retrieved June 3, 2021.
