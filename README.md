# 3D reconstruction from noisy 2D projections

Based on parallel beam geometry single axis tilting.

The 3D volume can be regarded as several slices (2D images) along one direction.

## Projection

The Parallel beam projection consists of a vector of parallel line integrals, we compute the projection by rotating several angles spacing between (0,180), with density $f(x)=f(x,y)$ at every angle $θ$ and displacement $r$.

Mathematically: the projection $P_θ(r)$ of an image can be described by the Radon Transform, is a function on the space of line integral along each line as, $P_θ (r) = \int_{l}^{} f (x) dx$. $s$ is the distance between the line and the origin. $θ$ is the angle between the x-axis and s-axis.

Basically, we take the image (which is just a matrix of intensities in MATLAB), rotate it, and sum up the intensities (the sum of each column).  

In MATLAB this is easily accomplished with the 'imrotate' and 'sum' commands.  

First, we zero pad the image so we don't lose anything when we rotate (the images are rectangular so the distance across the diagonal is longer than the distance on a side).  

Then we rotate the image $θ$ degrees (so that the projection is lined up in the columns) using the 'imrotate' command, and finally summed up the columns using the 'sum' command.

Gaussian noise or amplifier noise is a form of Gaussian distribution that is usually added during an image acquisition. The amplifier noise is caused by sensors, high temperature and electronic circuit components.



## Simple back projection

Back projection are just the summation of all projection according to the angle.

In the central area of the FFT image (called a **power spectrum**) corresponding to the low frequencies, the various spokes are close to each other and thus, there is no "empty space". However, at the periphery (the high frequencies area), the spokes are distant. In other words, we have much more missing data at the periphery than at the center. In conclusion, we have more data in low frequencies than in high frequencies.

This uneven distribution between the low and high frequencies is responsible for the blurring effect observed in the simple back-projection algorithm

First, it produces an image which has a high density in the center. This is due to the fact that many different images are being overlapped in this area.

Secondly, the resulting image is severely blurred. This effect comes from **the overlapping of the Fourier-transformed images** around the low frequency region. To control these effects, it is clear that a filter is needed during reconstruction of the projections.



## Filtered back projection

### Filter and Window

**Ramp filter**

Low Pass Filter (as known as smoothing), helps to remove high spatial frequency noise from the signal. High Pass Filter (as known as sharpening), helps to make an image appear sharper.

**The window function** is a mathematical function that is applied in an input signal to **minimize the spectral leakage**. In signal processing and statistics, when another function or waveform is multiplied by a window function, the product is zero-value outside in the line interval.

The typical operation of the filters is defined by: 

$q_{θ}(t)=\frac{1}{2π}\int P_{θ}(Ω)|Ω|W(Ω)e^{iΩ(xcosθ+ysinθ)}dθ$

$q_θ(t) = P_θ(t) ∗ h(t)$,  $h(t)$ is the convolution of filter and window function.



**Spatial aliasing**

**Zero padding uses avoid spatial aliasing c.f. circular convolution and linear convolution.**

In the Fourier domain, we derive the window function symmetric with respect to 0.

### Filtered Back Projection

We can find the one-dimension FFT of each angle of projection, and then it to integrates into the two-dimensional Fourier transform.

**Fourier slice theoren:** we have 1D fourier transform of the projection slice is equal to 2D fourier transform of the object functions at angle $θ$.

After accounting all the slides, we noticed the medial portion around the edges overlapped and converted into more noise and blurriness, which is inaccurate. 

We ought to adopt a filter function to declining to weigh the medial portion around the edges. 

We multiplied each of the FFT with the filter function, which is called **Filtered Back Projection**. As we cannot use only a filter for the Filtered back-projection inside the system as it integrates into infinity. Therefore we multiply the window function along with the filter to prevent the possible FFT numerical issue.

In principle we cannot recover the original image during reconstruction, there is always some data lost when resampled twice.

Then we can simply take the 2D inverse Fourier Transform and use different interpolation have our reconstructed image.

