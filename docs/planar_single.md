---
title: Planar instance
layout: default
---


<figure>
	<img src="drawings/pvtol_single.png" alt="PVTOL schematics" width="360">
	<figcaption>Fig1. - PVTOL schematics.</figcaption>
</figure>

## 1: Model 
<button type="button" class="btn btn-info" data-toggle="collapse" data-target="#newton_derivation">Newton derivation</button>
<div id="newton_derivation" class="collapse derivation" markdown="1">

Applying Newton's second law to our vehicle in ground frame, we get

$$\begin{align*}
m\ddot{x} &= (f_1+f_2)\sin(\theta) \\
m\ddot{z} &= -mg + (f_1+f_2)\cos(\theta) \\
J\ddot{\theta} &= d(-f_1+f_2)
\end{align*}$$


</div>

<button type="button" class="btn btn-info" data-toggle="collapse" data-target="#lagrange_derivation">Lagrange derivation</button>
<div id="lagrange_derivation" class="collapse derivation" markdown="1">

  * Generalized coordinates
 
$$q=\begin{pmatrix}x&z&\theta\end{pmatrix}^T $$

  * Kinetic energy

$$ T = \frac{1}{2} m(\dot{x}^2+\dot{z}^2) + \frac{1}{2} J \dot{\theta}^2$$

  * Potential energy

$$ V = mgz$$

  * Lagrangian

$$\begin{align*} 
L &= T-V \\
L &= \frac{1}{2} m(\dot{x}^2+\dot{z}^2) + \frac{1}{2} J \dot{\theta}^2 -  mgz
\end{align*}$$

### Partial derivatives

<table>
<tr><td>
$$
\begin{cases}
\frac{\partial{L}}{\partial{x}} =  0 \\
\frac{\partial{L}}{\partial{z}} =  - m g\\
\frac{\partial{L}}{\partial{\theta}} = 0 \\ 
\end{cases}
$$
</td><td>
$$
\begin{cases}
\frac{\partial{L}}{\partial{\dot{x}}} = m\dot{x} \\
\frac{\partial{L}}{\partial{\dot{z}}} = m\dot{z} \\
\frac{\partial{L}}{\partial{\dot{\theta}}} = J \dot{\theta}
\end{cases}
$$
</td></tr>
</table>

### Lagrange equations

 *
 
$$
\frac{d}{dt}\left( \frac{\partial{L}}{\partial{\dot{x}}} \right) - \frac{\partial{L}}{\partial{x}} = F_x
$$

$$
m\ddot{x} = -(f_1+f_2) \sin{\theta}
$$

 *

$$
\frac{d}{dt}\left( \frac{\partial{L}}{\partial{\dot{z}}} \right) - \frac{\partial{L}}{\partial{z}} = F_z
$$

$$
m\ddot{z} + mg = (f_1+f_2) \cos{\theta}
$$

 *

$$
\frac{d}{dt}\left( \frac{\partial{L}}{\partial{\dot{\theta}}} \right) - \frac{\partial{L}}{\partial{\theta}} = M_{\theta}
$$

$$
J\ddot{\theta} = d \left( -f_1+f_2 \right)
$$

 </div>
 
 
### State Space Representation

Using $$X = \begin{pmatrix}x&z&\theta&\dot{x}&\dot{z}&\dot{\theta}\end{pmatrix}^T$$ as state and $$U = \begin{pmatrix}f_1 & f_2 \end{pmatrix}^T$$ as input, a state space represeantation can be obtained as:


$$
\begin{equation}
\dot{X} = f(X,U) = \begin{pmatrix}
  \dot{x} \\
  \dot{z} \\
  \dot{\theta} \\
  -\frac{1}{m}  \sin{\theta} \left( f_1+f_2 \right) \\
  -g + \frac{1}{m}  \cos{\theta} \left( f_1+f_2 \right)\\
  \frac{d}{J} \left( -f_1+f_2 \right)
\end{pmatrix}
\end{equation}
$$

The following input variable change:
$$ 
U' = \begin{pmatrix}u_t\\u_d\end{pmatrix} = \begin{pmatrix}\frac{1}{m}(f_1+f_2) \\ \frac{d}{J}(-f_1+f_2)\end{pmatrix}
$$
leads to the simpler state space representation

$$
\begin{equation}
\dot{X} = f'(X,U') = \begin{pmatrix}
  \dot{x} \\
  \dot{z} \\
  \dot{\theta} \\
  -\sin{\theta}. u_t \\
  -g + \cos{\theta}. u_t\\
  u_d
\end{pmatrix}
\end{equation}
$$

[code](https://github.com/poine/these_ricardo/blob/main/src/single.py)


{%comment%}
### 2: Planning <button type="button" class="btn btn-info" data-toggle="collapse" data-target="#single_pvtol_planning">show</button>
<div id="single_pvtol_planning" class="collapse exemple" markdown="1">

</div>
{%endcomment%}



### 2: Control 

#### 2.1: Full State Feedback Regulation

<button type="button" class="btn btn-info" data-toggle="collapse" data-target="#single_pvtol_control">show</button>
<div id="single_pvtol_control" class="collapse exemple" markdown="1">

</div>

<figure>
	<img src="plots/single_step_x.apng" alt="step x" width="512">
	<figcaption>Fig1. - step x.</figcaption>
</figure>
[code](https://github.com/poine/these_ricardo/blob/main/src/single_test_2.py)

#### 2.2: Trajectory tracking
<button type="button" class="btn btn-info" data-toggle="collapse" data-target="#single_pvtol_control2">show</button>
<div id="single_pvtol_control2" class="collapse exemple" markdown="1">
</div>

<figure>
	<img src="plots/single_circle_tracking.apng" alt="circle tracking" width="512">
	<figcaption>Fig1. - circle tracking.</figcaption>
</figure>
[code](https://github.com/poine/these_ricardo/blob/main/src/single_test_3.py)
