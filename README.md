# Three-Body-Problem

<h4>Authors of this code:</h4>

- Rodrigo Eduardo Esparza Rivera
- Javier Alejandro Pérez Garza
- Hugodiego Jimenez Estudillo
- Gustavo Lozano Villarreal
- Oliver Drogi Tello

## Code Explanation
In this Matlab code we solve the three body problem with a Sun-Jupiter-asteroid system. We start from the Universal Law of Gravitation in its vectorial form,

$$
\vec{F}(\vec{r})  = -\frac{GMm}{|\vec{r}|^ 2}\hat{r}.
$$

Simplifying the expression we can get:

$$
m_i\ddot{\vec{r}}  = -\frac{Gm_im_j}{|\vec{r}|^ 3}\vec{r}.
$$

Now, modelling the forces of the three gravitationally interacting bodies we can set the following system of equations:

$$\vec{F_{1}}=-\vec{F_{12}}-\vec{F_{13}},$$
$$\vec{F_{2}}=-\vec{F_{21}}-\vec{F_{23}},$$
$$\vec{F_{3}}=-\vec{F_{31}}-\vec{F_{32}}.$$

By definition, $\vec{r_{ij}}=\vec{r_i}-\vec{r_j}$. Finally, the system of equations which we solved using Verlet algorithm is in the form of:

$$
\ddot{\vec{r_1}}  = -\frac{Gm_2}{|\vec{r_1}-\vec{r_2}|^ 3}\left(\vec{r_1}-\vec{r_2}\right)-\frac{Gm_3}{|\vec{r_1}-\vec{r_3}|^ 3}\left(\vec{r_1}-\vec{r_3}\right),
$$

$$
\ddot{\vec{r_2}}  = -\frac{Gm_1}{|\vec{r_2}-\vec{r_1}|^ 3}\left(\vec{r_2}-\vec{r_1}\right)-\frac{Gm_3}{|\vec{r_2}-\vec{r_3}|^ 3}\left(\vec{r_2}-\vec{r_3}\right),
$$

$$
\ddot{\vec{r_3}}  = -\frac{Gm_1}{|\vec{r_3}-\vec{r_1}|^ 3}\left(\vec{r_3}-\vec{r_1}\right)-\frac{Gm_2}{|\vec{r_3}-\vec{r_2}|^ 3}\left(\vec{r_3}-\vec{r_2}\right).
$$

