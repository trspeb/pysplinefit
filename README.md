PYSPLINEFIT
===========

Python library to fit a general spline to data in python/scipy.

Presentation
-------

Let's imagine you have experimental data in two arrays or matrix: x and y.
Let's imagine you would like to fit a lookup table on the data
on the knots (1, 2, 3, 4, 5).

Step 1: Create the basis. We choose a hat functions basis.

```
basis = Hatbasis((1, 2, 3, 4, 5))
```

Step 2: Fit the data.

```
spline = Optspline(basis, x, y)
```

Step 3: Use your spline as a function.

```
spline(2.5) # evaluate @ x = 2.5
```

Use cases
----

See the [handspinner article](res/article.pdf).
