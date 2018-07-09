### What is This?

An implementation of the [UPGMA](https://en.wikipedia.org/wiki/UPGMA) algorithm! It takes as input a matrix of traits, and computes the corresponding UPGMA tree and difference matrices for each step of computation. For displaying the tree, it uses the [ETE Toolkit](http://etetoolkit.org/), a Python framework for working with trees. It was initially built to make grading UPGMA problems on exams less tedious, and is thus not a particularly efficient implementation (cubic in the number of objects). 

### Usage Instructions

The main dependencies are Python 3, [NumPy](http://www.numpy.org/), [SciPy](https://www.scipy.org/), and [ETE3](http://etetoolkit.org/). It's probably easiest to use by setting up a local Conda environment and running it without having to worry about installing things globally. Start by [installing Miniconda](https://conda.io/docs/user-guide/install/index.html). Then, in whatever folder you have `upgma.py` in, set up and activate your environment by running
```shell
conda env create -f environment.yml
source activate env
```

This will create and activate a local environment named `env` with all the necessary dependencies.
Now, put the matrix of traits into a file named `traitmatrix.txt` in the same folder, and run

```
python upgma.py
```

Your tree should display, and your matrices should print to console and also be stored in `out.txt`

If you ever want to deactivate the local environment, run

```shell
source deactivate env
```


###TODOs

Here are things that still need doing! If you find yourself bored or interested, please do them! 

- Tree branch lengths need to be implemented. This isn't hard at all, just adding it in as per [Wikipedia](https://en.wikipedia.org/wiki/UPGMA). 
- I'm not sure how packages and dependencies and general usage of this program would work on Windows. If someone knows how, add it into the README and whatnot?
- This is ambitious, but what'd be super cool is if we can standardize a format for giving UPGMA problems and then use OCR to automatically read the student's trait matrix and grade the problem. Then, UPGMA problems could be graded by just scanning. Of course, it would still require a cross-check at the end, but it'd be far faster than our current system


