# gni_predictors

This implementation is an alpha version library that will include different functional gene predictors. Only the tool NewGOA is provided. To use the gni_predictor tool we need to download from the GitHub repository. Then, into the repository, It is recommended to create a 'build' folder before compile and into it execute the following commands (make sure that you have "Armadillo" installed):
```
cmake ..
make
```
Now you should have a executable called `GNI_predictors`.

To check the different functions use `help`:

```
./GNI_predictors --help
```