## parallel_computation
OpenMP with C/C++ language 


1. Install clang-omp using homebrew:
    ```
    brew install llvm
    ```

2. Add llvm binaries to your path using:
    ```
    echo 'export PATH="/usr/local/opt/llvm/bin:$PATH"' >> ~/.bash_profile
    ```

3. Create a new Xcode project.

4. Under Build Settings
    * Add a new user-defined setting CC with the value
      ```
      /usr/local/opt/llvm/bin/clang
      ```
    * Add to Other C Flags
      ```
      -fopenmp
      ```
    * Add to Header Search Paths
      ```
      No
      ```

5. Under Build Phases: Add file to _Link Binary With Libraries_ from
```
/usr/local/Cellar/libomp/12.0.1/lib/libomp.dylib
```

6. In build settings > build options > Enable Index-While-Building Functionality to `No`
