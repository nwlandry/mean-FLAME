# mean-FLAME

Compiling for Mac: Using HomeBrew (after installing, of course), run the following in Terminal
```
brew install cmake
brew install boost
brew install gsl
```

## life-cycle
From the top-level directory, run the following commands:
```
cd life-cycle
mkdir build
cd build
cmake ..
cmake --build .
```
so

g++ -I /opt/homebrew/Cellar/boost/1.85.0/include -I /opt/homebrew/Cellar/gsl/2.8/include tevol_source_diff.cpp -o main -lboost_system -std=c++17



