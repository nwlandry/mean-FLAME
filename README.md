# mean-FLAME

Compiling for Mac: Using HomeBrew (after installing, of course), run the following in Terminal
```
brew install boost
brew install gsl
```

Now type
```
brew info boost
brew info gsl
```
to get <dur1> and <dir2>

```
g++ -I <dir1>/include -I <dir2>/include tevol_source_diff.cpp -o main -lboost_system
```

In my case,
dir1 is "/opt/homebrew/Cellar/boost/1.85.0"
dir2 is "/opt/homebrew/Cellar/gsl/2.8"

so

g++ -I /opt/homebrew/Cellar/boost/1.85.0/include -I /opt/homebrew/Cellar/gsl/2.8/include tevol_source_diff.cpp -o main -lboost_system -std=c++17


brew install cmake