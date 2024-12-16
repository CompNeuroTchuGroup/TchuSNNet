# Features
## v1.0

## v.1.1

## v2.0

## v.2.1

## v3.0

## v.3.1

## v3.2
(In development)
- (Fixed) Appending suffixes in title when using iterate parameters will be available for the user.
- (Fixed) Random number selection improved in speed for RandomConnectivity.
- (Fixed) Now if you write `#` next to a parameter value without a space to separate them, that parameter value will still be read.
- (Fixed) A bug in the code made it so that unitary traversal of iterate parameters in Parameter.txt with all iterate parameter sets being composed of non-unitary parameters (parameters with more than 1 number to substitute), it would duplicate simulations, possibly crashing when going above the intended number of simulations due to overflowing the array.
  - Fix: Forcing title if there are iterate parameter sets, which is a unitary parameter.

- Other small improvements (string constants added, avoiding unnecessary copies).

- Still not done: 
  - Convert as many shared ptrs to unique as possible.
  - Fill this file for the other versions, and this one too.


