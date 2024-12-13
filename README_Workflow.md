# Workflow
Hi, I am Antoni, one of the previous contributors to this code-base. I refactored the code-base to remove bugs and improve performance. Another reason I did it is to simplify adding new features and code for newcomers. Now that all is plain and simple, I want to make sure it stays that way for as long as possible.

## Editing code
Before editing a single line or code or making a contribution, you should:
- Know C++ up to C++20, preferentially familiar with C++23. Not by heart, search [here](https://en.cppreference.com) for features you do not understand. Go [here](https://www.learncpp.com/) to learn C++ from scratch.
- Understand the code. Maybe not the entire code-base (although this is highly recommended), but at least the functions of the class trees you are working on. And the best way to understand it is to read it.

## Code formatting

### Formatter setup
The formatting that will be used across the entire code-base is stored [here](./NeuralNetworkCode/clang_format_config.txt). 
You can use any formatter you want that adheres to this clang configuration. The easiest way is to put this configuration in your VSCode settings `C_Cpp: Clang_format_style`.

### Variable, function and class naming
- Use camelCase for functions, same but first letter capital for classes, templates and structs.
- Use snake case for variables. Yes, this has currently not been enforced until now, but it is the guiding direction.

## Branches
The branches are set up to act as a path to eventually move your desired changes to `master` while minimizing the damage you may do changing the code.

### dev branch
- All the changes that will go into a single merge into `master` must gather here, and conflicts must be resolved here.
- The changes must be pushed from other branches, not this one.

### testing branch
- From `dev`, we merge into this branch to test the code. Once the code has passed these tests, you are free to merge into `master`.
- The tests must reproduce to a very high degree of precision. There currently are no standards for this, but if the behaviour is different, something is different. I recommend running the code with MVSC in Debug mode, examine the data structures. Look to what creates the divergence and whether the behaviour is not a bug.

### master branch
- The `master` branch is protected. No commits must be made here. If a fix is necessary, look into the `fix` branch.
- Once you merge here, you should merge `master` into every other branch.
- Make sure to make a `git tag` for the version.

### fix branch
- Nothing should break if you follow this process. If something does, make the fixes in this branch and merge into `master`.

### model-dev branch
- Same as `dev`, but for models. Forward all changes to `dev` before going into testing.

### refactor branch
- A branch where we refactor code. Must be pushed to `testing`.

### Paper branch
- Any branch that has a paper-like name, do not touch. It is the code-base's `master` branch where the paper's code was run.

## Merging good practices
`master <- testing <- dev` is the main pipeline. Tests must be passed in the `testing` branch.


## General good practices

### Ownership
- Avoid using shared pointers. 
- The class that has ownership over the object must have a `std::unique_ptr<T>` over the owned object.
- If another class needs to look into an already owned object, use a regular pointer (`T*`).
- If you have an array of objects that you want in the heap, do not make a `std::vector<std::unique_ptr<T>>`, a `vector<T>` is enough. Remember that a `std::vector<T>` only contains a pointer in its `.data` to the heap with the actual array.


### Variable naming
Use explanatory names in all your variables, no matter how minor they may be. This code-base currently follows the principle of "The code is the documentation". Abbreviations can be done, but it must be easy to come by an explanation for it for another user. Same goes for functions, classes, structs, templates...