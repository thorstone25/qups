# Contributing to QUPS
### Identifying a Bug
Not every unexpected result or error is a bug - before confirming that you've found a bug, check for the following gotchas:
- [ ] Are the units correct? For example, are *all* units in meters / seconds / Hertz ? Hint: use `UltrasoundSystem.plot()` as a quick sanity check!
- [ ] Do the data sizes correspond? For example, does `us.rx.numel == chd.N` and `us.seq.numPulse == chd.M` ? 
- [ ] Is this an unimplemented feature? For example, the f-k Stolt's migration beamformer expects a `TransducerArray` not a `TransducerConvex`.
- [ ] Is there a theoretical limitation? For example, a CFL > 0.3 may lead to instability or a sampling frequency < Nyquist may lead to aliasing.
- [ ] Are you running out of memory? This ... is probably not a bug, but can usually be handled by a for-loop, or using a `bsize` (block size) input if available.
- [ ] Is it a parallel pool incompatibility? This may be intentionally maintained for forward compatibility, as later MATLAB versions have [relaxed restrictions on parallel pools](https://www.mathworks.com/help/parallel-computing/release-notes.html).

### Submitting an Issue
If none of the above addresses the bug, please submit an issue, and include:
* What did you intend to do?
* What output did you expect?
* What output did you receive?
* (Optional) What steps have you taken to fix or circumvent the bug?

### Submitting a Patch
If you also have a patch that fixes a bug, after submitting an issue, feel free to open a pull request.

The patch does not need to follow all stylistic coding conventions, but any changes must pass all of the existing tests in QUPS. It is recommended but not required to add additional tests within an appropriate test class that addresses the bug.
Please reference the issue when submitting a pull request.

### New Features
If you are interested in a new feature, please feel free to open an issue describing what such a feature should accomplish.
Feature requests may be implemented in future versions as feasible. Simpler features may likely arrive sooner than more complex features. 
In any case, there is no gaurantee that a request will be implemented.

If you have implemented a new feature and would like for it to be incorporated, feel free to open a pull request.
For a new feature to be accepted, it should conform to the appropriate [class structure](src/README.md#class-structure) and [data format](src/README.md#data-format) conventions to maximize compatibility.

### Examples & How-tos 
QUPS is a toolbox, not a library: for most applications there is some assembly required. 
Example scripts that clarify syntax, applications, principles, or concepts are more than welcome! 
If you believe an example would be useful to the community, please feel free to create an issue to solicit an example, or a pull request to contribute one.

Example scripts should clearly demonstrate one or two ideas and complete in a few minutes at most.
If required, contributed examples should either create data from within QUPS (e.g. using `UltrasoundSystem.greens()`) or reference a downloadable public dataset.
Submissions which include data will be rejected.

