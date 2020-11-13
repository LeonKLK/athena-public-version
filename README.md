# Athena++

Athena++ is a strophysical magnetohydrodynamics (MHD) code using C++. Our project focused on the the interaction between a planet and the accretion disk. Thus most of our changes comes from the `disk.cpp` in `pgen`(and the corresponding input file`athinput.disk_cyl`). Orginally, `disk.cpp` contains only an revolving accretion disk. We then modified different initial conditions, including velocity, pressure, density distribution, etc.

Besides initial conditions, in the `EnrollUserExplicitSourceFunction` we added our source function to investigate the interaction, including interaction between the disk and planet, Leapforg integrator for updating position of planet and cooling function at last.

If you have any questions regarding my code, feel free to ask me!

Remark: I DON'T OWN THIS CODE. Athena++ is developed by a remarkable team in Princeton University. If you have any questions regarding the performance of the code or any other problems unrelated to my above changes, please ask them in their Github directly. Thank you.



<!-- Jenkins Status Badge in Markdown (with view), unprotected, flat style -->
<!-- In general, need to be on Princeton VPN, logged into Princeton CAS, with ViewStatus access to Jenkins instance to click on unprotected Build Status Badge, but server is configured to whitelist GitHub -->
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

<!--[![Public GitHub  issues](https://img.shields.io/github/issues/PrincetonUniversity/athena-public-version.svg)](https://github.com/PrincetonUniversity/athena-public-version/issues)
[![Public GitHub pull requests](https://img.shields.io/github/issues-pr/PrincetonUniversity/athena-public-version.svg)](https://github.com/PrincetonUniversity/athena-public-version/pulls) -->

