jdock
=====
![build workflow](https://github.com/stcmz/jdock/actions/workflows/build.yml/badge.svg)
![release workflow](https://github.com/stcmz/jdock/actions/workflows/release.yml/badge.svg)

jdock is an extended variant of idock which was originally developed by [@HongjianLi](https://github.com/HongjianLi) and is distributed under the same license.

>idock is a standalone tool for structure-based [virtual screening] powered by fast and flexible ligand docking. It was inspired by [AutoDock Vina], and is hosted on GitHub at https://GitHub.com/HongjianLi/idock under [Apache License 2.0]. idock is also available as a web server at [istar].

jdock keeps full compatibility with idock and the idock branding remains unchanged across the source code. The only naming change is the name of output binary and GitHub repository.

Features
--------

jdock inherits all features of idock in branch v2.x and tracks changes in that branch. Along with some bug fixes, jdock also adds its own functionalities including but not limited to:
* non-standard residue elimination,
* protonation with pKa values,
* per residue energy summarization and emission,
* alternate location indicator choosing,
* non-compulsory RF-score calculation,
* precision mode to avoid the use of grid maps,
* scoring and docking in a single run,
* compatibility with all kinds of line feedings.


Supported operating systems and compilers
-----------------------------------------

* Linux x86_64 and g++ 8.3.1
* Mac OS X x86_64 and clang 8.0.0
* Windows x86_64 and Visual Studio 2019


Compilation from source code
----------------------------

### Get Boost

jdock depends on the `Program Options` and `Asio` components in [Boost C++ Libraries]. Boost 1.75.0 was tested. There are several ways to get Boost.

#### With `vcpkg` on Windows, macOS or Linux:
```
# Note: this will download and build from source
vcpkg install boost-program-options boost-asio
```

#### With `nuget` on Windows:
```
nuget install boost_program_options-vc142
```

#### With `apt` on Ubuntu/Debian:
```
sudo apt install libboost-program-options-dev
```

#### With `Homebrew` on macOS:
```
brew install boost
```

#### With `PowerShell` on Windows:
```
$Url = "https://sourceforge.net/projects/boost/files/boost-binaries/1.75.0/boost_1_75_0-msvc-14.2-64.exe"
(New-Object System.Net.WebClient).DownloadFile($Url, "$env:TEMP\boost.exe")
Start-Process -Wait -FilePath "$env:TEMP\boost.exe" "/SILENT","/SP-","/SUPPRESSMSGBOXES"
```

#### Build from Source on Linux or Windows:

Download [Boost 1.75] and unpack the archive to `boost_1_75_0/include`.

Build on Linux run:
```
cd boost_1_75_0/include
./bootstrap.sh
./b2 --build-dir=../build/linux_x64 --stagedir=../ -j 8 link=static address-model=64
```

Or, on Windows run:
```
cd boost_1_75_0\include
bootstrap.bat
b2 --build-dir=../build/win_x64 --stagedir=../ -j 8 link=static address-model=64
```

Then add the path of the `boost_1_75_0` directory the to the BOOST_ROOT environment variable.

### Build with CMake

This project uses cross-platform build system CMake to build from source. It detects your environment and decides the most appropriate compiler toolset. The minimum version of CMake required is `3.20`. To build, simply run
```
cmake -B build
cmake --build build
```

The generated objects and executable will be placed in the `build` folder.

Optionally, on Linux or macOS one may install the output binary to the system (usually `/usr/local/bin`) by running
```
sudo cmake --install build
```

On Windows, the script should be run without sudo but under Administrator. The executable will be copied to an individual directory under `Program Files`.


### Build with Visual Studio

Visual Studio 2019 solution and project files are provided. To compile, simply run
```
msbuild /t:Build /p:Configuration=Release
```

Or one may open `idock.sln` in Visual Studio 2019 and do a full rebuild.

The generated objects will be placed in the `obj` folder, and the generated executable will be placed in the `..\bin` folder.


Usage
-----

First add jdock to the PATH environment variable.

To display a full list of available options, simply run the program without arguments
```
jdock
```

The `examples` folder contains several use cases. For example, to dock the ligand TMC278 against HIV-1 RT of PDB ID 2ZD1,
```
jdock -r receptors/2ZD1.pdbqt -l ligands/T27 -x 49.712 -y -28.923 -z 36.824 --size_x 18 --size_y 18 --size_z 20
```

Or one can instruct jdock to load the options from a configuration file
```
cd examples/2ZD1/T27
jdock --config idock.conf
```


Change Log
----------

### 2.2.3c (2021-05-17)

* Switched build system from Make to cross-platform CMake
* Added GitHub Actions workflows for building and releasing binaries

### 2.2.3b (2020-04-20)

* Merged changes from upstream HongjianLi/idock/v2.x for adding support for ten rare chemical elements
* Added --ignore_errors program option

### 2.2.3a (2019-06-14)

* Initial public release based on idock v2.2.3
* Added computation of residue energy contribution
* Added support for protonation via pka
* Added per residue output for gauss1/gauss2/gauss/repulsion/steric/hydrophobic/hbonding/nonsteric/total scoring combinations
* Added program option shortcuts
* No RF score computation by default (use --rf_score to revert)
* Removed precompiled binaries in conformance with the US export control regulation
* More changes see the Feature section above

Reference
---------

### jdock
**Maozi Chen**, Zhiwei Feng, Siyi Wang, Weiwei Lin, Xiang-Qun Xie. **MCCS, a scoring function-based characterization method for protein-ligand binding**. *Briefings in Bioinformatics*. 2020, October 14. [DOI: 10.1093/bib/bbaa239](https://doi.org/10.1093/bib/bbaa239) [PubMed: 33051641](https://pubmed.ncbi.nlm.nih.gov/33051641/)

### idock
Hongjian Li, Kwong-Sak Leung, and Man-Hon Wong. **idock: A Multithreaded Virtual Screening Tool for Flexible Ligand Docking**. *2012 IEEE Symposium on Computational Intelligence in Bioinformatics and Computational Biology (CIBCB)*, pp.77-84, San Diego, United States, 9-12 May 2012. [DOI: 10.1109/CIBCB.2012.6217214]


Author
--------------

### jdock
[Maozi Chen]

### idock
[Jacky Lee]


[virtual screening]: http://en.wikipedia.org/wiki/Virtual_screening
[docking]: http://en.wikipedia.org/wiki/Docking_(molecular)
[AutoDock Vina]: http://vina.scripps.edu
[Apache License 2.0]: http://www.apache.org/licenses/LICENSE-2.0
[istar]: http://istar.cse.cuhk.edu.hk/idock
[Boost C++ Libraries]: http://www.boost.org
[doxygen]: http://www.doxygen.org
[semantic versioning]: http://semver.org
[CodePlex]: http://idock.codeplex.com
[DOI: 10.1109/CIBCB.2012.6217214]: http://dx.doi.org/10.1109/CIBCB.2012.6217214
[Jacky Lee]: http://www.cse.cuhk.edu.hk/~hjli
[Maozi Chen]: https://www.linkedin.com/in/maozichen/
[Boost 1.75]: https://www.boost.org/users/history/version_1_75_0.html
