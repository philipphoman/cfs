# Reversal learning of threat under continuous flash suppression

# Getting Started
This repository contains all the data and analysis
code to reproduce the manuscript "Affective flexibility without
perceptual awareness". These instructions describe how to obtain a copy
of the project up and running on your local machine for reproducing the
analysis described in the manuscript. The repository contains a Makefile
which reflects the dependencies of the analysis; analysis, figures and
manuscript can be produced by simply typing 'make' from the Unix command
line.

## Installing
Once the prerequisites are met, running the analyses should be
straightforward. The analysis file is in src/. 

## Prerequisites
The project was developed and tested on a Linux Ubuntu 17.10 machine. To
produce the manuscript from the command line, GNU emacs and texlive must
be installed together with some additional emacs and texlive
packages. Specifically, the following software and associated libraries
are required:

### texlive
-   xetex
-   latexmk

### emacs
-   org-mode
-   ess
-   org-ref

### R
-   pacman
-   R.matlab
-   cowplot
-   ellipse
-   astsa
-   magick
-   here
-   png
-   grid
-   flexmix
-   DescTools
-   boot
-   lme4
-   lsmeans
-   graphics
-   tibble
-   tidyr
-   ggplot2
-   dplyr
-   MASS

# Running the analysis
Change to the cfs directory and run 'make analysis'.

# Producing the figures
Change to the cfs directory and run 'make figures'. The figures can then
be found in output/figures.

# Producing the manuscript
Change to the cfs directory and run 'make manuscript'. The manuscript
will be in src/cfs<sub>ms.pdf</sub>

# License
This project is licensed under the MIT License - see the
[LICENSE.md](LICENSE.md) file for details

# Acknowledgments
We thank Patrik Vuilleumier who created and shared the spider
stimuli. This work was supported in part through the computational
resources and staff expertise provided by Scientific Computing at the
Icahn School of Medicine at Mount Sinai. Funding was provided by NIMH
grant MH105515 and a Klingenstein-Simons Fellowship Award in the
Neurosciences to Daniela Schiller.; ERC Advanced Grant XSPECT-DLV-692739
to David Carmel (Co-I); and Swiss National Science Foundation grant SNF
161077 to Philipp Homan. The funding source had no role in the design
and conduct of the study; collection, management, analysis, and
interpretation of the data; preparation, review, or approval of the
manuscript; and decision to submit the manuscript for publication.

# Built With
Org-mode 9.1.9.

