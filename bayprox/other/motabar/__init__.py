"""
Package motabar.
   
Provides classes for regression/interpolation/smoothing using the 
MoTaBaR method (Moving Taylor Bayesian Regression).

Related publications:
  - J. Heitzig: Moving Taylor Bayesian Regression. Working Paper. 
    Potsdam Institute for Climate Impact Research (PIK), 2011.
    
B{Author}: Jobst Heitzig <heitzig@pik-potsdam.de>

B{Requires}: Python 2.6.5+, Numpy, Scipy, progressbar

B{Version}: 0.9

B{See}: U{The B{MoTaBaR} webpage<http://www.pik-potsdam.de/members/heitzig/motabar>}

B{To do}:
  - A lot - See current product backlog.
  - Choose appropriate License!
  
B{Bugs}:
  - Remain to be discovered.

B{License}: MIT Open Source License 

B{Copyright}: (c) 2010-2011 Jobst Heitzig.

B{Contributors}:
  - U{Jobst Heitzig<mailto:heitzig@pik-potsdam.de>}

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""


__name__ = "motabar"
"""The project name."""

__author__ = "Jobst Heitzig <heitzig@pik-potsdam.de>"
"""The primary author of pygeonetwork."""

__copyright__ = "Copyright (c) 2010-2011 Jobst Heitzig."
"""The copyright holder of pygeonetwork."""

__license__ = "MIT Open Source License"
"""The license governing the use and distribution of pygeonetwork."""

__url__ = "http://www.pik-potsdam.de/members/heitzig/motabar"
"""The URL for the package's homepage."""

__version__ = 0.9
"""Version number."""

__date__ = "2011-08-10"
"""The release date of this version."""

__docformat__ = "epytext en"
"""The epydoc documentation format for this file."""

#
#  Import classes
#
#from function import Function # B. Goswami: commented this out@23/05/2013
from .function import Function   # B. Goswami: add this for Python 3 @ 13/06/2019


