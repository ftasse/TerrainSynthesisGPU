

INTRODUCTION

  This is a SpeedIT Classic toolkit, a GPU accelerated library with solvers
  for sparse linear systems of algebraic equations. 

  SpeedIT Classic library requires Nvidia CUDA 2.3 

  In the current version (1.0) library provides Conjugate Gradient solver and 
  sparse matrix dense vector multiply routine. 

LICENCE TERMS

Copyright (C) 2010  Vratis Ltd.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    	
    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
    or see <http://www.gnu.org/licenses/>. 	
    
    For any inquires regarding the toolkit contact us at
    Vratis Ltd.
    Muchoborska 18
    PL54424 Wroclaw, Poland
    or support <at> vratis.com

REQUIREMENTS

  SpeedITCLassic library requires Nvidia CUDA 2.3 installed. To obtain 
  CUDA please visit http://developer.nvidia.com/object/cuda_2_3_downloads.html


INSTALATION

1. Download the library

  Libray can be downloaded from svn

  svn checkout https://62.87.249.40/repos/speedIT/branches/1.0/SpeedIT_Classic/


  ON LINUX SYSTEMS
  
    2. Go into the following directory
    
      cd SpeedIT_Classic/
    
    3. Compile
    
      make
    
    4. If you wish you can copy libSpeedIT_Classic.so to system wide directory. 
    On Linux system this is usually /lib or /usr/lib.
    You may also copy header file si_classic.h to /usr/include
  
  
  ON WINDOWS SYSTEMS
  
    2. Open file SpeedIT_Classic/SpeedIT_Classic.sln with Visual Studio 
    
    3. Compile the project
