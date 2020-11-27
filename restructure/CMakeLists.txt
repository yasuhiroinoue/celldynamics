cmake_minimum_required(VERSION 2.8)
add_library(
  restructure STATIC
  cellrearrangement.cpp
  cellrearrangement2.cpp
  celldivision.cpp
  cellintersection.cpp
  vertexoverlap.cpp
  removal.cpp
  ../_parameters.h
)
