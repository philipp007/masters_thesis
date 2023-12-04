#!/bin/sh

echo "Cleaning..."

rm -fv core

rm -fv \
    precice-*.log \
    ./*.vtk

rm -fv \
    electrostatics/precice-electrostatics-events.json \
    electrostatics/precice-*.log \
    thermal/precice-thermal-events.json \
    thermal/precice-*.log \

echo "Cleaning parallel logs"

rm -fv \
    precice-electrostatics-events.json \
    precice-*.log \
    precice-thermal-events.json \
    precice-*.log \

rm -rfv precice-run


echo "Cleaning thermal + adapter configs"

rm -fv thermal/dolfin_files/*
rm -fv thermal/output/*
rm -fv thermal/precice_adapter_configs/*

echo "Cleaning electrostatics"

rm -fv electrostatics/output/*
rm -fv electrostatics/dolfin_files/*

echo "Cleaning complete!"
#------------------------------------------------------------------------------
