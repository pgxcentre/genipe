#!/usr/bin/env bash

# Creating a directory for the skeleton
mkdir -p skeleton
pushd skeleton

# Creating the skeleton
conda skeleton pypi genipe

# The different python versions and platforms
python_versions="3.4 3.5"
platforms="linux-32 linux-64 osx-64"

# Building
for python_version in $python_versions
do
    # Building
    conda build --python $python_version genipe &> log.txt
    filename=$(egrep "^# [$] anaconda upload \S+$" log.txt | cut -d " " -f 5)

    # Converting
    for platform in $platforms
    do
        conda convert -p $platform $filename -o ../conda_dist
    done
done

popd
rm -rf skeleton

# Indexing
pushd conda_dist
conda index *
popd
