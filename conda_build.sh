#!/usr/bin/env bash

# Getting genipe's version to build
genipe_version=$1
if [ -z $genipe_version ]
then
    echo "usage: $0 VERSION" 1>&2
    exit 1
fi

# Creating a directory for the build module
mkdir -p conda_dist

# Creating a directory for the skeleton
mkdir -p skeleton
pushd skeleton

# Creating the skeleton
conda skeleton pypi genipe --version $genipe_version

# Checking that fetching genipe was successful
if [ $? -ne 0 ]
then
    echo "Error when creating skeleton for genipe version $genipe_version" 1>&2
    exit 1
fi

# The different python versions and platforms
python_versions="3.4 3.5 3.6"
platforms="linux-32 linux-64 osx-64"

# Building
for python_version in $python_versions
do
    # Building
    conda build --python $python_version genipe &> log.txt

    # Checking the build was completed
    if [ $? -ne 0 ]
    then
        cat log.txt
        echo "Error when building genipe $genipe_version (python" \
             "$python_version)" 1>&2
        exit 1
    fi

    # Fetching the file name of the build
    filename=$(grep -oP "anaconda upload \K(\S+)$" log.txt)

    # Checking the file exists
    if [ -z $filename ]||[ ! -e $filename ]
    then
        echo "Problem fetching file $filename" 1>&2
        exit 1
    fi

    # Converting to the different platforms
    for platform in $platforms
    do
        conda convert -p $platform $filename -o ../conda_dist &> build_log.txt

        # Checking the conversion was completed
        if [ $? -ne 0 ]
        then
            echo "Problem converting genipe $genipe_version (python" \
                 "$python_version) to $platform" 1>&2
            exit 1
        fi

        # Checking if a conversion was skipped due to same platform
        if egrep --quiet Skipping build_log.txt; then
            # Finding which platform was skipped
            missing=$(
                egrep Skipping build_log.txt |
                egrep -o "'([[:alnum:]]+-[[:digit:]]+)'" |
                uniq |
                sed -e "s/'//g"
            )
            mkdir -p ../conda_dist/$missing
            cp $filename ../conda_dist/$missing
        fi

    done

    # Purging
    conda build purge
done

popd
rm -rf skeleton

# Indexing
pushd conda_dist
conda index *
popd
