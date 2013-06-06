Cytoscape
=========

This project contains metadata for working with all the Cytoscape core subprojects, each of which have their own Git
repository.  Most people won't need to clone this repository.  Keep reading below to learn about how to work with
Cytoscape's source code.

If you are interested in building Cytoscape apps, you don't need to build Cytoscape from source.  You can follow the
guide here:

http://wiki.cytoscape.org/Cytoscape_3/AppDeveloper

# Getting Started with Cytoscape Source

## Requirements

You need the following software to build Cytoscape 3:

* JDK 6 or later
* Maven 3
* Git

If you want to clone all sub-projects at once, you also need:

* [Repo](http://code.google.com/p/git-repo/)

## Cytoscape 3 Core Subprojects
* parent
* api
* impl
* support
* gui-distribution
* headless-distribution
* app-developer
* samples

## Cloning the Cytoscape 3 Subprojects

```
repo init -u git@github.com:cytoscape/cytoscape.git
repo sync
```

### Choosing a Branch

This command will clone all of the subprojects but won't select a particular branch.  If you want to work with the
source of the latest stable release, checkout the `master` branch.  For example, if you're using the bash shell or
something similar, you can use the following command:

```
for PROJECT in api app-developer gui-distribution headless-distribution impl parent samples support
do
    pushd ${PROJECT}
    git checkout -b master origin/master
    popd
done
```

To work with the source of the latest integration build, you'll need to checkout the `develop` branch for
each subproject:

```
for PROJECT in api app-developer gui-distribution headless-distribution impl parent samples support
do
    pushd ${PROJECT}
    git checkout -b develop origin/develop
    popd
done
```

To see which branch you're on, you can run `repo status` on the top-level clone:

```
$ repo status
project api/                                    branch master
project app-developer/                          branch master
project gui-distribution/                       branch master
project headless-distribution/                  branch master
project impl/                                   branch master
project parent/                                 branch master
project samples/                                branch master
project support/                                branch master
```

If `repo status` gives you this output:
```
nothing to commit (working directory clean)
```
...that means you haven't checkout out a branch.

While it's possible to clone the subprojects manually without using `repo`, we don't recommend this.

## Building the core
From the top directory, type:
```
mvn -fae clean install
```

The option `-fae` is short for "fail at end", which allows the build to continue even if unit tests fail.  When Maven
is done, you can find the application in `gui-distribution/assembly/target/cytoscape`.
