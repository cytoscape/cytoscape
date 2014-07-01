Cytoscape
=========

This project contains metadata for working with all the Cytoscape core subprojects, each of which have their own Git repository.  Most people won't need to clone this repository.  Keep reading below to learn about how to work with Cytoscape's source code.

If you are interested in building Cytoscape apps, you don't need to build Cytoscape from source.  You can follow the guide here:

http://wiki.cytoscape.org/Cytoscape_3/AppDeveloper


# Getting Started with Cytoscape Source

## Requirements

You need the following tools to build Cytoscape 3:

* JDK 6 or later
* Maven 3.0.x series
    * The latest version of maven is 3.1.0, but we have not tested the build system with 3.1 as of 9/11/2013
* Git
* [git-flow](https://github.com/nvie/gitflow)
* [cy](https://github.com/cytoscape/cytoscape-scripts/releases/tag/1.2.0) - Utility script for Cytoscape developers

## Cytoscape 3 Core Subprojects
* [parent](https://github.com/cytoscape/cytoscape-parent)
* [api](https://github.com/cytoscape/cytoscape-api)
* [impl](https://github.com/cytoscape/cytoscape-impl)
* [support](https://github.com/cytoscape/cytoscape-support)
* [gui-distribution](https://github.com/cytoscape/cytoscape-gui-distribution)
* [headless-distribution](https://github.com/cytoscape/cytoscape-headless-distribution)
* [app-developer](https://github.com/cytoscape/cytoscape-app-developers)

## Optional Projects
* [samples](https://github.com/cytoscape/cytoscape-samples)

## Cloning the Cytoscape 3 Subprojects

1. Install JDK, maven, git, and git-flow.
1. Download latest version of [cy script](https://github.com/cytoscape/cytoscape-scripts/releases/) and unzip it to your local disk.
1. **cy** command is a shell script.  You need to change permission to execute it.
1. Execute the following command:

#### Core Developers

```
cy init
```

#### Other Developers

```
cy -r init
```

This clones read-only repository from github.


Now you can find a new directory named **cytoscape**.  It should contains the following:


- README.md
- api
- app-developer
- gui-distribution
- headless-distribution
- impl
- parent
- pom.xml
- support


### Choosing a Branch
Switching branches is easy with **cy** script.  Simply go to the top level directory and type:

```
cy switch BRANCH_NAME
```

where **BRANCH_NAME** is the name of the branch you want to switch.  All Cytoscape sub-projects are following git-flow style branching scheme.  *Master* is used only for releases, and *develop* is the latest development branch.

## Building the Core
From the top directory, type:
```
mvn -fae clean install
```

The option `-fae` is short for "fail at end", which allows the build to continue even if unit tests fail.  When Maven
is done, you can find the application in `gui-distribution/assembly/target/cytoscape`.


## Notes for Windows Developers
Windows implementations of Git and other tools differ slightly from the above.

- The Windows Git installer creates two application shortcuts: Git Bash and Git GUI. You should use Git Bash for command line operations.
- When executing the ''' cy init ''' script, be user that your current path contains no blanks. The '''cy''' script's path parser does not understand blanks.
- 
