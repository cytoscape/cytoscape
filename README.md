# Cytoscape Core Distribution: Building Guide

This repository contains top-level pom file and utility script for building entire Cytoscape core distribution.  Most App developers won't need to clone this repository.  Keep reading below to learn about how to work with Cytoscape's source code.

If you are interested in building Cytoscape apps, you don't need to build Cytoscape from source.  You can follow the guide here:

http://wiki.cytoscape.org/Cytoscape_3/AppDeveloper


## Getting Started with Cytoscape Source Code


### Requirements

You need the following tools to build latest development version of Cytoscape 3:

* [JDK 8](http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html) or later
* [Maven 3](https://maven.apache.org/) (Tested with 3.3.3)
* [Git](https://git-scm.com/)
* _cy.sh_ - Utility script for building Cytoscape core distribution. 

### Cytoscape 3 Core Sub Projects
Cytoscape source code is maintained in several GitHub repositories:

* [parent](https://github.com/cytoscape/cytoscape-parent)
* [api](https://github.com/cytoscape/cytoscape-api)
* [impl](https://github.com/cytoscape/cytoscape-impl)
* [support](https://github.com/cytoscape/cytoscape-support)
* [app](https://github.com/cytoscape/cytoscape-app)
* [gui-distribution](https://github.com/cytoscape/cytoscape-gui-distribution)
* [app-developer](https://github.com/cytoscape/cytoscape-app-developers)

Instead of cloning each repository one-by-one, you can use utility script in this repository to initialize your workspace at once.

### Optional Projects
* [samples](https://github.com/cytoscape/cytoscape-samples)

## Building Development Version of Cytoscape 3
Here is the step-by-step guide to build Development version of Cytoscape.


### Clone all sub-projects
1. Install required tools: JDK, maven, git, and git-flow.
1. Clone this repository: ```git clone https://github.com/cytoscape/cytoscape.git```
1. CD to the cloned directory: ```cd ./cytoscape```
1. Execute the following command: 
    - Create new folder under current working directory: ```./cy.sh init```
    - (Optional) Specify target directory: ```./cy.sh init /path/to/new/cytoscape/source/code```
1. Now you can see a new directory named **cytoscape**:

```
    .
├── README.md
├── cy.sh
├── cytoscape
│   ├── README.md
│   ├── api
│   ├── app
│   ├── app-developer
│   ├── cy.sh
│   ├── gui-distribution
│   ├── impl
│   ├── parent
│   ├── pom.xml
│   └── support
└── pom.xml

8 directories, 6 files
```

### Building Cytoscape
1. Go into the _cytoscape_ directory ```cd ./cytoscape```
1. Run Maven: ```mvn clean install```
1. Have a coffee break...  It depends on your machine specification and internet connection speed, but will take 5-20 minutes. 
1. ```cd gui-distribution/assembly/target/cytoscape```
1. Run development version of Cytoscape: 
   - Mac/Linux: ```./cytoscape.sh```
   - Windows: ```./cytoscape.bat```

---

## New in 3.3.0: Core Apps
___Core Apps___ are Cytoscape apps originally from the core distribution.  They have their own GitHub repositories and ___app___ directory is a placeholder for them.  This means projects under _app_ are ___submodules___, or references to specific commits.  Because of this, you need to follow the special procidure to update contents in those directories.

### Updating Core Apps
Assume you are in the top level directory of Cytoscpae project.

1. ```cd app```
1. Go into the directory.  For example, ```cd cyREST```
1. Create new feature branch: ```git checkout -b my-new-feature```
1. Write your code
1. Add your changes to the branch: ```git add -A```
1. Commit new changes to the branch: ```git commit -m "YOUR COMMENTS HERE..."```
1. Go back to master: ```git checkout master```
1. Merge the branch: ```git merge my-new-feature```
1. Delete the feature branch: ```git branch -d my-new-feature```
1. Push your changes to upstream: ```git push```
1. ```cd ../```
1. Update the pointer to the new commit: ```git commit -am "YOUR COMMENTS HERE..."```
1. Push the changes to upstream: ```git push```

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

* The Windows Git installer creates two application shortcuts: Git Bash and Git GUI. You should use Git Bash for command line operations.
 
* When executing the `cy init` script, be user that your current path contains no blanks. The `cy` script's path parser does not understand blanks.

* You can follow the Git SSH instructions to create your SSH key, but when you start the SSH agent, use `eval $(ssh-agent)` instead of `eval 'ssh-agent' -s`.

* When running `cy init`, if you get "flags: FATAL unable to determine getopt version" somewhere in the output, you must be sure to put `getopt` in your PATH. The default location for `getopt` is `C:\Program Files (x86)\GnuWin32\bin`.

## Notes for All Developers

* Note that the `cy init` script accepts a path as a parameter. The path specifies where Cytoscape projects should be installed. Omitting the path reverts to the current working directory.

* Be sure you have installed Java JDK, not Java JRE.

* If you are developing on a virtual machine, be sure to configure around 8GB RAM and 50GB disk.

* To create a Cytoscape project in Eclipse (once you have run `cy init`), select File | Import, and then select Maven | Existing Maven Projects. Browse to the Cytoscape directory created by `cy init`, and note that all pom.xml files are found. To finish the import, wait for all projects to be created and compiled. This may take several minutes.

* To debug Cytoscape, follow this video: http://opentutorials.cgl.ucsf.edu/index.php/Tutorial:Remote_Execution_for_Debugging. To add all Cytoscape sources, use the Source tab in the Debug Configurations dialog, click the Add button, choose the Java Project container, and select all projects.

* To edit-compile-run, make your changes in the project you're working in. From Eclipse, you can Run As ... Maven Install. Eclipse will build the .class files automatically, so Maven's job is to create the .jar and promote it to private Maven repository. An unresolved compile issue will show in the Cytoscape console window when you run ... Maven doesn't complain, and Eclipse complains visually. Alternative: in Git Bash, set pwd to project directory (e.g., welcome-impl) and do `mvn clean install`.

* Valuable additional information: http://wiki.cytoscape.org/Cytoscape_3/CoreDevelopment


