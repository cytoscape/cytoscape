# Cytoscape Core Distribution: Building Guide

This repository contains top-level pom file and utility script for building entire Cytoscape core distribution.  Most App developers won't need to clone this repository.  Keep reading below to learn about how to work with Cytoscape's source code.

If you are interested in building Cytoscape apps, you don't need to build Cytoscape from source.  You can follow the guide here:

http://www.cytoscape.org/documentation_developers.html


## Getting Started with Cytoscape Source Code


### Requirements

You need the following tools to build latest development version of Cytoscape 3:

* [JDK 8](http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html) or later
* [Maven 3](https://maven.apache.org/) (Tested with 3.3.3)
* [Git](https://git-scm.com/)
* [_cy.sh_](https://github.com/cytoscape/cytoscape-scripts) - Utility script for building Cytoscape core distribution. 

While you can use any IDE to maintain Cytoscape 3, a popular IDE for this is Eclipse, which has its own Maven and Git support, too. However, for the initial repository clones and builds, we recommend that you follow the command line-based procedure below, and then switch to whichever IDE you prefer.

### Cytoscape 3 Core Main and Sub Projects
Cytoscape source code is maintained in several GitHub repositories, and is organized into a main project and several sub-projects:

* [parent](https://github.com/cytoscape/cytoscape-parent)
* [api](https://github.com/cytoscape/cytoscape-api)
* [impl](https://github.com/cytoscape/cytoscape-impl)
* [support](https://github.com/cytoscape/cytoscape-support)
* [gui-distribution](https://github.com/cytoscape/cytoscape-gui-distribution)
* [app-developer](https://github.com/cytoscape/cytoscape-app-developers)

Instead of cloning each sub-project's repository one-by-one, you can use the utility script in the Cytoscape repository to initialize your workspace all at once (see below).

### Optional Projects
* [samples](https://github.com/cytoscape/cytoscape-samples)

## Building Development Version of Cytoscape 3
Here is the step-by-step guide to build a development version of Cytoscape.

### Clone the Main Project
1. Install required tools: JDK, maven, git, and git-flow.
1. Clone this repository: ```git clone https://github.com/cytoscape/cytoscape.git```
    - Now you can see a new directory named **cytoscape**:

```
cytoscape     <= parent level directory
├── README.md
├── cy.sh
└── pom.xml
```

### Clone the SubProjects
1. CD to the cloned main project directory: ```cd ./cytoscape```
1. Execute the following command: 
    - Create new folder under current working directory: ```./cy.sh init```
    - (Optional) Specify target directory: ```./cy.sh init /path/to/new/cytoscape/source/code```
1. Now you can see a new sub-directory (also) named **cytoscape**, which contains the subprojects:

```
cytoscape     <= parent level directory
├── README.md
├── cy.sh
├── cytoscape     <= sub-project sub-directory
│   ├── README.md
│   ├── api
│   ├── app-developer
│   ├── cy.sh
│   ├── gui-distribution
│   ├── impl
│   ├── parent
│   ├── pom.xml
│   └── support
└── pom.xml
```

### Building Cytoscape
1. Go into the **cytoscape** sub-directory ```cd ./cytoscape```
1. Run Maven: ```mvn clean install```
    - Option: use ```mvn -fae clean install``` (... see below)
1. Have a coffee break...  It depends on your machine specification and internet connection speed, but will take 5-120 minutes. 
1. ```cd gui-distribution/assembly/target/cytoscape```
1. Run development version of Cytoscape: 
   - Mac/Linux: ```./cytoscape.sh```
   - Windows: ```./cytoscape.bat```

The option `-fae` is short for "fail at end", which allows the build to continue even if unit tests fail. (Current unit tests
do fail during the first few compiles, but eventually pass.) When Maven
is done, you can find the application in `gui-distribution/assembly/target/cytoscape`.

### Choosing a Branch
If your Cytoscape project has Git branches, you can switch branches easily with **cy** script.  Simply go to the main project *cytoscape* directory and type:

```
cy switch BRANCH_NAME
```

where **BRANCH_NAME** is the name of the branch you want to switch.  All Cytoscape sub-projects are following git-flow style branching scheme.  *Master* is used only for releases, and *develop* is the latest development branch.

### Managing Nested Git Repos with a GUI
The Cytoscape project is organized as a nested set of Git repositories. This provides for modularity and flexibility, but at the cost of greater complexity. The script at the parent level repository helps to manage ```git pull``` commands, but beyond that it can be challenging to manage the state of each repository. Fortunately, there are GUIs that can help:

* Eclipse: [Setup instructions](../../wiki/Managing-Cytoscape-Git-Repos-in-Eclipse)
* SourceTree: https://www.sourcetreeapp.com/


## New in 3.3.0: Core Apps
___Core Apps___ are Cytoscape apps originally from the core distribution.  They are located in their own separate GitHub repositories. Cytoscape depends on the latest version of each core app deployed to the Nexus repository, so you don't need to build core apps to build Cytoscape core.

### Checking out core apps
Assume you are in the top level directory of Cytoscape project. ```./cy.sh apps``` will check out every core app into the ```apps``` subdirectory. Each is hosted in its own GitHub repository, and changes can be committed directly to each directory. To test changes, simply install the JAR using the App Manager or copy to the ~/CytoscapeConfiguration/3/apps/installed directory

### Updating core apps
When a core app update is ready to be pushed to the App Store and the development version of Cytoscape, there are a few steps you need to follow. Before releasing, the version number must be a non-SNAPSHOT version, and the core app must not depend on SNAPSHOT APIs.

1. Deploy to Nexus using ```mvn deploy``` - note that you will need to have our repository properly configured in ~/.m2/settings.xml to do this.
1. Update the gui-distribution/assembly/pom.xml file in Cytoscape core to depend on the new version of the core app.
1. Submit the new core app to the App Store using the web-based submission process at apps.cytoscape.org.
1. Update the core app's version for future development - i.e. if the last release was 3.3.1, update it to 3.3.2-SNAPSHOT.

## Notes for Windows Developers
Windows implementations of Git and other tools differ slightly from the above.

* The Windows Git installer creates two application shortcuts: Git Bash and Git GUI. You should use Git Bash for command line operations.
 
* When executing the `cy init` script, be sure that your current path contains no blanks. The `cy` script's path parser does not understand blanks.

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




