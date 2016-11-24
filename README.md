# Cytoscape Core: Building Guide

![](https://avatars1.githubusercontent.com/u/956141?v=3&s=200)

#### Status:
- 11/14/2016 - Updated for Cytoscape 3.5 release

## Introduction
Cytoscape is a fairly complex application and its core distribution has multiple repositories for managing its code.  This repository contains top-level pom file and utility script for building Cytoscape core distribution.  Most App developers won't need to clone this repository.  Keep reading below to learn about how to work with Cytoscape's source code.

### Target Audience
This document is a guide for developers who want to build the entire Cytoscape core distribution from scratch.  If you are interested in building Cytoscape apps, you don't need to build Cytoscape from source.  You can follow the guide here to learn more about Cytoscape app development:

* [Cytoscape Developer Documents](http://www.cytoscape.org/documentation_developers.html)

## Requirements
You need the following tools to build latest development version of Cytoscape 3:

* Computer with Windows, Mac, or Linux
* [JDK 8](http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html) or later
* [Maven 3](https://maven.apache.org/) (Tested with 3.3.x)
* [Git](https://git-scm.com/)
* [_cy.sh_](https://github.com/cytoscape/cytoscape-scripts) - Utility script for building Cytoscape core distribution.

While you can use any IDE to maintain Cytoscape 3, a popular IDE for this is Eclipse, which has its own Maven and Git support, too. However, for the initial repository clones and builds, we recommend that you follow the command line-based procedure below, and then switch to whichever IDE you prefer.

## Cytoscape 3 Core Directory Structure
Cytoscape source code is maintained in several GitHub repositories, and is organized into a main project and several sub-projects.

### The Core
```
── cytoscape
   ├── README.md
   ├── api
   ├── app-developer
   ├── cy.sh
   ├── gui-distribution
   ├── impl
   ├── parent
   ├── pom.xml
   └── support
```

* [parent](https://github.com/cytoscape/cytoscape-parent) - Contains several shared variables for building Cytoscape sub projects
* [api](https://github.com/cytoscape/cytoscape-api) - Public API bundles
* [impl](https://github.com/cytoscape/cytoscape-impl) - Implementation bundles
* [support](https://github.com/cytoscape/cytoscape-support) - Misc. bundles for building core
* [gui-distribution](https://github.com/cytoscape/cytoscape-gui-distribution) - Actual distribution created from core projects and core apps
* [app-developer](https://github.com/cytoscape/cytoscape-app-developers) - API JAR file for app developers

Instead of cloning each sub-project's repository one-by-one, you can use the utility script in the Cytoscape repository to initialize your workspace all at once (see below).

### New in 3.3.0: Core Apps
___Core Apps___ are Cytoscape apps originally from the core distribution.  They are located in their own separate GitHub repositories. Cytoscape depends on the latest version of each core app deployed to the Nexus repository, so you don't need to build core apps to build Cytoscape core.

**Note that each core app has its own repository and there is no parent-child relationship between Cytoscape Core and the Apps.  This means, in the core building process, local core apps will not be used in the process.**

As of Cytoscape 3.5.0 (November 2016), Cytoscape core distribution comes with the following core apps:

```
── apps
   ├── biopax
   ├── command-dialog
   ├── core-apps-meta
   ├── cyREST
   ├── datasource-biogrid
   ├── idmapper
   ├── json
   ├── network-analyzer
   ├── network-merge
   ├── opencl-cycl
   ├── opencl-layout
   ├── psi-mi
   ├── sbml
   ├── webservice-biomart-client
   ├── webservice-psicquic-client
   └── welcome
```

### Optional Projects
This repository contains sample code for app developers and it will not be included in the core distribution.  **You don't have to build this repository if you just want to build the Cytoscape Core distribution.**

* [samples](https://github.com/cytoscape/cytoscape-samples)

----

## Building Development Version of Cytoscape 3
Here is the step-by-step guide to build a development version of Cytoscape.

### Branch Management
#### Cytoscape Core
For the core projects, development version always uses the branch named **develop**.  **Master** branch is only for the final release.  If you want to build the latest development version of Cytoscape, you should use **develop** branch for all sub-projects.

#### Core Apps branch management
Since core apps have their own release cycles, they have different branching scheme.  Usually, features are developed in feature branches, and there is only one common branch called **master**.  Head of the master branch is always the latest development version of the core app.    

### Step 1: Clone the Main Project
1. Install required tools: JDK, Maven, and Git. On some systems, these may be preinstalled - you can use those versions if they are relatively recent, though we would recommend Oracle's JDK over OpenJDK.

1. Add JDK, Maven, and Git to your system PATH if necessary. On some platforms, this is done automatically on installation - try running mvn, git, or java at a command line to check this.

 If you need to add tools to the PATH, the steps you should follow vary by operating system. On Windows, this can be done in the Environment Variables dialog - open the file browser, right click on "Computer" or "This PC", select "Advanced system settings", then click "Environment Variables" - you should be able to edit the the PATH by selecting the PATH environment variable and clicking Edit. On Mac or Linux, you would need to edit the .profile or .bashrc file in your home directory to set environment variables, depending on the type of shell you are using. It may be helpful to add the following to .profile and set the PATH in .bashrc so that all shells will read the same values:

 ```
 if [ -f ~/.bashrc ]; then
   source ~/.bashrc
 fi
 ```

 Then, in .bashrc, add something like the following:

 ```
 export PATH=/path/to/java/bin:/path/to/maven/bin:/path/to/git/bin:$PATH
 ```

 Use the directories where binaries are located, as this wiill ensure that the command line knows where to find them.

1. Set the JAVA_HOME environment variable to the JDK 1.8 installation directory. This is only necessary if you have multiple versions of Java installed - if JDK 1.8 is the only one, Cytoscape will be able to automatically find it without the need for an environment variable. To do this, follow the same instructions as above, but for JAVA_HOME instead of PATH. On Windows, you may have to click the "Add..." button under System Variables if JAVA_HOME does not already exist. On Mac/Linux, you would add an additional line to .bashrc (or .profile if you set environment variables there) like the following.

 ```
 export JAVA_HOME=/path/to/java
 ```

 On Mac, you can use```$(/usr/libexec/java_home -v 1.8)``` instead of the actual path to automatically specify the latest 1.8 JVM installed.

1. If you are using Git for the first time, you need to set your name and e-mail address. To do this, use the following commands:
 ```
 git config user.name "Your Name"
 git config user.email yourname@yourname.com
 ```
 Substitute your actual name and e-mail address in the commands.

1. Generate an ssh key and set it up on GitHub. To do this, you will first need to be added to the Cytoscape GitHub project by one of the core developers. Then, click the arrow in the top-right hand corner of any GitHub page and choose "Settings". Click "SSH and GPG keys", then follow the instructions on the linked guide to generate an SSH key on your particular operating system. Once you have generated an SSH key, return to the original "SSH and GPG keys" page and add the generated key using "New SSH key".

1. Clone this repository: ```git clone https://github.com/cytoscape/cytoscape.git```
    - Now you can see a new directory named **cytoscape**:

```
cytoscape     <-- parent level directory
├── README.md
├── cy.sh
└── pom.xml
```

#### Shortcut
If you want to skip the following steps, you can use this command to clone and build all core projects and core apps:

```
./cy.sh init-all
```
Once finished, you can skip to Step 4.

### Step 2: Clone the Sub Projects
1. _cd_ to the cloned main project directory: ```cd ./cytoscape```
1. Execute of of the following command:
    - Create new folder under current working directory: ```./cy.sh init```
    - (Optional) Specify target directory: ```./cy.sh init /path/to/new/cytoscape/source/code```
1. Now you can see a new subdirectory (also) named **cytoscape**, which contains the sub projects:

```
cytoscape     <-- parent level directory
├── README.md
├── cy.sh
├── cytoscape     <-- subproject directory
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

### Step 3: Building Cytoscape
1. Go into the **cytoscape** subproject directory ```cd ./cytoscape```
1. Run Maven: ```mvn clean install```
    - Option: use ```mvn -fae clean install``` (... see below)
1. Have a coffee break...  It depends on your machine specification and internet connection speed, but will take 5-120 minutes.  When you build Cytoscape for the first time, it will take a long time because maven downloads all dependencies from the remote server.

### Step 4: Run the new build
Now you are ready to run the new
1. ```cd gui-distribution/assembly/target/cytoscape```
1. Run development version of Cytoscape:
   - Mac/Linux: ```./cytoscape.sh```
   - Windows: ```./cytoscape.bat```

The option `-fae` is short for "fail at end", which allows the build to continue even if unit tests fail. (Current unit tests
do fail during the first few compiles, but eventually pass.) When Maven
is done, you can find the application in `gui-distribution/assembly/target/cytoscape`.

Note that if you want to test the new build with a clean slate, we recommend to remove the entire ```~/CytoscapeConfiguration``` directory.

----

## Building Core Apps
All of the core apps are maintained in their own repository and if you want to try the latest version of the core app, you need to build them separately.

### Step 1: Checking out core apps
Assuming you are in the subproject directory of Cytoscape project (not the parent level), then ```./cy.sh apps``` will check out every core app into the ```apps``` subdirectory. Each is hosted in its own GitHub repository, and changes can be committed directly to each directory.  All of the core apps are hosted under this org account:

* [Cytoscape Consortium GitHub Repository](https://github.com/cytoscape)

### Step 2: Building Core Apps
All of core apps are independent to each other and there is no dependency among those.  To build the apps, you can just _cd_ to the app's directory and run:

```
mvn clean install
```

to build the latest version.  You can also use the following command from top-level directory to build all core apps:

```
./cy.sh build-apps
```

This command simply runs ```mvn clean install``` for each core app directory.

### Step 3: Install the new build
To test changes, simply install the JAR using the App Manager or copy to the ```~/CytoscapeConfiguration/3/apps/installed``` directory.

----

# Adding a new Core App
If you need to add an new core apps, you need to follow these steps:

## Step 1: Update directory structure
All of the Cytoscape core apps should respect the following conventions for consistency.

1. **Organize the directory**
  - The top level directory should contain only the minimal set of files.
  ```
.
├── .gitignore
├── README.md
├── pom.xml
└── src
          ├── main
          └── test
  ```
  - Create proper _.gitignore_ file
  - Add _README.md_ to briefly describe the app
  - Make sure _src_ contains only main and test sub directories
1. **Use standard group and package name**
  - _org.cytoscape_ should be used as the group ID
  - Package name also needs to be updated to follow the standard, like: ```org.cytoscape.YOUR_APP_NAME```

## Step 2: Update pom.xml
1. **Use standard name and ID**
    ```
    <groupId>org.cytoscape</groupId>
	  <artifactId>my-app</artifactId>
	  <version>3.5.0</version>
	  <packaging>bundle</packaging>
	  <name>My APP</name>
    ```
1. **Add all required plugin instructions**
  - The core app should work as an independent repository.  This means you should not depends on any parent pom file.
  - All instructions, including compiler arguments and parameters for BND plugin should be provided
  - Make sure you have proper information for the following sections:
    - repositories
    - distributionManagement
    - scm
    - build
1. **Update version number**
  - Use standard versioning convention.  For example, if your app is introduced for Cytoscape 3.5.0, set your app's version to 3.5.0.  **DO NOT USE _-SNAPSHOT_ SUFFIX**.

## Step 3: Add proper set of unit tests
If you don't have test suite, you should add it.  This is important to detect errors caused by future changes.

## Step 4: Test the new build locally
Once you finish all of the changes above, build it locally.  If everything looks OK, install it to Cytoscape and test the app.

## Step 5: Move your repository to the Cytoscape org account
Move your app's repository to Cytoscape org account: https://github.com/cytoscape

## Step 6: Setup Jenkins CI job
You need to setup a Jenkins job for CI.  You need an ID and password to setup a new job.

* Go to http://code.cytoscape.org/jenkins/
* Login to the server as an admin
* Copy the job from any of the existing coreapp jobs
* Name the job with standard prefix, e.g.,  _coreapp-Your-App-Name_
* Update the source code repository URL

## Step 7: Push your changes to the remote repository
If Jenkins job is ready, it should start the build once you push your changes to the remote GitHub repository.  If the build result is _success_, then Jenkins automatically deploy the new JAR file to the [Nexus repository](http://code.cytoscape.org/nexus/index.html#view-repositories;public~browsestorage).

## Step 8: Add your app to _assembly_ pom
Once the new app is in Nexus repository, you can use it from the core.  Add your app as a new dependency in **cytoscape-gui-distribution/assembly/pom.xml**:

* https://github.com/cytoscape/cytoscape-gui-distribution/blob/develop/assembly/pom.xml

It should be a new entry under _**artifactItems**_ section.

## Step 9: Build Cytoscape
Now you can build Cytoscape.  Run the new version and make sure your app is part of the distribution.

## Step 10: Bump up the app's version number
Change the version number of your app to 3.x.x-SNAPSHOT.


----
## Misc. Instruction for core developers

### Creating a new branch
At some time, it may be necessary to create a new branch for a Cytoscape release. This is the case if it is necessary to carry on multiple threads of development simultaneously (i.e. 3.5 is in a release candidate stage, but we want to continue development on features destined for 3.6.)  

In this case, you will need to create a new branch for each repository in Cytoscape based on the current working branch. Change to the root directory of the checked-out code, and run the following command to create a new branch and push this to GitHub (make sure all changes are committed beforehand).

```
git checkout -b 3.x.x
```

Substitute the desired name of the new branch for "3.x.x" (the version number of the release is recommended).

Then, repeat for each individual sub-repository (those being api, app-developer, gui-distribution, impl, parent, and support). Then, run the following in the root directory and each subdirectory to push the new branch to GitHub:

```
git push
```

Your new branch has been created and pushed.

### Choosing a Branch
If your Cytoscape project has Git branches, you can switch branches easily with **cy** script.  Simply go to the parent level  *cytoscape* directory and type:

```
cy switch BRANCH_NAME
```

where **BRANCH_NAME** is the name of the branch you want to switch.  All Cytoscape subprojects are following git-flow style branching scheme.  *Master* is used only for releases, and *develop* is the latest development branch.

### Managing Nested Git Repos with a GUI
The Cytoscape project is organized as a nested set of Git repositories. This provides for modularity and flexibility, but at the cost of greater complexity. The script at the parent level repository helps to manage ```git pull``` commands, but beyond that it can be challenging to manage the state of each repository. Fortunately, there are GUIs that can help:

* [GitHub Desktop](https://desktop.github.com/)
* Eclipse: [Setup instructions](../../wiki/Managing-Cytoscape-Git-Repos-in-Eclipse)
* [SourceTree](https://www.sourcetreeapp.com/)

----

# Building a New Release
This section is for core developers only.

## Introduction
There are several steps involved in building a release. Some of these can be automated using the [cytoscape-scripts](https://github.com/cytoscape/cytoscape-scripts) - however, these may not completely replicate the steps below and have not been updated in a while. As such, it is recommended to do the following steps manually.

## Updating version numbers
When preparing to make a new release (or start a new development branch), it is necessary to update the version numbers in Cytoscape to reflect this release. This is typically done at the release-candidate state, and also when updating the development branch for the next development version.

As the version number is repeated many times throughout the Cytoscape build configuration, it is easiest to make use of the Maven versions plugin to help make this change. To do this, change to the root directory in the Cytoscape directory structure. Then, run the following command:

```
mvn versions:set -DnewVersion=3.x.x
```

Substitute the desired version for "3.x.x" - if you are starting a new development branch, you would use a version number ending in -SNAPSHOT (i.e. 3.6.0-SNAPSHOT).

Then, change to the "parent" directory and repeat this command. When you are satisfied with the updated version number, run ```mvn versions:commit``` in both directories (run ```mvn versions:revert``` to undo the changes).

Though this will update most instances of the version number, it doesn't get them all. **You will need to manually update the version number in the following places:**

* cytoscape.sh (in gui-distribution/assembly/src/bin)
* cytoscape.bat (in gui-distribution/assembly/src/bin)
* parent/pom.xml (look for taglets)
* pom.xml in app-developer, gui-distribution, impl, support (look for properties tag)
* pom.xml in src/main/resources/archetype-resources/pom.xml (for each archetype subdirectory in support/archetypes)
* pom.xml in event-impl/it, model-impl/it, session-impl/impl, session-impl/integration-test, viewmodel-impl-parent/it, vizmap-impl-parent/it, work-swing-impl-parent/it

You can edit by hand, or use grep/sed to update these numbers.  Push and commit all changes to GitHub when you are done updating the version numbers.

## Releasing unreleased updates to core apps
Typically, updates to core apps will be released separately from the Cytoscape core development cycle, and the development branch will be updated to use any new updates as they are released. However, if a core app depends on an unreleased API in the development version of Cytoscape, this won't work. In that case, we have to release the unreleased core app when Cytoscape is being released.

After updating the version numbers to a release version (at the release candidate stage), check the gui-distribution/assembly/pom.xml to see if there are any -SNAPSHOT dependencies in core apps. If so, you will want to update these core apps to use non-SNAPSHOT API dependencies and a non-SNAPSHOT version number. Then, commit/push your changes and do an ```mvn deploy``` to deploy the new version. Finally, update the gui-distribution/assembly/pom.xml core app dependency to the new (non-SNAPSHOT) version number.  The updated app should be released to the App Store right before Cytoscape is released (take care to submit the same version as was deployed to Nexus - rebuilding will result in a non-matching checksum).

## Deploying to Nexus
When we are ready to release, we need to deploy artifacts for each bundle to Nexus. To build and deploy all artifacts, run the following command from the Cytoscape top-level directory:

```
mvn clean deploy
```
Note that **you will need to configure the Nexus server in ~/.m2/settings.xml before doing this.** Deploying to Nexus will always rebuild Cytoscape, so each deployment will have a different timestamp/SHA hash than the last deployment (even if nothing has changed).

## Building Installers
To build installers, we use a proprietary tool called install4j. This can be downloaded [here](https://www.ej-technologies.com/download/install4j/files) - the license key should be available to core team members. After installing install4j, you may need to update the path in the build configuration, which is set in the install4j.executable property in gui-distribution/packaging/pom.xml. The default value will work with the default install path on a Mac - this will need to be changed if you are building on a Windows/Linux system or install to a non-default path on Mac (it is advised not to commit changes to this property, as it is intended to be set locally).

To build installers for the most recently-built code, run the following command from the ```gui-distribution/packaging``` directory:

```
mvn install
```

### Signing Mac Installer
After that finishes, you will need to run the following command to sign the Mac DMG. This requires the Mac App Store certificate 'Developer ID Application' to be installed:

```
./sign-dmg.sh 'Developer ID Application'
```

After this is done, you would upload the built installers (in the target/install4j subdirectory of packaging, with the exception of the signed DMG which will be in target/install4j/signed) to the desired host (this would generally be chianti). If you are building a full release, also upload the swing-app-api JAR (in api/swing-app/api/target under the Cytoscape build root) and the API Javadocs (in app-developer/target/API) in the same directory.

## Merging a new release into the master branch
When a new release is cut, the release branch (or develop branch if that is being used) needs to be merged into the master branch. To do this, first make sure all your changes are checked in. Then, switch to the master branch. You can use the cy.sh script for this - to do that, run from the Cytoscape build root:

```
./cy.sh switch master
```

Then, in the root directory and each component subdirectory (api, app-developer, gui-distribution, impl, parent, and support), run the following command:

```
git merge branch_name
```

(substitute the branch you're merging into master for branch_name)
This should merge cleanly - if there are any conflicts you will need to resolve them. When you are satisfied, you will need to push all the changes - once again, you can use cy.sh from the build root:

```
./cy.sh push
```

### Tagging a Release
After merging a release into the master branch, it should be tagged with the release version number. To do this, we start in the Cytoscape build root with master checked out in all projects and execute the following commands:
```
git tag x.y.z
git push origin x.y.z
```
(substitute the version number for x.y.z)

This will tag the repository and push it to GitHub. Then, repeat this for all the sub-repositories (api, app-developer, gui-distribution, impl, parent, and support).

### Finalizing a Release
There are a few other steps that need to be completed when building a release. These should be done at the very end of the process (i.e. right before sending out the announcement and updating the website/release notes).

1. Update the version number in news.html, and add the announcement for a new release. This file is located at
http://chianti.ucsd.edu/cytoscape-news/news.html
1. Tag the manual to correspond with the new release. The manual is now a GitHub repository (located [here](https://github.com/cytoscape/cytoscape-manual)), and tagging it will create a new version of the document on ReadTheDocs. This is referenced by the Cytoscape application (using its internal version to determine the URL) - when tagging, the version number should not include any prefix or suffix and should always have three digits and two decimal places (so 3.6 should be "3.6.0").


### Updating core apps
When a core app update is ready to be pushed to the App Store and the development version of Cytoscape, there are a few steps you need to follow. Before releasing, the version number must be a non-SNAPSHOT version, and the core app must not depend on SNAPSHOT APIs.

1. Deploy to Nexus using ```mvn deploy``` - note that you will need to have our repository properly configured in ~/.m2/settings.xml to do this.
1. Update the gui-distribution/assembly/pom.xml file in Cytoscape core to depend on the new version of the core app.
1. Submit the new core app to the App Store using the web-based submission process at apps.cytoscape.org.
1. Update the core app's version for future development - i.e. if the last release was 3.3.1, update it to 3.3.2-SNAPSHOT.

Note: If a core app update depends on an unreleased version of the Cytoscape API, it cannot be released to the App Store. In that case, to add it to the development version you will want to follow the same steps as above, but skip the steps relating to App Store submission. You will leave the core app's version as a -SNAPSHOT release, and use that when updating the gui-distribution/assembly/pom.xml. When preparing to release Cytoscape, at that time you would remove the -SNAPSHOT and update the API dependencies to release versions (though the app shouldn't be submitted to the App Store until you are ready to release).

### The core apps meta-app
The core apps meta-app is a special core app - it contains no code itself besides a CyActivator, but depends on all the other core apps. The purpose of this is to allow us to add new core apps to each Cytoscape installation in between releases. This app should be updated whenever a core app is added to add the new app dependency to the Cytoscape-App-Dependencies section of the pom.xml. It can also be updated when a release is made (or even between releases) to update all the dependencies to the latest versions, but it is unclear that this is necessary as Cytoscape app dependencies refer to that version or greater.

----
## Tips and Notes for Core Developers

### Notes for Windows Developers
Windows implementations of Git and other tools differ slightly from the above.

* The Windows Git installer creates two application shortcuts: Git Bash and Git GUI. You should use Git Bash for command line operations.
* When executing the `cy init` script, be sure that your current path contains no blanks. The `cy` script's path parser does not understand blanks.
* You can follow the Git SSH instructions to create your SSH key, but when you start the SSH agent, use `eval $(ssh-agent)` instead of `eval 'ssh-agent' -s`.
* When running `cy init`, if you get "flags: FATAL unable to determine getopt version" somewhere in the output, you must be sure to put `getopt` in your PATH. The default location for `getopt` is `C:\Program Files (x86)\GnuWin32\bin`.

# Windows Notes

The Windows command line does not execute **cy** scripts, so you must execute it from within the Git Bash window within the Git GUI. However, the Git GUI doesn't present the Git Bash menu item until you have already created a repository. To get around this, create a new repository Test, and then start Git Bash.

The **cy init** script attempts to fetch repositories into the current directory. If the current directory has an element containing a blank (e.g., c:\users\My Account), **cy** will fail.

### Notes for All Developers
* Note that the `cy init` script accepts a path as a parameter. The path specifies where Cytoscape projects should be installed. Omitting the path reverts to the current working directory.
* Be sure you have installed Java JDK, not Java JRE.
* If you are developing on a virtual machine, be sure to configure around 8GB RAM and 50GB disk.
* To create a Cytoscape project in Eclipse (once you have run `cy init`), select File | Import, and then select Maven | Existing Maven Projects. Browse to the Cytoscape directory created by `cy init`, and note that all pom.xml files are found. To finish the import, wait for all projects to be created and compiled. This may take several minutes.
* To debug Cytoscape, follow this video: http://opentutorials.cgl.ucsf.edu/index.php/Tutorial:Remote_Execution_for_Debugging. To add all Cytoscape sources, use the Source tab in the Debug Configurations dialog, click the Add button, choose the Java Project container, and select all projects.
* To edit-compile-run, make your changes in the project you're working in. From Eclipse, you can Run As ... Maven Install. Eclipse will build the .class files automatically, so Maven's job is to create the .jar and promote it to private Maven repository. An unresolved compile issue will show in the Cytoscape console window when you run ... Maven doesn't complain, and Eclipse complains visually. Alternative: in Git Bash, set pwd to project directory (e.g., welcome-impl) and do `mvn clean install`.
* Valuable additional information: http://wiki.cytoscape.org/Cytoscape_3/CoreDevelopment
