This is the primary Maven project for building the Cytoscape Desktop application.

To report bugs in this or other Cytoscape Desktop sub-projects, please use the bug report form [here](https://cytoscape.org/bug-report.html).

# Cytoscape Core: Building Guide

![](https://avatars1.githubusercontent.com/u/956141?v=3&s=200)

#### Status:
- 4/2023 - Updated for 3.10.0 release
- 3/2025 - Updated for 3.10.3 release

## Introduction
Cytoscape is a fairly complex application and its core distribution has multiple repositories for managing its code.  This repository contains top-level pom file and utility script for building Cytoscape core distribution.  Most App developers won't need to clone this repository.  Keep reading below to learn about how to work with Cytoscape's source code.

### Target Audience
This document is a guide for developers who want to build the entire Cytoscape core distribution from scratch.  If you are interested in building Cytoscape apps, you don't need to build Cytoscape from source.  You can follow the guide here to learn more about Cytoscape app development:

* [Cytoscape Developer Documents](http://www.cytoscape.org/documentation_developers.html)

## Requirements
You need the following tools to build latest development version of Cytoscape 3.10:

* Computer with Windows, Mac, or Linux
* [JDK 17](https://adoptium.net/temurin/releases/?version=17) — but see **Which JDK** below if you also build the core apps
* [Maven 3](https://maven.apache.org/)
* [Git](https://git-scm.com/)
* _cy.sh_ - Utility script for building Cytoscape core distribution (available in this repository).

### Which JDK

**The core and the core apps do not build with the same JDK.** There is no single version that does both, so pick according to what you are building and switch ```JAVA_HOME``` when you move between them.

| Command | JDK required | Basis |
| --- | --- | --- |
| ```pull```, ```pull-apps```, ```switch```, ```switch-apps```, ```branch```, ```branch-apps```, ```status```, ```push```, ```reset``` | any (none) | These only run _git_. No Java is involved. |
| ```build``` (the core) | **17 or newer** | _parent/pom.xml_ compiles with ```<release>17</release>```, which Maven cannot honour on an older JDK. Verified on JDK 17. |
| ```build-apps``` (the core apps) | **11**, and **must be older than 16** | The apps compile for Java 11, and they pin _maven-bundle-plugin_ 4.1.0, which fails with a ```ConcurrentModificationException``` on JDK 16 and newer. Confirmed failing on both JDK 17 and JDK 21. |

So a normal core-development setup needs only **JDK 17**. If you also intend to run ```./cy.sh build-apps```, install **JDK 11** alongside it and point ```JAVA_HOME``` at 11 for that command only:

```
JAVA_HOME=/path/to/jdk-11 ./cy.sh build-apps
```

On Mac, ```/usr/libexec/java_home -V``` lists the JDKs you have installed and ```JAVA_HOME=$(/usr/libexec/java_home -v 11)``` selects one.

This is a limitation of the core apps' own build configuration, not of _cy.sh_ — a plain ```mvn clean install``` inside an app directory fails the same way on a modern JDK. Fixing it properly means bumping _maven-bundle-plugin_ in each app's own repository.

While you can use any IDE to maintain Cytoscape 3, a popular IDE for this is Eclipse, which has its own Maven and Git support, too. However, for the initial repository clones and builds, we recommend that you follow the command line-based procedure below, and then switch to whichever IDE you prefer.

### Related repositories
All build scripts are moved to the following repository:

* [Cytoscape Admin Scripts](https://github.com/cytoscape/cytoscape-admin-scripts)

## Services

Note that some Cytoscape functions rely on code deployed as services available on web servers. Generally, such services are callable by Cytoscape or directly by non-Cytoscape clients (e.g., Python) in the larger bioinformatics community. Some services are provided by other organizations (e.g., PSICQUIC for importing public networks), while others are provided by Cytoscape developers (e.g., Diffusion) and are located in or rely on other GitHub repositories. Here is a list of known external repositories containing services called by Cytoscape and maintained by Cytoscape core developers:

* [CXMate](https://github.com/cxmate/cxmate) - adapters that simplify service writing
* [Diffusion](https://github.com/idekerlab/heat-diffusion) - called by Diffusion core app

Each repo contains information on how to build and deploy the service.

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
* [app-developer](https://github.com/cytoscape/cytoscape-app-developer) - API JAR file for app developers

Instead of cloning each sub-project's repository one-by-one, you can use the utility script in the Cytoscape repository to initialize your workspace all at once (see below).


## Building Development Version of Cytoscape 3
Here is the step-by-step guide to build a development version of Cytoscape.

#### tl;dr
```
git clone https://github.com/cytoscape/cytoscape.git
cd cytoscape
./cy.sh pull
mvn -fae install -U -Dmaven.test.skip=true
./gui-distribution/assembly/target/cytoscape/cytoscape.sh
````
[Eclipse Users](https://github.com/cytoscape/cytoscape/wiki/Importing-Git-Repos-in-Eclipse) - Eclipse Import Instructions

*NOTE: Build order matters, and ```./cy.sh build``` takes care of it for you — it builds api/event-api, then api, then support, then impl, and only then the whole project. This is not optional: api/pom.xml and impl/pom.xml bind maven-source-plugin's ```aggregate``` goal, which forks a separate build that resolves dependencies from your local maven repository instead of from the reactor. On a first build nothing is in that repository yet, so a plain ```mvn install``` from the top dies at the second of 122 modules looking for event-api. Building the inner pieces first puts them in the repository so those forked builds can find them.*

*ANOTHER NOTE: Test compilation cannot be skipped — 34 of the poms depend on another module's test-jar. Use ```-DskipTests``` (compiles tests, does not run them), never ```-Dmaven.test.skip=true``` (does not compile them at all, so the test-jars are never produced and everything needing one fails). ```./cy.sh build``` already does this.*

*ANOTHER NOTE: If you see errors about "Could not transfer artifact" and "Blocked mirror for repositories," then you may be running a newer version of Maven that doesn't work for our repo at this time. Try Maven v3.6.0*

### cy.sh commands
Cytoscape core is spread over seven repositories, and the core apps over twenty more. The _cy_ script drives all of them together, so you rarely run _git_ or _maven_ in a single repository by hand. Every command below is run from inside your clone of this repository and acts on **all** of the repositories in its column:

| | Core repositories | Core apps |
| --- | --- | --- |
| Get them locally, and update them later | ```./cy.sh pull``` | ```./cy.sh pull-apps``` |
| Switch to an existing branch | ```./cy.sh switch BRANCH``` | ```./cy.sh switch-apps BRANCH``` |
| Create a new branch | ```./cy.sh branch NEW_BRANCH``` | ```./cy.sh branch-apps NEW_BRANCH``` |
| Build | ```./cy.sh build``` | ```./cy.sh build-apps``` |

The core repositories sit in the root of your clone; the core apps sit under ```apps/```. The two columns are independent of each other — the core repositories track **develop** while the core apps track **master** — so a command in one column never touches the other.

If you are unsure whether you want ```switch``` or ```branch```, see [Optional: Working on a different branch](#optional-working-on-a-different-branch) — the short version is that ```switch``` only moves to a branch that already exists and ```branch``` only creates one that does not.

Also available: ```status```, ```push```, ```reset```, ```run-all``` and ```validate-apps```.

### Branch Management
#### Cytoscape Core
For the core projects, development version always uses the branch named **develop**.  **Master** branch is only for the final release.  If you want to build the latest development version of Cytoscape, you should use **develop** branch for all sub-projects.  All of the core repositories are expected to be on the same branch at all times.

#### Core Apps branch management
Since core apps have their own release cycles, they have different branching scheme.  Usually, features are developed in feature branches, and there is only one common branch called **master**.  Head of the master branch is always the latest development version of the core app.    

Because the apps track **master** rather than _develop_, the two columns of the command table are never in step with each other, and that is expected.  One thing to know about ```branch-apps```: it requires all of the core apps to be on the same branch to begin with, since that shared branch is what the new one is cut from — run ```switch-apps``` first if they have drifted apart.

### Step 1: Clone the Main Project
1. Install required tools: JDK, Maven, and Git. On some systems, these may be preinstalled - you can use those versions if they are relatively recent. Released versions of Cytoscape use
JDKs from [Eclipse Adoptium](https://adoptium.net/).

2. Add JDK, Maven, and Git to your system PATH if necessary. On some platforms, this is done automatically on installation - try running mvn, git, or java at a command line to check this.

 If you need to add tools to the PATH, the steps you should follow vary by operating system. On Windows, this can be done in the Environment Variables dialog - open the file browser, right click on "Computer" or "This PC", select "Advanced system settings", then click "Environment Variables" - you should be able to edit the the PATH by selecting the PATH environment variable and clicking Edit. On Mac or Linux, you would need to edit the .profile, .bashrc or .zshrc file in your home directory to set environment variables, depending on the type of shell you are using. It may be helpful to add the following to .profile and set the PATH in .bashrc so that all shells will read the same values:

 ```
 if [ -f ~/.bashrc ]; then
   source ~/.bashrc
 fi
 ```

 Then, in .bashrc, add something like the following:

 ```
 export PATH=/path/to/java/bin:/path/to/maven/bin:/path/to/git/bin:$PATH
 ```

 Use the directories where binaries are located, as this will ensure that the command line knows where to find them.

3. Set the JAVA_HOME environment variable to the JDK 17 installation directory. This is only necessary if you have multiple versions of Java installed - if JDK 17 is the only one, Cytoscape will be able to automatically find it without the need for an environment variable. To do this, follow the same instructions as above, but for JAVA_HOME instead of PATH. On Windows, you may have to click the "Add..." button under System Variables if JAVA_HOME does not already exist. On Mac/Linux, you would add an additional line to .bashrc (or .profile if you set environment variables there) like the following.

 ```
 export JAVA_HOME=/path/to/java
 ```

 On Mac, you can use ```/usr/libexec/java_home -V``` to show the paths of all installed JVMs.

4. MAVEN_HOME, and M2_HOME to your environment variables. On some platforms, this is done automatically on installation. These are environment variables that can be set using the same methods as JAVA_HOME and PATH, and should point to the Maven installation directory (example: `/path/to/Maven/apache-maven-3.6.3`. Use the relevant echo command to test these (example: `echo $MAVEN_HOME` for Ubuntu/Mac or `echo %MAVEN_HOME%` on Windows). 

5. If you are using Git for the first time, you need to set your name and e-mail address. To do this, use the following commands:
 ```
 git config user.name "Your Name"
 git config user.email yourname@yourname.com
 ```
 Substitute your actual name and e-mail address in the commands.

6. Generate an ssh key and set it up on GitHub. To do this, you will first need to be added to the Cytoscape GitHub project by one of the core developers. Then, click the arrow in the top-right hand corner of any GitHub page and choose "Settings". Click "SSH and GPG keys", then follow the instructions on the linked guide to generate an SSH key on your particular operating system. Once you have generated an SSH key, return to the original "SSH and GPG keys" page and add the generated key using "New SSH key".

7. Clone this repository: ```git clone https://github.com/cytoscape/cytoscape.git```
    - Now you can see a new directory named **cytoscape**:

```
cytoscape     <-- your working directory from here on
├── README.md
├── cy.sh
└── pom.xml
```

Cloning this repository is always the first step: the _cy_ script lives in it, so you need it on disk before you can run any of the commands below. Everything from here on is run from inside this directory.

### Step 2: Clone the Sub Projects
The core of Cytoscape is split across six more repositories — _api_, _impl_, _parent_, _support_, _gui-distribution_ and _app-developer_. They are not submodules of this one, so cloning this repository does not bring them with it. The _cy_ script clones them for you.

1. _cd_ into the directory you cloned in Step 1: ```cd ./cytoscape```
1. Run: ```./cy.sh pull```
1. The six sub projects are now folders inside it, alongside _pom.xml_:

```
cytoscape     <-- your working directory
├── README.md
├── cy.sh
├── api
├── app-developer
├── gui-distribution
├── impl
├── parent
├── pom.xml
└── support
```

Run ```./cy.sh pull``` again whenever you want to bring everything up to date: it clones whatever is missing and pulls whatever is already there, so it is both the setup command and the update command.

It works over HTTPS or SSH — the protocol is taken from the ```origin``` remote of the clone you are in, so the sub projects are fetched the same way you fetched this one. Set ```CY_GIT_URL``` to override that, e.g. ```CY_GIT_URL=git@github.com:my-fork/ ./cy.sh pull```.

The sub projects must sit directly inside your clone, as shown above. The Maven build and the other _cy_ commands all expect that layout, and ```pull``` is what produces it.

#### Optional: Working on a different branch
```pull``` brings the sub projects down on **develop**, and the core repositories are expected to stay on the same branch as each other. Do not check out a branch in one repository by hand: use the **cy** script so that all of them move together.

**Which command you need depends on whether the branch already exists.**

| The branch you want to work on... | Command | What it does |
| --- | --- | --- |
| **already exists** — _develop_, _master_, a release branch someone else created | ```./cy.sh switch BRANCH_NAME``` | Checks that existing branch out in every repository. |
| **does not exist yet** — you are the one starting it | ```./cy.sh branch NEW_BRANCH_NAME``` | Creates the new branch in every repository, branching off the branch they are all currently on, and leaves them checked out on it. |

The two are not interchangeable, and each one refuses the other's job:

* ```switch``` only checks out a branch that is already there — it cannot create one. Pointed at a name that does not exist, it reports that it could not check it out and carries on to the next repository, which can leave the project half switched.
* ```branch``` only creates. If the name is already taken in any repository it stops and creates nothing, rather than quietly reusing what is there.

So to start new work, first ```switch``` to the branch you want to base it on (if you are not already on it), then ```branch``` to create the new one from that point. For example, to start a 3.11 release branch from _develop_:

```
./cy.sh switch develop
./cy.sh branch release/3.11
```

```branch``` validates every repository before it changes anything, and stops without creating a single branch if a sub project is missing, if they are not all on the same branch, or if the new name is already taken. See [Creating a new branch](#creating-a-new-branch) for the full details.

Build from whichever branch you end up on — the steps below are the same either way.

### Step 3: Building Cytoscape
1. Stay in your clone — the directory containing _pom.xml_ and the sub projects
1. Run: ```./cy.sh build```
1. Have a coffee break...  It depends on your machine specification and internet connection speed, but will take 3-60 minutes.  When you build Cytoscape for the first time, it will take a long time because maven downloads all dependencies from the remote server.

Use ```./cy.sh build``` rather than calling Maven yourself. A plain ```mvn install``` from the top **fails on a first build** — see the note about build order above — whereas ```cy.sh build``` runs the sub-builds in the order that works and passes the right flags. It is safe to re-run at any time; it works the same on a fresh clone and on an up-to-date one.

If you do want to drive Maven directly, the equivalent of the final step is ```mvn -fae install -U -DskipTests```, but it only works once ```api```, ```support``` and ```impl``` have already been installed.

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

### Step 5: Continue with Eclipse project steps
If you are developing in Eclipse, continue to set up with [these steps](https://github.com/cytoscape/cytoscape/wiki/Importing-Git-Repos-in-Eclipse).
You can also configure Eclipse to [debug Cytoscape](https://github.com/cytoscape/cytoscape/wiki/Launching-Cytoscape-from-Eclipse).

----

## New from 3.10.0: Core Apps
___Core Apps___ are Cytoscape Apps that are included with the Cytoscape distribution.  They are located in their own separate GitHub repositories. Cytoscape depends on the latest version of each core app deployed to the Nexus repository, so you don't need to build core apps to build Cytoscape core.

As of Cytoscape 3.10.0 (March 2023), Cytoscape core distribution comes with the following core apps:

https://github.com/cytoscape/cytoscape-gui-distribution/blob/develop/assembly/pom.xml#L447


### Optional Projects
This repository contains sample code for app developers and it will not be included in the core distribution.  **You don't have to build this repository if you just want to build the Cytoscape Core distribution.**

* [samples](https://github.com/cytoscape/cytoscape-samples)

----

## Building Core Apps
All of the core apps are maintained in their own repository and if you want to try the latest version of the core app, you need to build them separately.

### Step 1: Checking out core apps
From your clone of the main project, ```./cy.sh pull-apps``` will check out every core app into the ```apps``` subdirectory. Each is hosted in its own GitHub repository, and changes can be committed directly to each directory.

Run it again any time to bring the apps up to date — like ```pull``` for the core repositories, it clones whatever is missing and pulls whatever is already there, so it is both the setup command and the update command.  All of the core apps are hosted under this org account:

* [Cytoscape Consortium GitHub Repository](https://github.com/cytoscape)

### Step 2: Building Core Apps
All of core apps are independent to each other and there is no dependency among those.  To build the apps, you can just _cd_ to the app's directory and run:

```
mvn clean install -U
```

to build the latest version.  You can also use the following command from top-level directory to build all core apps:

```
./cy.sh build-apps
```

**Use JDK 11 for this.** The core apps do not build on JDK 16 or newer — see [Which JDK](#which-jdk). If you have been building the core, you are on JDK 17 and will need to switch for this step:

```
JAVA_HOME=$(/usr/libexec/java_home -v 11) ./cy.sh build-apps
```

This command simply runs ```mvn clean install``` for each core app directory.

### Step 3: Install the new build
To test changes, install the JAR from the Cytoscape menu at Apps > App Store > Install Apps from File, or copy to the ```~/CytoscapeConfiguration/3/apps/installed``` directory.

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
    - repositories ( You can just copy these from an existing pom, like CyREST: https://github.com/cytoscape/cyREST/blob/3e1ef7fcf867dbcd80749618c4134173e02688e9/pom.xml#L35)
    - distributionManagement ( You can just copy these from an existing pom, like CyREST: https://github.com/cytoscape/cyREST/blob/3e1ef7fcf867dbcd80749618c4134173e02688e9/pom.xml#L35)
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

## Step 6: Setup Travis CI job
You can set up a Travis CI.  You need to be registered with Travis to setup a new job, and have a username and password for https://nrnb-nexus.ucsd.edu/ to deploy to nexus.

* Go to https://travis-ci.com/github/cytoscape
* Find and select the relevant repository in the repositories list
* Follow the steps to enable Travis CI for that repository
* You can use the Travis CI setup as well as the ```.travis.yml``` from the cytoscape build as a guide to building and deploying your app.

## Step 7: Push your changes to the remote repository
If you have set up Travis CI, it will automatically deploy changes to nexus.

If you haven't set up Travis CI, you will need to deploy manually follow the steps below:

Add the code below to your maven settings.xml file (location varies between OS's, check the Maven documentation for location).

   If you use brew to install maven, then the settings file should be in: ```/usr/local/Cellar/maven/<version>/libexec/conf```

```
<servers>
  <server>
    <id>releases</id>
    <username>deployment</username>
    <password>deploy</password>
  </server>
  <server>
    <id>snapshots</id>
    <username>deployment</username>
    <password>deploy</password>
  </server>
  <server>
    <id>thirdparty</id>
    <username>deployment</username>
    <password>deploy</password>
  </server>
  <server>
    <id>cytoscape_releases</id>
    <username>cytoscape_deployer</username>
    <password>turtlesallthewaydown</password>
  </server>
  <server>
    <id>cytoscape_snapshots</id>
    <username>cytoscape_deployer</username>
    <password>turtlesallthewaydown</password>
  </server>
  <server>
    <id>cytoscape_thirdparty</id>
    <username>cytoscape_deployer</username>
    <password>turtlesallthewaydown</password>
  </server>
</servers>
```
And then run:
```
mvn clean deploy
```

This step requires that you have appropriate credentials for https://nrnb-nexus.ucsd.edu/ and that they are added to your Maven [settings.xml](https://maven.apache.org/settings.html#servers) file.


## Step 8: Add your app to _assembly_ pom
Once the new app is in Nexus repository, you can use it from the core.  Add your app as a new dependency in **cytoscape-gui-distribution/assembly/pom.xml**:

* https://github.com/cytoscape/cytoscape-gui-distribution/blob/develop/assembly/pom.xml

It should be a new entry under _**artifactItems**_ section.

## Step 9: Build Cytoscape
Now you can build Cytoscape.  Run the new version and make sure your app is part of the distribution.

## Step 10: Add the new core app in js4cytoscape repo
The core apps list is hard-coded in js4cytoscape repo develop branch for the use of app manager. The new core app should be added at the end of this [list](https://github.com/cytoscape/js4cytoscape/blob/29b12daf67adb5bb1f4d45042bc0c7803260abe3/packages/js4cytoscape/appmanager.html#L167).

## Step 11: Bump up the app's version number
Change the version number of your app to 3.x.x.


## Step 12: Update description on the App Store
Update app description on App store (https://apps.cytoscape.org/) page.

* Add 'core app' tag in categories.

* Add the following paragraph on the app page:
```
Note: This app is pre-installed in Cytoscape 3.x or later
You should not download this app from the store yourself. If Cytoscape detects that a new version is available in the App Store, it will notify you as it's starting and will give you a chance to download it directly.

You should not download this app if your Cytoscape is earlier than 3.x. An earlier version of this app is pre-installed in your Cytoscape, and this app cannot replace it.
```
# Rules for updating Core Apps

## 1. Keep _master_ branch always releasable
If you are working on a new features, create feature branch from the master.  You can push the changes anytime to the feature branch.  Once it's ready and if you want to release a new version to the store, merge it back to the master and tag it.  For example, if you are working on _feature/a_ branch and want to release version 2.1.0, do the following:

```
git checkout master
git merge feature/a
mvn clean install -U (To test the build)
git tag 2.1.0
git push --tags
```

## 2. Tagged version should be deployed to the app store
Once you tagged your commit, submit the new JAR to the store

## 3. Update _gui-distribution/assembly/pom.xml_
Core apps can be released any time, but if you want to include the new version of your core app in the next release of Cytoscape, you have to change the dependency in the pom file above.

The dependency section starts from [here](https://github.com/cytoscape/cytoscape-gui-distribution/blob/develop/assembly/pom.xml#L202).

----
## Misc. Instruction for core developers

### Creating a new branch
At some time, it may be necessary to create a new branch for a Cytoscape release. This is the case if it is necessary to carry on multiple threads of development simultaneously (i.e. 3.5 is in a release candidate stage, but we want to continue development on features destined for 3.6.)  

In this case, you will need to create a new branch for each repository in Cytoscape based on the current working branch. The **cy** script does this for all of them at once. Change to the root directory of the checked-out code (make sure all changes are committed beforehand) and run:

```
cy branch 3.x.x
```

Substitute the desired name of the new branch for "3.x.x" (the version number of the release is recommended). This creates the branch in the root repository and in every sub-repository (api, app-developer, gui-distribution, impl, parent, and support), branching each one off the branch it is currently on, and leaves them all checked out on the new branch.

Before creating anything, the command validates the whole project and **aborts without making a single change** if:

* any sub-repository is not present as a folder in the current directory — run ```cy pull``` first to set up the local repositories;
* the repositories are not all on the same branch — run ```cy switch BRANCH_NAME``` first to line them up;
* a branch of the new name already exists in any repository, locally or on _origin_.

The new branches are **local only**. To push them all to GitHub:

```
cy run-all "git push -u origin 3.x.x"
```

Note that plain ```cy push``` cannot be used for this first push: it runs ```git push -u origin``` with no branch name, which fails for a branch that does not yet have an upstream. Once the branches have been pushed once, ```cy push``` works normally.

Your new branch has been created and pushed.

### Choosing a Branch
If your Cytoscape project has Git branches, you can switch branches easily with **cy** script.  Simply go to the parent level  *cytoscape* directory and type:

```
cy switch BRANCH_NAME
```

where **BRANCH_NAME** is the name of the branch you want to switch.  All Cytoscape subprojects are following git-flow style branching scheme.  *Master* is used only for releases, and *develop* is the latest development branch.

To *create* a new branch across all of the subprojects rather than switch to an existing one, use ```cy branch NEW_BRANCH_NAME``` instead — see [Creating a new branch](#creating-a-new-branch) above.

### Managing Nested Git Repos with a GUI
The Cytoscape project is organized as a nested set of Git repositories. This provides for modularity and flexibility, but at the cost of greater complexity. The script at the parent level repository helps to manage ```git pull``` commands, but beyond that it can be challenging to manage the state of each repository. Fortunately, there are GUIs that can help:

* [GitHub Desktop](https://desktop.github.com/)
* Eclipse: [Setup instructions](../../wiki/Managing-Cytoscape-Git-Repos-in-Eclipse)
* [SourceTree](https://www.sourcetreeapp.com/)

----

# Building a New Release
This section is for core developers only.

## Introduction
These are general instructions on how to build a Cytoscape release locally. The Cytoscape Build Server (https://cytoscape-builds.ucsd.edu/) contains scripts to automate this process and should be used to generate releases. However, the general steps of the release build are outlined here to provide additional reference material.

## Updating version numbers
When preparing to make a new release (or start a new development branch), it is necessary to update the version numbers in Cytoscape to reflect this release. This is typically done at the release-candidate state, and also when updating the development branch for the next development version.

As the version number is repeated many times throughout the Cytoscape build configuration, it is easiest to make use of the Maven versions plugin to help make this change. To do this, change to the root directory in the Cytoscape directory structure. Then, run the following command:

```
mvn versions:set -DnewVersion=3.x.x
```

Substitute the desired version for "3.x.x" - if you are starting a new development branch, you would use a version number ending in -SNAPSHOT (i.e. 3.6.0-SNAPSHOT).

Though this will update most instances of the version number, it doesn't get them all. **You will need to manually update the version number in the following places:**

* cytoscape.sh (in gui-distribution/assembly/src/main/bin)
* cytoscape.bat (in gui-distribution/assembly/src/main/bin)
* parent/pom.xml (look for taglets)
* pom.xml in app-developer, gui-distribution, impl, support (look for properties tag)
* pom.xml in src/main/resources/archetype-resources/pom.xml (for each archetype subdirectory in support/archetypes)
* pom.xml in event-impl/it, model-impl/it, model-impl/performance, session-impl/impl, session-impl/integration-test, viewmodel-impl/it, vizmap-impl/it, work-swing-impl/it

You can edit by hand, or use grep/sed to update these numbers.  Push and commit all changes to GitHub when you are done updating the version numbers.

#### Tip: To find all instances of a version, run ```grep -ri "3.X.X-SNAPSHOT" .``` from the parent directory.

## Releasing unreleased updates to core apps
Typically, updates to core apps will be released separately from the Cytoscape core development cycle, and the development branch will be updated to use any new updates as they are released. However, if a core app depends on an unreleased API in the development version of Cytoscape, this won't work. In that case, we have to release the unreleased core app when Cytoscape is being released.

After updating the version numbers to a release version (at the release candidate stage), check the gui-distribution/assembly/pom.xml to see if there are any -SNAPSHOT dependencies in core apps. If so, you will want to update these core apps to use non-SNAPSHOT API dependencies and a non-SNAPSHOT version number. Then, commit/push your changes and do an ```mvn deploy``` to deploy the new version. Finally, update the gui-distribution/assembly/pom.xml core app dependency to the new (non-SNAPSHOT) version number.  The updated app should be released to the App Store right before Cytoscape is released (take care to submit the same version as was deployed to Nexus - rebuilding will result in a non-matching checksum).

## Deploying to Nexus
When we are ready to release, we need to deploy artifacts for each bundle to Nexus. To build and deploy all artifacts, run the following command from the Cytoscape top-level directory:

```
mvn clean deploy
```
Note that **you will need to configure the Nexus server in ~/.m2/settings.xml before doing this.** Deploying to Nexus will always rebuild Cytoscape, so each deployment will have a different timestamp/SHA hash than the last deployment (even if nothing has changed). The Cytoscape Project POM does not deploy and will throw an error. This is expected.

## Building Installers
To build installers, we use a proprietary tool called install4j. This can be downloaded [here](https://www.ej-technologies.com/download/install4j/files) - the license key should be available to core team members. After installing install4j, you may need to update the path in the build configuration, which is set in the install4j.executable property in gui-distribution/packaging/pom.xml. The default value will work with the default install path on a Mac - this will need to be changed if you are building on a Windows/Linux system or install to a non-default path on Mac (it is advised not to commit changes to this property, as it is intended to be set locally).

To build installers for the most recently-built code, run the following command from the ```gui-distribution/packaging``` directory:

```
mvn install
```
### How to build release

These are general instructions on how to build a Cytoscape release locally. The Cytoscape Build Server (https://cytoscape-builds.ucsd.edu/) contains scripts to automate this process and should be used to generate releases. However, the general steps of the release build are outlined here for reference.

1. Make sure you have already installed install4j
1. CD to Cytoscape project's top directory
1. Switch to release branch: `cy switch release/3.x.x`
1. Run `cy pull` to synchronize local repository to remote
1. Run `mvn clean install -U` to make sure you can build all bundles without problems
1. CD to `gui-distribution/packaging`
1. Run `mvn clean install -U`
1. CD to `target/install4j` and check you have installers for each platform
1. (Optional you need Apple developer account and Mac to do this!) CD to `gui-distribution/packaging` and run `sign-dmg.sh 'your account'` and check you have signed dmg in `signed` directory

Alternative way to build test releases is updating this shell script:

* https://github.com/cytoscape/cytoscape-scripts/blob/develop/deploy_installers.sh

This script contains machine-specific hard-coded values, and you need to understand and modify the code to run on your machine.

### Signing Mac Installer

After the installers are built, you will need to sign the dmg installers (found in the ```gui-distribution/packaging/target/media``` directory) on a Mac machine with Xcode installed as well as the Mac App Store certificate 'Developer ID Application'. To sign the Mac DMGs run the following code, substituting a valid Cytoscape Mac developer ID and password, as well as the appropriate VERSION and UNIQUE_NOTARIZATION ID:

```
xcrun notarytool submit --apple-id myappleid@blah.com --team-id {CY_TEAM_ID} --password "app-specific-password" ./Cytoscape_{VERSION}_macos_aarch64.dmg
```

The notarization command uses an "App Specific Password"
* To create an app specific password go to https://account.apple.com/
* Go to "Sign-In and Security"
* Select "App Specific Password"
* Generate a password for "notarytool"
* You must be logged into the Mac used to notarize with the same apple account.

After this is done, draft a github release at https://github.com/cytoscape/cytoscape/releases and upload the built installers. If you are building a full release, also upload the swing-app-api JAR (in api/swing-app/api/target under the Cytoscape build root) and the API Javadocs (in app-developer/target/API).


### Tagging a Release
After merging a release into the master branch, it should be tagged with the release version number. To do this, we start in the Cytoscape build root with master checked out in all projects and execute the following commands:
```
git tag x.y.z
git push origin x.y.z
```
(substitute the version number for x.y.z)

This will tag the repository and push it to GitHub. Then, repeat this for all the sub-repositories (api, app-developer, gui-distribution, impl, parent, and support).

### Finalizing a Release
There are a few other steps that need to be completed when building a release. These should be done at the very end of the process (i.e., right before sending out the announcement and updating the website/release notes).

1. Tag the manual to correspond with the new release. The manual is now a GitHub repository (located [here](https://github.com/cytoscape/cytoscape-manual)), and tagging it will create a new version of the document on ReadTheDocs. This is referenced by the Cytoscape application (using its internal version to determine the URL) - when tagging, the version number should not include any prefix or suffix and should always have three digits and two decimal places (so 3.6 should be "3.6.0").
1. Make release and create tag in the process: https://github.com/cytoscape/cytoscape/releases.
1. Update the Cytoscape.org web site: create a new release notes html file and update: [releasenotes.html](https://github.com/cytoscape/cytoscape.github.com/blob/master/releasenotes.html), [download.html](https://github.com/cytoscape/cytoscape.github.com/blob/master/download.html), [download.js](https://github.com/cytoscape/cytoscape.github.com/blob/master/js/download.js). [setup_page.js](https://github.com/cytoscape/cytoscape.github.com/blob/master/js/setup_page.js), [roadmap.html](https://github.com/cytoscape/cytoscape.github.com/blob/master/roadmap.html).
1. Update the [system-checker scripts](https://github.com/cytoscape/cytoscape-admin-scripts/tree/master/system-checker): OS, Java, and Cytoscape versions.
1. Update the version number in news.html, and add the announcement for a new release. This file is located at
http://chianti.ucsd.edu/cytoscape-news/news.html (in web.ideker.ucsd.edu:/var/www/vhosts/chianti.ucsd.edu/cytoscape-news/news.html).
1. Send an announcement e-mail to cytoscape-helpdesk, cytoscape-app-dev and cytoscape-announce.



### Updating core apps
When a core app update is ready to be pushed to the App Store and the development version of Cytoscape, there are a few steps you need to follow.

1. Update the version number in the App's pom.xml file. The version number must be a non-SNAPSHOT version, and the core app must not depend on SNAPSHOT APIs.
1. Merge your changes to the main repo.
1. Clone the repository and make sure it compiles.
1. Tag this place in the repository with the version. ```git tag MAJOR.MINOR.PATCH```. Push this tag to the remote repository using ```git push origin MAJOR.MINOR.PATCH``` to push the current this tag or ```git push origin --tags``` to push all tags in the local repo.
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
* When executing the `cy` script, be sure that your current path contains no blanks. The `cy` script's path parser does not understand blanks.
* You can follow the Git SSH instructions to create your SSH key, but when you start the SSH agent, use `eval $(ssh-agent)` instead of `eval 'ssh-agent' -s`.
* When running `cy pull`, if you get "flags: FATAL unable to determine getopt version" somewhere in the output, you must be sure to put `getopt` in your PATH. The default location for `getopt` is `C:\Program Files (x86)\GnuWin32\bin`.


### Notes for All Developers
* Note that the `cy` script always works on the clone it is run from — the sub projects are cloned into that directory. There is no option to install them somewhere else; clone this repository where you want the projects to live.
* Be sure you have installed Java JDK, not Java JRE.
* If you are developing on a virtual machine, be sure to configure around 8GB RAM and 50GB disk.
* To create a Cytoscape project in Eclipse (once you have run `cy pull`), select File | Import, and then select Maven | Existing Maven Projects. Browse to your clone of this repository, and note that all pom.xml files are found. To finish the import, wait for all projects to be created and compiled. This may take several minutes.
* To add all Cytoscape sources, use the Source tab in the Debug Configurations dialog, click the Add button, choose the Java Project container, and select all projects.
* To edit-compile-run, make your changes in the project you're working in. From Eclipse, you can Run As ... Maven Install. Eclipse will build the .class files automatically, so Maven's job is to create the .jar and promote it to private Maven repository. An unresolved compile issue will show in the Cytoscape console window when you run ... Maven doesn't complain, and Eclipse complains visually. Alternative: in Git Bash, set pwd to project directory (e.g., welcome-impl) and do `mvn clean install -U`.

