Cytoscape
=========

Entire Cytoscape 3 projects.

# How to build Cytoscape 3

## Requirments

You need the following software to build Cytoscape 3:

* JDK 6 or later
* Maven 3
* Git

If you are a core developer, you also need:

* [git-flow](https://github.com/nvie/gitflow)
  * Cytoscape 3 Core Projects are managed by following this branching rule:
    * [A successful Git branching model](http://nvie.com/posts/a-successful-git-branching-model/)

If you want to clone all sub-projects at once, you also need:

* [Repo](http://code.google.com/p/git-repo/)

## Cytoscape 3 Project Structure
* Parent
* API
* Impl
* Support
* GUI-Distribution
* Headless-Distribution
* APP-Developer
* Samples

## Clone All Cytoscape 3 Sub-Projects

```
repo init -u git@github.com:cytoscape/cytoscape.git
repo sync
```

You can clone Cytoscape 3 project manually if you want.

## Build
From the top directory, type:
```
mvn clean install
```


