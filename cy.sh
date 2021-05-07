#!/bin/bash
#
# @(#) cy version 3.0.1 11/7/2017
#
#  USAGE:
#    init
#
# DESCRIPTION:
#   Cytoscape 3 repository management utility.
#   This script is only for core developers.
#
# Requirments:
#   - git
#
# By Keiichiro Ono (kono at ucsd edu)
#
###############################################################################

# Command Name
CMDNAME=$(basename $0)

# Error Message
ERROR_MESSAGE="Usage: $CMDNAME [-h] [action]"

# Help
HELP='Cytoscape build helper script'

# Git base URL
BASE_URL='git@github.com:cytoscape/cytoscape-'

# Core Apps URL
APP_URL='git@github.com:cytoscape/'

# Cytoscape repository names
REPOSITORIES=(. parent api impl support gui-distribution app-developer)

# List of Core Apps
CORE_APPS=(amatreader analyzer biopax core-apps-meta cyREST \
json idmapper network-merge opencl-cycl opencl-layout \
psi-mi sbml webservice-psicquic-client webservice-biomart-client \
cx diffusion cy-ndex-2 copycat-layout cyBrowser file-transfer)

#######################################
# Handling command-line arguments     #
#######################################
while getopts 'hd:' OPT
do
  case $OPT in
    h)  FLG_H=1
        echo "$HELP: $ERROR_MESSAGE"
        exit 0
        ;;
    ?)  echo $ERROR_MESSAGE 1>&2
        exit 1 ;;
  esac
done

shift $(($OPTIND - 1))

COMMAND=$1
TARGET_DIR=$2

if [[ -z $COMMAND ]]; then
  echo "COMMAND is required. $ERROR_MESSAGE" 1>&2
  exit 1
fi


###############################################################################
# Functrions
###############################################################################
function reset {
  echo "This command resets all of your local changes!"
  confirm

  for REPO in "${REPOSITORIES[@]}"; do
    echo "\n - Resetting local changes: $REPO"
    pushd $REPO
    git clean -f -d
    git reset --hard
    popd
  done
}

function pull {
	echo "------------------------------------------------------------------------"
  for REPO in "${REPOSITORIES[@]}"; do
    pushd $REPO > /dev/null
    echo "Downloading changes from upstream: $REPO"
    git pull
    popd > /dev/null
		echo "------------------------------------------------------------------------"
  done
}

function push {
  echo "- Sending all local commits to upstream..."
  for REPO in "${REPOSITORIES[@]}"; do
    pushd $REPO
    git push -u origin
    popd
  done
}

function status {
	echo "------------------------------------------------------------------------"
  for REPO in "${REPOSITORIES[@]}"; do
    pushd $REPO > /dev/null || { echo Could not find subproject; exit 1; }
    echo "- $REPO:"
		echo
    git status
    popd > /dev/null
		echo "------------------------------------------------------------------------"
  done

}

function switch {
  TARGET="${TARGET_DIR}"
  if [[ -z $TARGET ]]; then
    echo "Branch name is required: cy switch BRANCH_NAME" 1>&2
    exit 1
  fi

  for REPO in "${REPOSITORIES[@]}"; do
    echo "\n - Switching to ${TARGET}: $REPO"
    pushd $REPO || { echo Could not find subproject; exit 1; }

    # Switch
    git checkout $TARGET || { echo Could not checkout branch $TARGET; }
    popd
  done
}

# Not finished yet.
function resetAll {
  git checkout master
  git reset --hard $(git BRANCH -av | grep "remotes/origin/master" | awk '{ print $2 }')
  git clean -d -f

  git checkout develop
  git reset --hard $(git BRANCH -av | grep "remotes/origin/develop" | awk '{ print $2 }')
  git clean -d -f

  git checkout $BRANCH
  git reset --hard $(git BRANCH -av | grep "remotes/origin/$BRANCH" | awk '{ print $2 }')
  git clean -d -f
  git BRANCH -avv
}

function confirm {
  printf 'Do you really want to continue?  Type [yes] to proceed: '
  read answer

  if [[ $answer != 'yes' || -z $answer ]]; then
    echo "Abort\n"
    exit 0
  fi
}

#################################################################################
#
# Project initializer.
#
#   This function clones top-level pom.xml and then clones all sub-projects
#   into the top-level Cytoscape directory.
#
#   Use -d option to specify directory.  Otherwise, it will create new project
#   in the current directory.
#
#################################################################################
function init {
  # By default, clone everything in current directory.

  echo "Target directory = $TARGET_DIR"

  if [[ -z "$TARGET_DIR" ]]; then
    TARGET_DIR=$(pwd)
  elif ! [ -e "$TARGET_DIR" ]; then
    echo "No such dir: $TARGET_DIR"
    exit 1
  fi

  echo "Cytoscape project will be cloned to: ${TARGET_DIR}"

  cd $TARGET_DIR || { echo Could not find target directory: $TARGET_DIR; exit 1; }

  # Clone cy
  let LENGTH=${#BASE_URL}-1

  CY3REPO_URL=$(echo $BASE_URL | cut -c 1-$LENGTH)
  git clone "${CY3REPO_URL}.git" || { echo Could not clone remote repository: $TARGET_DIR; exit 1; }

  cd cytoscape

  for REPO in "${REPOSITORIES[@]}"; do
    if [ ${REPO} != . ]
    then
        REPO_URL="$BASE_URL$REPO.git"
        echo "Cloning: $REPO (URI = $REPO_URL)"
        git clone $REPO_URL $REPO || { echo Could not clone remote repository: $TARGET_DIR; exit 1; }
        pushd $REPO
        # git checkout master
        # git flow init -d
        git checkout develop
        popd
    fi
  done

  echo "\n\n - Finished: here is the current status:\n"
  status
}


# Utility command to clone all Core Apps as independent repositories.
function apps {
  mkdir ./apps
  cd apps

  for app in "${CORE_APPS[@]}"; do
    echo "\n - Cloning $app"
    REPO_URL="$APP_URL$app.git"
    git clone $REPO_URL || { echo Could not clone remote repository: $REPO_URL; exit 1; }
  done
  cd ..
}

function update-apps {
  cd apps

  for app in "${CORE_APPS[@]}"; do
    cd $app
    echo "- Pulling changes for $app"
    git pull || { echo Could not pull changes: $REPO_URL; exit 1; }
    cd ..
  done
  cd ..
}

function validate-apps {
  cd apps

  for app in "${CORE_APPS[@]}"; do
    cd $app
    (mvn validate | grep Building \
    | awk '{for (i=3; i<NF; i++) printf $i " "; print $NF}') || { echo Could not validate: $REPO_URL; exit 1; }
    cd ..
  done
  cd ..
}

function build-apps {
  cd apps
  echo "Core apps: $CORE_APPS"

  for app in "${CORE_APPS[@]}"; do
    cd $app
    echo "- Building $app"
    mvn clean install || { echo Could not pull changes: $REPO_URL; exit 1; }
    cd ..
  done
  cd ..
}

function switch-apps {
  TARGET="${TARGET_DIR}"
  if [[ -z $TARGET ]]; then
    echo "Branch name is required: cy switch-apps BRANCH_NAME" 1>&2
    exit 1
  fi

  cd apps
  for app in "${CORE_APPS[@]}"; do
    cd $app
    echo "- Switching to ${TARGET}: $app"
    git checkout $TARGET || { echo Could not checkout branch $TARGET; }
    cd ..
  done
  cd ..
}

function init-all {
  cur_dir=`pwd`

  init
  cd ./cytoscape
  mvn clean install || { echo Failed to build Cytoscape; }

  cd ${cur_dir}
  apps
  cd ./apps
  build-apps
  cd -
}

function run-all {
  echo "------------------------------------------------------------------------"
  echo "Executing command: $TARGET_DIR"
  for REPO in "${REPOSITORIES[@]}"; do
    echo "--in $REPO"
		pushd $REPO > /dev/null
		$TARGET_DIR
    popd > /dev/null
		echo "------------------------------------------------------------------------"
  done
}

###############################################################################
# Main workflow
###############################################################################

# Save current directory location
START_DIR=$(pwd)

case $COMMAND in
  init )    init ;;
  init-all )    init-all ;;
  reset )    reset ;;
  push )    push ;;
  pull )    pull ;;
  switch )  switch ;;
  status )  status ;;
  apps )    apps ;;
  pull-apps )    update-apps ;;
  build-apps )    build-apps ;;
  switch-apps )  switch-apps ;;
  validate-apps )  validate-apps ;;
  run-all ) run-all ;;
  * )      echo "Invalid command $COMMAND: $ERROR_MESSAGE"
          exit 1;;
esac
