#!/bin/sh
#
# @(#) cy version 2.0.0 9/15/2015
#
#  USAGE:
#    init
#
# DESCRIPTION:
#   Cytoscape 3 repository management utility.
#   This script is only for core developers.
#
# Reqiirments:
#   - git
#   - git-flow
#
# By Keiichiro Ono (kono at ucsd edu)
#
###############################################################################

# Command Name
CMDNAME=$(basename $0)

# Error Message
ERROR_MESSAGE="Usage: $CMDNAME [-h] [action]"

# Help
HELP='Cytoscape repository management tool'

# Git base URL
BASE_URL='git@github.com:cytoscape/cytoscape-'

# Core Apps URL
APP_URL='git@github.com:cytoscape/'

NON_CORE_URL='git://github.com/cytoscape/cytoscape-'

# Cytoscape repository names
REPOSITORIES=(. parent api impl support app gui-distribution app-developer)

# List of Core Apps
CORE_APPS=(biopax command-dialog datasource-biogrid network-analyzer network-merge psi-mi sbml welcome webservice-psicquic-client webservice-biomart-client)


#######################################
# Handling command-line arguments     #
#######################################
while getopts 'hrd:' OPT
do
  case $OPT in
    r)  FLG_R=1
        echo " - Using read-only repository."
        ;;
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
    popd ..
  done
}

function pull {
  for REPO in "${REPOSITORIES[@]}"; do
    pushd $REPO
    echo "Downloading changes from upstream: $REPO"
    git pull
    popd
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
  for REPO in "${REPOSITORIES[@]}"; do
    pushd $REPO || { echo Could not find subproject; exit 1; }
    echo "\n- $REPO:"
    git status
    popd
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
        git checkout master
        git flow init -d
        git checkout develop
        if [[ $REPO == "app" ]]; then
          echo "Initializing Core App Submodules..."
          git submodule init
          git submodule update
        fi
        popd
    fi
  done

  echo "\n\n - Finished: here is the current status:\n"
  status
}


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

###############################################################################
# Main workflow
###############################################################################

# Save current directory location
START_DIR=$(pwd)

case $COMMAND in
  init )    init ;;
  reset )    reset ;;
  push )    push ;;
  pull )    pull ;;
  switch )  switch ;;
  status )  status ;;
  apps )    apps ;;

  * )      echo "Invalid command $COMMAND: $ERROR_MESSAGE"
          exit 1;;
esac
