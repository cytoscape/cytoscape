#!/bin/bash
#
# @(#) cy version 3.0.1 11/7/2017
#
#  USAGE:
#    Clone this project first, then run these from inside your clone.
#    Each command has a core-repository form and a core-app form:
#
#                            CORE REPOSITORIES     CORE APPS
#    clone and update        pull                  pull-apps
#    check out a branch      switch BRANCH         switch-apps BRANCH
#    create a new branch     branch NEW_BRANCH     branch-apps NEW_BRANCH
#    build                   build                 build-apps
#
#    Also: status, push, reset, run-all, validate-apps
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

# Git base URL.  Both of these are replaced at runtime by detect-git-urls with
# the protocol this workspace is actually using; the SSH form is only a default.
BASE_URL='git@github.com:cytoscape/cytoscape-'

# Core Apps URL
APP_URL='git@github.com:cytoscape/'

# GitHub organization hosting the project
GITHUB_ORG='cytoscape'

# Cytoscape repository names
REPOSITORIES=(. parent api impl support gui-distribution app-developer)

# Branch a freshly cloned repository is expected to be on.  The core repositories
# follow git-flow and develop on 'develop'; the core apps have their own release
# cycles and develop on 'master'.  A clone lands on whatever the remote's default
# branch happens to be, which is not always these, so 'pull' and 'pull-apps'
# check out these explicitly to keep the result deterministic.
CORE_BRANCH='develop'
APPS_BRANCH='master'

# List of Core Apps.  These are repository names under the GitHub org, which do
# not always match the app's maven artifactId - file-transfer-app builds the
# 'file-transfer' artifact, for instance.
CORE_APPS=(amatreader analyzer biopax core-apps-meta cyREST \
json idmapper network-merge opencl-cycl opencl-layout \
psi-mi sbml webservice-psicquic-client webservice-biomart-client \
cx diffusion cy-ndex-2 copycat-layout cyBrowser file-transfer-app)

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

#################################################################################
#
# Decide whether to talk to GitHub over SSH or HTTPS, and set BASE_URL/APP_URL.
#
#   Cloning fails outright if the protocol does not match how you authenticate
#   to GitHub, so it is worked out automatically rather than hard coded:
#
#     1. CY_GIT_URL, if set, is used verbatim as the organization URL prefix.
#        Use this for a fork or a mirror, e.g.
#          CY_GIT_URL=git@github.com:my-fork/ ./cy.sh pull
#     2. Otherwise, if the current directory is already a clone, the protocol of
#        its own 'origin' remote is reused, so the sub-projects are cloned the
#        same way the top-level project was.
#     3. Otherwise SSH is used if it can actually authenticate to GitHub, and
#        HTTPS if it cannot.
#
#################################################################################
function detect-git-urls {
  ORG_URL="$CY_GIT_URL"

  if [[ -z $ORG_URL ]]; then
    ORIGIN_URL=$(git config --get remote.origin.url 2>/dev/null)

    case "$ORIGIN_URL" in
      https://*|http://*)  ORG_URL="https://github.com/$GITHUB_ORG/" ;;
      git@*|ssh://*)       ORG_URL="git@github.com:$GITHUB_ORG/" ;;
    esac
  fi

  if [[ -z $ORG_URL ]]; then
    # Nothing to learn from, so ask GitHub whether SSH would work at all.
    if ssh -T -o BatchMode=yes -o ConnectTimeout=10 \
           -o StrictHostKeyChecking=accept-new git@github.com 2>&1 \
         | grep -q 'successfully authenticated'; then
      ORG_URL="git@github.com:$GITHUB_ORG/"
    else
      ORG_URL="https://github.com/$GITHUB_ORG/"
    fi
  fi

  BASE_URL="${ORG_URL}cytoscape-"
  APP_URL="$ORG_URL"
}

#################################################################################
#
# Put a freshly cloned repository on the branch this script expects.
#
#   A clone lands on whatever the remote's default branch is, and those are not
#   consistent across the Cytoscape repositories - cyREST defaults to 'develop'
#   while the other core apps default to 'master', for instance.  Cloning is
#   therefore not enough on its own: check where the clone actually landed, move
#   it if it is somewhere else, and fail loudly if the expected branch does not
#   exist at all rather than silently leaving the repository somewhere unexpected.
#
#################################################################################
function require-branch {
  CLONED_REPO="$1"
  WANT_BRANCH="$2"

  GOT_BRANCH=$(git -C "$CLONED_REPO" symbolic-ref --short -q HEAD)

  # The branch has to exist, not just be the name HEAD points at: a clone whose
  # remote HEAD is dangling reports a branch that was never actually created.
  if [[ $GOT_BRANCH == $WANT_BRANCH ]] \
     && git -C "$CLONED_REPO" show-ref --verify --quiet "refs/heads/$WANT_BRANCH"; then
    return
  fi

  if ! git -C "$CLONED_REPO" show-ref --verify --quiet "refs/remotes/origin/$WANT_BRANCH"; then
    echo "FAILED: $CLONED_REPO has no '$WANT_BRANCH' branch (it cloned onto '$GOT_BRANCH')." 1>&2
    echo "Every repository here is expected to have '$WANT_BRANCH'.  Either the repository" 1>&2
    echo "has changed its branching scheme or it no longer belongs in this script's list." 1>&2
    exit 1
  fi

  echo "  cloned onto '$GOT_BRANCH', checking out '$WANT_BRANCH'"
  git -C "$CLONED_REPO" checkout "$WANT_BRANCH" \
    || { echo "FAILED: could not check out '$WANT_BRANCH' in $CLONED_REPO" 1>&2; exit 1; }
}

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

#################################################################################
#
# Bring every repository up to date, cloning any that are not here yet.
#
#   This is the command to run in a freshly cloned top-level project to get all
#   of the sub-projects locally.  Newly cloned sub-projects are put on 'develop';
#   use 'switch' afterwards if you want them all on some other branch.
#
#################################################################################
function pull {
  detect-git-urls

	echo "------------------------------------------------------------------------"
  for REPO in "${REPOSITORIES[@]}"; do
    if [ ${REPO} != . ] && [ ! -e ${REPO}/.git ]; then
      if [ -d ${REPO} ]; then
        echo "Cannot clone $REPO: '$REPO' already exists but is not a git repository." 1>&2
        echo "Remove or rename it, then run this again." 1>&2
        exit 1
      fi

      REPO_URL="$BASE_URL$REPO.git"
      echo "Cloning missing sub-project: $REPO (URI = $REPO_URL)"
      git clone $REPO_URL $REPO || { echo Could not clone remote repository: $REPO_URL; exit 1; }
      require-branch "$REPO" "$CORE_BRANCH"
			echo "------------------------------------------------------------------------"
      continue
    fi

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

#################################################################################
#
# Create a new branch in every Cytoscape repository.
#
#   The "origin branch" is the branch currently checked out in the top-level
#   repository.  Every repository listed in REPOSITORIES must exist as a child
#   folder of the current directory and must already be on that same branch.
#
#   Everything is validated up front: if any repository is missing, is on a
#   different branch, or already has a branch of the new name, then nothing is
#   modified at all.  This command never clones anything - run 'pull' first to
#   set up the local repositories.
#
#   New branches are created locally only.  Use 'push' to publish them.
#
#################################################################################
function branch {
  NEW_BRANCH="${TARGET_DIR}"
  if [[ -z $NEW_BRANCH ]]; then
    echo "Branch name is required: cy branch NEW_BRANCH_NAME" 1>&2
    exit 1
  fi

  git check-ref-format --branch "$NEW_BRANCH" > /dev/null 2>&1 \
    || { echo "Invalid branch name: $NEW_BRANCH" 1>&2; exit 1; }

  # The origin branch is whatever the top-level repository is currently on.
  ORIGIN_BRANCH=$(git symbolic-ref --short -q HEAD)
  if [[ -z $ORIGIN_BRANCH ]]; then
    echo "Could not determine the current branch of the top-level repository (detached HEAD?)" 1>&2
    exit 1
  fi

  if [[ $NEW_BRANCH == $ORIGIN_BRANCH ]]; then
    echo "Already on branch $ORIGIN_BRANCH: nothing to do" 1>&2
    exit 1
  fi

  echo "Origin branch: $ORIGIN_BRANCH"
  echo "New branch:    $NEW_BRANCH"
  echo "------------------------------------------------------------------------"

  # Validation pass.  Nothing is modified here.
  ERRORS=0
  MISSING_REPOS=0
  MISMATCHED_REPOS=0

  for REPO in "${REPOSITORIES[@]}"; do
    if [[ ! -d $REPO ]]; then
      printf '  %-20s %s\n' "$REPO" "MISSING: no such directory"
      let ERRORS=ERRORS+1
      let MISSING_REPOS=MISSING_REPOS+1
      continue
    fi

    # Test for .git in the folder itself.  Without it, git would walk up the
    # tree and quietly report on the enclosing top-level repository instead.
    if [[ ! -e $REPO/.git ]] || ! git -C "$REPO" rev-parse --git-dir > /dev/null 2>&1; then
      printf '  %-20s %s\n' "$REPO" "NOT A GIT REPO: directory exists but is not a git repository"
      let ERRORS=ERRORS+1
      let MISSING_REPOS=MISSING_REPOS+1
      continue
    fi

    CURRENT_BRANCH=$(git -C "$REPO" symbolic-ref --short -q HEAD)

    if [[ -z $CURRENT_BRANCH ]]; then
      printf '  %-20s %s\n' "$REPO" "DETACHED HEAD: expected $ORIGIN_BRANCH"
      let ERRORS=ERRORS+1
      let MISMATCHED_REPOS=MISMATCHED_REPOS+1
      continue
    fi

    if [[ $CURRENT_BRANCH != $ORIGIN_BRANCH ]]; then
      printf '  %-20s %s\n' "$REPO" "MISMATCH: on $CURRENT_BRANCH, expected $ORIGIN_BRANCH"
      let ERRORS=ERRORS+1
      let MISMATCHED_REPOS=MISMATCHED_REPOS+1
      continue
    fi

    if git -C "$REPO" show-ref --verify --quiet "refs/heads/$NEW_BRANCH"; then
      printf '  %-20s %s\n' "$REPO" "EXISTS: $NEW_BRANCH is already a local branch"
      let ERRORS=ERRORS+1
      continue
    fi

    if git -C "$REPO" show-ref --verify --quiet "refs/remotes/origin/$NEW_BRANCH"; then
      printf '  %-20s %s\n' "$REPO" "EXISTS: origin/$NEW_BRANCH already exists"
      let ERRORS=ERRORS+1
      continue
    fi

    printf '  %-20s %s\n' "$REPO" "$CURRENT_BRANCH   OK"
  done

  echo "------------------------------------------------------------------------"

  if [[ $ERRORS -ne 0 ]]; then
    echo "FAILED: $ERRORS of ${#REPOSITORIES[@]} repositories cannot be branched." 1>&2
    echo "No branches were created." 1>&2
    if [[ $MISSING_REPOS -ne 0 ]]; then
      echo "Run './$CMDNAME pull' first to set up all of the local repositories." 1>&2
    fi
    if [[ $MISMATCHED_REPOS -ne 0 ]]; then
      echo "Run './$CMDNAME switch $ORIGIN_BRANCH' to put every repository on the origin branch." 1>&2
    fi
    exit 1
  fi

  # Creation pass.  The origin branch is named explicitly as the starting point
  # so that every repository branches off it, and not off whatever HEAD happens
  # to be.  Validation has already established that they are the same thing.
  CREATED=""

  for REPO in "${REPOSITORIES[@]}"; do
    echo " - Creating $NEW_BRANCH from $ORIGIN_BRANCH: $REPO"
    git -C "$REPO" checkout -b "$NEW_BRANCH" "$ORIGIN_BRANCH" || {
      echo "Could not create branch $NEW_BRANCH in $REPO" 1>&2
      echo "These repositories were already switched to $NEW_BRANCH:$CREATED" 1>&2
      echo "To undo, run in each of them: git checkout $ORIGIN_BRANCH && git branch -d $NEW_BRANCH" 1>&2
      exit 1
    }
    CREATED="$CREATED $REPO"
  done

  echo "------------------------------------------------------------------------"
  echo "All repositories are now on $NEW_BRANCH."
  echo "The new branches are local only.  To publish them, run:"
  echo "  ./$CMDNAME run-all \"git push -u origin $NEW_BRANCH\""
  echo "('$CMDNAME push' cannot publish them: it runs 'git push -u origin' with no"
  echo " refspec, which fails for a branch that has no upstream yet.)"
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
# Bring every core app up to date, cloning any that are not here yet.
#
#   The apps counterpart of 'pull'.  They are cloned into an 'apps' directory
#   inside your clone of the main project, one directory per app.
#
#   Core apps are left on whichever branch they clone with - normally 'master'.
#   They have their own release cycles and do not use 'develop' the way the core
#   repositories do, so nothing is checked out for you here.
#
#################################################################################
function pull-apps {
  detect-git-urls
  echo "Cloning from: $APP_URL"

  mkdir -p ./apps

	echo "------------------------------------------------------------------------"
  for app in "${CORE_APPS[@]}"; do
    if [ ! -e apps/${app}/.git ]; then
      if [ -d apps/${app} ]; then
        echo "Cannot clone $app: 'apps/$app' already exists but is not a git repository." 1>&2
        echo "Remove or rename it, then run this again." 1>&2
        exit 1
      fi

      REPO_URL="$APP_URL$app.git"
      echo "Cloning missing core app: $app (URI = $REPO_URL)"
      git clone $REPO_URL apps/$app || { echo Could not clone remote repository: $REPO_URL; exit 1; }
      require-branch "apps/$app" "$APPS_BRANCH"
			echo "------------------------------------------------------------------------"
      continue
    fi

    echo "Downloading changes from upstream: $app"
    git -C apps/$app pull
		echo "------------------------------------------------------------------------"
  done
}

# Fails unless every core app is here, so that the commands below cannot end up
# running maven in the wrong directory.
function require-apps {
  if [ ! -d apps ]; then
    echo "No 'apps' directory here: run './$CMDNAME pull-apps' first" 1>&2
    exit 1
  fi

  for app in "${CORE_APPS[@]}"; do
    if [ ! -d apps/${app} ]; then
      echo "Core app is missing: apps/$app - run './$CMDNAME pull-apps' first" 1>&2
      exit 1
    fi
  done
}

function validate-apps {
  require-apps

  for app in "${CORE_APPS[@]}"; do
    pushd apps/$app > /dev/null
    (mvn validate | grep Building \
    | awk '{for (i=3; i<NF; i++) printf $i " "; print $NF}') || { echo Could not validate: $app; exit 1; }
    popd > /dev/null
  done
}

function build-apps {
  require-apps

  for app in "${CORE_APPS[@]}"; do
    pushd apps/$app > /dev/null
    echo "- Building $app"
    mvn clean install || { echo Could not build: $app; exit 1; }
    popd > /dev/null
  done
}

function switch-apps {
  TARGET="${TARGET_DIR}"
  if [[ -z $TARGET ]]; then
    echo "Branch name is required: cy switch-apps BRANCH_NAME" 1>&2
    exit 1
  fi

  if [ ! -d apps ]; then
    echo "No 'apps' directory here: run './$CMDNAME pull-apps' first" 1>&2
    exit 1
  fi

  for app in "${CORE_APPS[@]}"; do
    echo "- Switching to ${TARGET}: $app"
    git -C apps/$app checkout $TARGET || { echo Could not checkout branch $TARGET; }
  done
}

#################################################################################
#
# Create a new branch in every core app.
#
#   The apps counterpart of 'branch'.  The core apps have no top-level project to
#   take the origin branch from, so instead they must all be on the same branch
#   already, and that shared branch is the one the new branch is cut from.
#
#   Everything is validated up front: if any app is missing, is on a different
#   branch from the others, or already has a branch of the new name, then nothing
#   is modified at all.  Run 'pull-apps' first if the apps are not here yet.
#
#   New branches are created locally only.
#
#################################################################################
function branch-apps {
  NEW_BRANCH="${TARGET_DIR}"
  if [[ -z $NEW_BRANCH ]]; then
    echo "Branch name is required: cy branch-apps NEW_BRANCH_NAME" 1>&2
    exit 1
  fi

  git check-ref-format --branch "$NEW_BRANCH" > /dev/null 2>&1 \
    || { echo "Invalid branch name: $NEW_BRANCH" 1>&2; exit 1; }

  if [ ! -d apps ]; then
    echo "No 'apps' directory here: run './$CMDNAME pull-apps' first" 1>&2
    exit 1
  fi

  echo "New branch:    $NEW_BRANCH"
  echo "------------------------------------------------------------------------"

  # Survey pass.  Work out what branch the apps are on, without changing a thing.
  ERRORS=0
  MISSING_APPS=0
  BRANCHES=""

  for app in "${CORE_APPS[@]}"; do
    if [ ! -e apps/${app}/.git ] || ! git -C "apps/$app" rev-parse --git-dir > /dev/null 2>&1; then
      printf '  %-30s %s\n' "$app" "MISSING: not a local git repository"
      let ERRORS=ERRORS+1
      let MISSING_APPS=MISSING_APPS+1
      continue
    fi

    CURRENT_BRANCH=$(git -C "apps/$app" symbolic-ref --short -q HEAD)

    if [[ -z $CURRENT_BRANCH ]]; then
      printf '  %-30s %s\n' "$app" "DETACHED HEAD"
      let ERRORS=ERRORS+1
      continue
    fi

    printf '  %-30s %s\n' "$app" "$CURRENT_BRANCH"

    case " $BRANCHES " in
      *" $CURRENT_BRANCH "*)  ;;
      *)  BRANCHES="$BRANCHES $CURRENT_BRANCH" ;;
    esac
  done

  echo "------------------------------------------------------------------------"

  if [[ $ERRORS -ne 0 ]]; then
    echo "FAILED: $ERRORS core app(s) cannot be branched.  No branches were created." 1>&2
    if [[ $MISSING_APPS -ne 0 ]]; then
      echo "Run './$CMDNAME pull-apps' first to set up all of the core apps." 1>&2
    fi
    exit 1
  fi

  # All apps must agree on one branch: that is what the new branch is cut from.
  BRANCH_LIST=($BRANCHES)
  if [[ ${#BRANCH_LIST[@]} -ne 1 ]]; then
    echo "FAILED: the core apps are on ${#BRANCH_LIST[@]} different branches:$BRANCHES" 1>&2
    echo "They must all be on the same branch before a new one can be created." 1>&2
    echo "Run './$CMDNAME switch-apps BRANCH_NAME' to line them up first." 1>&2
    echo "No branches were created." 1>&2
    exit 1
  fi

  ORIGIN_BRANCH="${BRANCH_LIST[0]}"
  echo "Origin branch: $ORIGIN_BRANCH (all core apps agree)"

  # The new name must not be taken anywhere, locally or on the remote.
  for app in "${CORE_APPS[@]}"; do
    if git -C "apps/$app" show-ref --verify --quiet "refs/heads/$NEW_BRANCH"; then
      printf '  %-30s %s\n' "$app" "EXISTS: $NEW_BRANCH is already a local branch"
      let ERRORS=ERRORS+1
    elif git -C "apps/$app" show-ref --verify --quiet "refs/remotes/origin/$NEW_BRANCH"; then
      printf '  %-30s %s\n' "$app" "EXISTS: origin/$NEW_BRANCH already exists"
      let ERRORS=ERRORS+1
    fi
  done

  if [[ $ERRORS -ne 0 ]]; then
    echo "FAILED: $NEW_BRANCH already exists in $ERRORS core app(s)." 1>&2
    echo "No branches were created." 1>&2
    exit 1
  fi

  # Creation pass.  The origin branch is named explicitly as the starting point.
  CREATED=""

  for app in "${CORE_APPS[@]}"; do
    echo " - Creating $NEW_BRANCH from $ORIGIN_BRANCH: $app"
    git -C "apps/$app" checkout -b "$NEW_BRANCH" "$ORIGIN_BRANCH" || {
      echo "Could not create branch $NEW_BRANCH in $app" 1>&2
      echo "These core apps were already switched to $NEW_BRANCH:$CREATED" 1>&2
      echo "To undo, run in each of them: git checkout $ORIGIN_BRANCH && git branch -d $NEW_BRANCH" 1>&2
      exit 1
    }
    CREATED="$CREATED $app"
  done

  echo "------------------------------------------------------------------------"
  echo "All core apps are now on $NEW_BRANCH."
  echo "The new branches are local only.  To publish them, run:"
  echo "  for a in apps/*/; do git -C \"\$a\" push -u origin $NEW_BRANCH; done"
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

###########################################
# Core building functions for newbs
###########################################

#################################################################################
#
# Build the core.  Works from a completely empty local maven repository, so it
# can be used for a first build and for every build after that.
#
#   Two things make a plain 'mvn install' from the top fail on a first build,
#   and both are worked around here:
#
#   1. api/pom.xml and impl/pom.xml both bind maven-source-plugin's 'aggregate'
#      goal, which FORKS a generate-sources lifecycle over their modules.  A
#      forked lifecycle resolves dependencies from the local repository instead
#      of the reactor, so it looks for event-api before the reactor has built it
#      and the whole build dies at module 2 of 122.  Building the inner pieces
#      first seeds the repository so those forks can resolve - which is what the
#      old 'build-core-magic' was doing, and why the README says a first build
#      needs event-api, then api, then everything.  'support' is in the sequence
#      because impl needs its cmdline, cy-commons-cli and integration-test-support
#      artifacts, which the old build-core-magic never built - which is why that
#      command could not have worked from cold either.
#      Note '-Dmaven.source.skip=true' does NOT avoid this - the fork is planned
#      before the skip is evaluated.
#
#   2. '-DskipTests', never '-Dmaven.test.skip=true'.  34 poms depend on another
#      module's test-jar, and the app-developer archetype integration tests want
#      event-api's test-jar too.  'maven.test.skip' skips compiling tests at all,
#      so those test-jars are never produced and everything needing one fails.
#      '-DskipTests' still compiles and packages tests, it just does not run them.
#
#################################################################################
function build {
  echo "- Seeding event-api (needed before the api source:aggregate fork runs)"
  (cd api/event-api && mvn install -U -DskipTests) \
    || { echo "Failed to build api/event-api" 1>&2; exit 1; }

  echo "- Building api"
  (cd api && mvn install -U -DskipTests) \
    || { echo "Failed to build api" 1>&2; exit 1; }

  # '-Darchetype.test.skip=true' only for this seeding pass.  support's archetype
  # integration tests generate a project from each template and build it, and the
  # starter-app template depends on model-impl - so support cannot be fully built
  # before impl, while impl needs support's cmdline.  Skipping just the archetype
  # ITs here breaks that cycle; the final full build runs them for real.
  echo "- Building support (archetype integration tests deferred to the full build)"
  (cd support && mvn install -U -DskipTests -Darchetype.test.skip=true) \
    || { echo "Failed to build support" 1>&2; exit 1; }

  echo "- Building impl"
  (cd impl && mvn install -U -DskipTests) \
    || { echo "Failed to build impl" 1>&2; exit 1; }

  echo "- Building everything"
  mvn -fae install -U -DskipTests
}


###############################################################################
# Main workflow
###############################################################################

# Save current directory location
START_DIR=$(pwd)

case $COMMAND in
  init )    echo "'init' has been removed: it tried to clone the main project, which you" 1>&2
            echo "must already have cloned in order to be running this script." 1>&2
            echo "Run './$CMDNAME pull' from your clone to get the sub projects instead." 1>&2
            exit 1 ;;
  apps )    echo "'apps' has been merged into 'pull-apps', which clones any missing" 1>&2
            echo "core app and updates the rest.  Run './$CMDNAME pull-apps' instead." 1>&2
            exit 1 ;;
  build-core )  echo "'build-core' has been renamed to 'build'." 1>&2
                echo "Run './$CMDNAME build' instead." 1>&2
                exit 1 ;;
  build-core-magic )  echo "'build-core-magic' has been removed: 'build' now does the staged" 1>&2
                      echo "build itself.  Run './$CMDNAME build' instead." 1>&2
                      exit 1 ;;
  build-magic )  echo "'build-magic' has been removed: 'build' now does the staged build" 1>&2
                 echo "itself, and does it correctly.  Run './$CMDNAME build' instead." 1>&2
                 exit 1 ;;
  init-all )  echo "'init-all' has been removed.  Run the commands you actually want:" 1>&2
              echo "  ./$CMDNAME pull       &&  ./$CMDNAME build        (core repositories)" 1>&2
              echo "  ./$CMDNAME pull-apps  &&  ./$CMDNAME build-apps   (core apps)" 1>&2
              exit 1 ;;

  # Core repositories
  pull )    pull ;;
  switch )  switch ;;
  branch )  branch ;;
  build )   build ;;

  # Core apps
  pull-apps )      pull-apps ;;
  switch-apps )    switch-apps ;;
  branch-apps )    branch-apps ;;
  build-apps )     build-apps ;;
  validate-apps )  validate-apps ;;

  # Everything else
  status )       status ;;
  push )         push ;;
  reset )        reset ;;
  run-all )      run-all ;;
  * )      echo "Invalid command $COMMAND: $ERROR_MESSAGE"
          exit 1;;
esac
