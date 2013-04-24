#!/bin/sh
###############################################################################
#
# SCRIPT: cy
#
# DESCRIPTION: Cytoscape 3 Repository management utility.
# 	This script is only for core developers.
#
# By Keiichiro Ono (kono at ucsd edu)
#
###############################################################################

function reset {
	for project in "${REPOSITORIES[@]}"; do
		echo "\n - Resetting local changes: $project"
		cd ../$project
		resetAll
	done

}

function resetAll {
	# TODO: 
	
	git checkout master
	git reset --hard $(git BRANCH -av | grep "remotes/origin/master" | awk '{ print $2 }')
	git checkout develop
	git reset --hard $(git BRANCH -av | grep "remotes/origin/develop" | awk '{ print $2 }')
	git checkout $BRANCH
	git reset --hard $(git BRANCH -av | grep "remotes/origin/$BRANCH" | awk '{ print $2 }')
	git BRANCH -avv
}

function confirm {
	echo "Do you really want to continue?  Type [yes] to proceed: "
	read answer

	if [[ $answer != 'yes' || -z $answer ]]; then
		echo "Abort\n"
		exit 0
	fi
}

function init {
	echo "Please enter target directory [default is parent of $START_DIR]: "
	read TARGET_DIR
	if [[ -z "$TARGET_DIR" ]]; then
		cd ..
		TARGET_DIR=$(pwd)
	elif ! [ -e "$TARGET_DIR" ]; then
		echo "No such dir: $TARGET_DIR"
		exit 1
	else
		cd $TARGET_DIR || { echo Could not find target directory: $TARGET_DIR; exit 1; }
	fi

	echo "Target directory: $TARGET_DIR/cytoscape3"
	mkdir cytoscape3
	cd cytoscape3
	
	CY3_DIR=$(pwd)

	for REPO in "${REPOSITORIES[@]}"; do
		REPO_URL="$BASE_URL$REPO.git"
		echo "Cloning: $REPO (URI = $REPO_URL)"
		git clone $REPO_URL $REPO
		cd $REPO
		git checkout master
		git flow init -d
		git checkout develop
		cd ..
	done

	# Copy required files
	cd $START_DIR
	cp pom.xml $CY3_DIR/

	echo "\n\n - Finished: here is the current status:\n"

	cd $CY3_DIR
	for REPO in "${REPOSITORIES[@]}"; do
		cd $REPO
		echo "$REPO:\n $(git branch -vv)\n"
		cd ..
	done
}


# Command Name
CMDNAME=$(basename $0)

# Error Message
ERROR_MESSAGE="Usage: $CMDNAME [-h] command"

# Help
HELP='Cytoscape repository management command'

# Git base URL
BASE_URL='git@github.com:cytoscape/cytoscape-'

# Cytoscape Repositories
REPOSITORIES=(parent api impl support headless-distribution gui-distribution app-developer samples)

while getopts 'h' OPT
do
	case $OPT in
		h)	FLG_H=1
				echo "$HELP: $ERROR_MESSAGE"
				exit 0
				;;
		?)	echo $ERROR_MESSAGE 1>&2
				exit 1 ;;
	esac
done

shift $(($OPTIND - 1))

COMMAND=$1
echo "Command = $COMMAND"

if [[ -z $COMMAND ]]; then
	echo "COMMAND is required. $ERROR_MESSAGE" 1>&2
	exit 1
fi

echo " --> Command is: $COMMAND"


################### Main Workflow #############################################

START_DIR=$(pwd)

case $COMMAND in
	init )	init ;;
	* )			echo "Invalid command $COMMAND: $ERROR_MESSAGE"
					exit 1;;
esac


