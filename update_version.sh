#!/bin/sh

#####################################
# Cytoscape Version Updater 
#
#####################################

###############################################
# Checkout, pull, and build the latest version.
###############################################
function prepare {

	echo "\n\n - Available branches:\n\n"
	git branch -av || { echo Failed to get list of branches.; exit 1; }

	echo "\n\n - Checking out branch: $branch"
	git checkout $branch

	echo "\n\n - Pulling changes from github remote repository:\n\n$(git remote -v)\n\n"
	git pull || { echo Could not pull the changes from remote; exit 1; }

	echo "\n\nBuilding Cytoscape for testing...\n\n"
	mvn clean install || { echo Build Failed; exit 1; }
}

###########################################################
# Build modified (version updated) version for testing.
###########################################################
function buildNewVersion {
	echo "\n\nBuilding new version ($version) of Cytoscape for testing...\n\n"
	mvn clean install || { echo Build Failed; exit 1; }
	echo "\n\nBuild Success!\n\n"
}

function updateVersionNumbers {
	echo "\n\n - Updating version numbers...\n\n"
	
	mvn versions:set -DnewVersion=$version || { echo Could not update project version.; exit 1; }
	mvn versions:commit || { echo Could not commit change to the version number.; exit 1; }
}

function updateProperties {
	echo "\n\n - Updating version numbers in properties..."
	sed -E -e "s/-SNAPSHOT//g" pom.xml > pom.updated.xml
#		| sed -e 's/>/##/g' | sed -e 's/<//g' | awk -F"##" '{ print $1 " == " $2 }'

	echo "\n\n - Here is the changes:\n\n"
	diff pom.xml pom.updated.xml

	printf "\n\n Do you want to update pom.xml files? [Y to proceed]: "
	read proceed
	if [[ $answer != 'Y' || -z $answer ]]; then
		echo "Abort\n"
		exit 0
	fi
	
	mv pom.updated.xml pom.xml
}


################### Main Workflow #############################################


### Error check
if (($# != 2)); then
	printf "%b" "Usage: update_version.sh new_version branch_name\n" >&2
	exit 1
fi

# Target version number for the release.
version=$1
branch=$2

echo "\n\n================================================================"
echo "\n Cytoscape Release Builder:\n"

cd ../

echo " - Cytoscape project top-level directory: $(pwd)"
echo "This script builds new Cytoscape release version = $version, from branch $branch.  Type [Y] to proceed: "

read answer

if [[ $answer != 'Y' || -z $answer ]]; then
	echo "Abort\n"
	exit 0
fi


# API bundles
#cd api
#prepare
#updateVersionNumbers
#buildNewVersion

cd impl
echo "\n\n - Processing impl bundles...\n"
updateProperties


# IMPL bundles
#prepare
#updateVersionNumbers
#buildNewVersion
#cd ../gui-distribution/


echo "=========== Done! ============="
