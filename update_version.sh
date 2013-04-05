#!/bin/sh

#####################################
# Cytoscape Version Updater 
#
#####################################

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
	#grep -n -G '^\t*<cytoscape\.' < pom.xml | sed -E -e "s/-SNAPSHOT//g" | sed -e 's/>/##/g' | sed -e 's/<//g' | awk -F"##" '{ print $1 " == " $2 }'

	#sed -e "s/<cytoscape.api.version>$version-SNAPSHOT<\/cytoscape.api.version>/<cytoscape.api.version>$version<\/cytoscape.api.version>/g" pom.xml | 
	#sed -e "s/<cytoscape.impl.version>$version-SNAPSHOT<\/cytoscape.impl.version>/<cytoscape.impl.version>$version<\/cytoscape.impl.version>/g" | 
	#sed -e "s/<cytoscape.distribution.version>$version-SNAPSHOT<\/cytoscape.distribution.version>/<cytoscape.distribution.version>$version<\/cytoscape.distribution.version>/g" | tee pom.updated.xml

	#mv pom.updated.xml pom.xml
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
