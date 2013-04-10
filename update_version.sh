#!/bin/sh

#####################################
#
# Cytoscape Version Updater
#
#####################################

CMDNAME=$(basename $0)

##########################
# Reset all local changes.
##########################
function resetAll {

	git checkout master
	git reset --hard $(git branch -av | grep "remotes/origin/master" | awk '{ print $2 }')
	git checkout develop
	git reset --hard $(git branch -av | grep "remotes/origin/develop" | awk '{ print $2 }')
	git checkout $branch
	git reset --hard $(git branch -av | grep "remotes/origin/$branch" | awk '{ print $2 }')
	git branch -avv

}

ERROR_MESSAGE="Usage: $CMDNAME [-c branch_name] version_number branch_name"

while getopts 'c:' OPT
do
	case $OPT in
		c)	FLG_C=1
				branch="$OPTARG"
				;;
		?)	echo $ERROR_MESSAGE 1>&2
				exit 1 ;;
	esac
done

if [ "$FLG_C" ]; then
	echo " - Clear option is ON."
	# Cleanup and exit
	if [ -z $branch ]; then
		echo "Branch name is required to reset."
		exit 1
	fi

	cd ../api
	resetAll
	cd ../impl
	resetAll
	cd ../gui-distribution
	resetAll

	echo "\n\n - Everything had been reset to HEAD of remote branches."
	exit 0
fi

shift $(($OPTIND - 1))
branch=$2
version=$1

echo "Version = $version"
echo "Branch = $branch"

if [[ -z $branch ]]; then
	echo "Branch name is required. $ERROR_MESSAGE" 1>&2
	exit 1
fi

if [[ -z $version ]]; then
	echo "Version number is required. $ERROR_MESSAGE" 1>&2
	exit 1
fi

echo " --> Target version is: $version"
echo " --> Target branch is: $branch"

###############################################
# Checkout, pull, and build the latest version.
###############################################
function prepare {

	echo "\n\n - Available branches:\n\n"
	git branch -avv || { echo Failed to get list of branches.; exit 1; }

	echo "\n\n - Checking out branch: $branch"
	git checkout $branch

	echo "\n\n - Pulling changes from github remote repository:\n\n$(git remote -v)\n\n"
	git pull || { echo Could not pull the changes from remote; exit 1; }

	#echo "\n\nBuilding Cytoscape for testing...\n\n"
	#mvn clean install || { echo Build Failed; exit 1; }
}

###########################################################
# Build modified (version updated) version for testing.
###########################################################
function buildNewVersion {
	echo "\n\nBuilding new version ($version) of Cytoscape for testing...\n\n"
	mvn clean install || { echo Build Failed; exit 1; }
	echo "\n\nBuild Success!\n\n"
	# Cleanup
	mvn clean
}

function updateVersionNumbers {
	echo "\n\n - Updating version numbers...\n\n"
	
	mvn versions:set -DnewVersion=$version || { echo Could not update project version.; exit 1; }
	mvn versions:commit || { echo Could not commit change to the version number.; exit 1; }
}

function updateProperties {
	echo "\n\n - Updating version numbers in properties..."
	sed -E -e "s/-SNAPSHOT//g" pom.xml > pom.updated.xml

	echo "\n\n - Here is the changes:\n\n"
	diff pom.xml pom.updated.xml

	echo "\n\n Do you want to update pom.xml files?"
	confirm
	mv pom.updated.xml pom.xml
}

function updateFeatures {
	cd features/src/main/resources
	sed -E -e "s/-SNAPSHOT//g" features.xml > features.updated.xml
	
	echo "\n\n - Here is the changes:\n\n"
	diff features.xml features.updated.xml

	printf "\n\n Do you want to update features.xml files? [Y to proceed]: "
	confirm

	mv features.updated.xml features.xml
	cd -
}

function confirm {
	echo "\n\n - Do you want to continue?  Type [Y] to proceed: "
	read answer

	if [[ $answer != 'Y' || -z $answer ]]; then
		echo "Abort\n"
		exit 0
	fi
}

function mergeAndTag {

	git diff --stat
	confirm

	# First, commit the pom files with new version numbers
	git commit -am "Version numbers are updated to $version.  This will be tagged to $version release."

	git checkout master
	git merge --no-ff $branch -m "$branch merged.  This is $version release."

	# Delete tag is exists
	git tag -d $version
	# Create Tag
	git tag -a $version -m "$version Release."

	#git branch -d $branch

	git branch -avv
	git tag

}

################### Main Workflow #############################################


echo "\n\n================================================================"
echo "\n Cytoscape Release Builder:\n"

cd ../

echo " - Cytoscape project top-level directory: $(pwd)"
confirm

# API bundles
cd api
echo "\n\n - Processing API bundles...\n"
prepare
updateVersionNumbers
buildNewVersion
# Finish the branch
mergeAndTag

echo "\n - Finished API bundles."
confirm


# IMPL bundles
cd ../impl
echo "\n\n - Processing impl bundles...\n"
prepare
updateVersionNumbers
updateProperties
cd help-impl
mvn clean install
cd ../model-impl
mvn clean install
cd ..

buildNewVersion

mergeAndTag

echo "\n - Finished IMPL bundles."
confirm

# GUI distribution
cd ../gui-distribution/
echo "\n\n - Processing GUI-Distribution bundles...\n"
prepare
updateVersionNumbers
updateProperties
updateFeatures
buildNewVersion
mergeAndTag

echo "\n - Finished GUI-Distribution bundles."
echo "\n\n=========== Done! ============="
