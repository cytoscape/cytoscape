#!/bin/sh
###############################################################################
#
# SCRIPT: release.sh
#
# DESCRIPTION: 
#    Cytoscape Release Script:
#
#    Create a new release from hotfix branch.
#    This is an interactive script and user should answer
#
#  Note: This script updates local repository only.
#        For the actual release, you need to push the changes to the upstream.
#
# By Keiichiro Ono (kono at ucsd edu)
#
###############################################################################

##########################
#
# Reset all local changes.
#  This deletes all local changes!
#
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


CMDNAME=$(basename $0)
ERROR_MESSAGE="Usage: $CMDNAME [-c branch_name] release_version_number branch_name"

##### Target Repositories #####
REPOSITORIES=(api impl gui-distribution app-developer)


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


#######################################
# Reset everything                    #
#######################################
if [ "$FLG_C" ]; then
	echo " - Clear option is ON."
	# Cleanup and exit
	if [ -z $branch ]; then
		echo "Branch name is required to reset."
		exit 1
	fi

	for project in "${REPOSITORIES[@]}"; do
		echo "\n - Resetting local changes: $project"
		cd ../$project
		resetAll
	done

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
#
# Checkout, pull, and build the latest version.
#
###############################################
function prepare {

	# Display available branches.
	echo "\n\n - Available branches:\n\n"
	git branch -avv || { echo Failed to get list of branches.; exit 1; }

	# Switch to the target branch.
	echo "\n\n - Checking out branch: $branch"
	git checkout $branch

	# Get chenges from remote. 
	echo "\n\n - Pulling changes from github remote repository:\n\n$(git remote -v)\n\n"
	git pull || { echo Could not pull the changes from remote; exit 1; }

}


###########################################################
#
# Build modified (version updated) version for testing.
#
###########################################################
function buildNewVersion {

	echo "\n\nBuilding new version ($version) of Cytoscape for testing...\n\n"
	mvn clean install || { echo Build Failed; exit 1; }
	echo "\n\nBuild Success!\n\n"
	# Cleanup
	mvn clean
}


#############################################################
#
# Bunp the version numbers in pom files.
#
#############################################################
function updateVersionNumbers {

	# Use Maven versions plugin
	echo "\n\n - Updating version numbers...\n\n"
	mvn versions:set -DnewVersion=$version || { echo Could not update project version.; exit 1; }
	mvn versions:commit || { echo Could not commit change to the version number.; exit 1; }
}


#####################################################################
#
# Modify properties in the pom.
#
#####################################################################
function updateProperties {
	echo "\n\n - Updating version numbers in properties..."
	sed -E -e "s/-SNAPSHOT//g" pom.xml > pom.updated.xml
	mv pom.updated.xml pom.xml
}


#####################################################################
#
# Modify features.xml: Simply remove all SNAPSHOTS.
#
#####################################################################
function updateFeatures {
	cd features/src/main/resources
	sed -E -e "s/-SNAPSHOT//g" features.xml > features.updated.xml
	
	echo "\n\n - Here is the changes:\n\n"
	mv features.updated.xml features.xml
	cd -
}


#####################################################################
#
# Ask user to proceed or not.
#
#####################################################################
function confirm {
	echo "\n\n - Do you want to continue?  Type [Y] to proceed: "
	read answer

	if [[ $answer != 'Y' || -z $answer ]]; then
		echo "Abort\n"
		exit 0
	fi
}


#####################################################################
#
# Merge the branch to master and tag it.
#
#####################################################################
function mergeAndTag {

	#git diff --stat
	#confirm

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


###############################################################################
# Main Workflow                                                               #
###############################################################################
echo "\n Cytoscape Release Builder:\n"
cd ../
echo " - Cytoscape project top-level directory: $(pwd)"


#######################################
# Release API bundles
#######################################
cd api
echo "\n\n - Releasing API bundles...\n"

prepare
updateVersionNumbers
buildNewVersion

##### Create JavaDoc #####
cd swing-app-api
mvn javadoc:javadoc install || { echo Could not create JavaDoc.; exit 1; }
cd ..

mergeAndTag

echo "\n - API bundles released."

#######################################
# Release IMPL bundles
#######################################
cd ../impl
echo "\n\n - Releasing impl bundles...\n"

prepare

# Make sure everything works fine.
mvn clean install || { echo Could not build SNAPSHOT version.  Fix this problem first before release.; exit 1; }

updateVersionNumbers
updateProperties

##### Update Help from Wiki #####
cd help-impl
ant update
mvn clean install

##### Need to commit new manual #####
git add src/docbkx/manual.xml
git commit -m "Manual XML document updated for $version release."

cd ..

buildNewVersion
mergeAndTag

echo "\n - IMPL bundles released."

#######################################
# Release GUI-Distribution
#######################################
cd ../gui-distribution/
echo "\n\n - Releasing GUI-Distribution...\n"

prepare
updateVersionNumbers
updateProperties
updateFeatures

##### Create new splash screen #####
cd splash-launcher
mvn clean install
cp target/classes/images/CytoscapeSplashScreen.png ../packaging/src/main/images/

##### Commit changes to local repository
git add ../packaging/src/main/images/CytoscapeSplashScreen.png
git commit -m "Splash Screen image file had been updated for $version release."
cd ..

buildNewVersion

##### Build installers #####
cd packaging
mvn clean install || { echo Could not create installers.; exit 1; }
cd ..

mergeAndTag

echo "\n - GUI-Distribution released."

#######################################
# Release App-developer bundle
#######################################
cd ../app-developer || { echo Could not find app developer directory.; exit 1; }
echo "\n\n - Releasing app developer...\n"

prepare

# Special case:  Need to use specific profile to include packaging
mvn -P release versions:set -DnewVersion=$version || { echo Could not update project version.; exit 1; }
mvn -P release versions:commit || { echo Could not commit change to the version number.; exit 1; }

updateProperties

buildNewVersion

##### ZIP Javadoc #####

mergeAndTag

cd ..

echo "\n\n=========== Cytoscape $version is ready to release. ============="
echo "=========== Don't forget to deploy bundles to Nexus and push changes to upstream! ============="
