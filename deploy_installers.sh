#!/bin/sh
################################################################################
#
# Build & Upload Installers 
#   Currently this script works only on Kei's Mac.
#
# Options:
#
#  -a: Upload API docs
#  -u: Update API and IMPL bundles
#
################################################################################

CMDNAME=$(basename $0)
ERROR_MESSAGE="Usage: $CMDNAME [-a] [-u] branch_name"

while getopts 'au' OPT
do
	case $OPT in
		a)	FLG_A=1 ;;
		u)	FLG_U=1 ;;
		?)	echo $ERROR_MESSAGE 1>&2
				exit 1 ;;
	esac
done

if [ "$FLG_A" ]; then
	echo " - API Upload option is on: JavaDoc will be uploaded later."
fi

if [ "$FLG_U" ]; then
	echo " - Update option is on: Build API and IMPL first."
fi

shift $(($OPTIND - 1))

branch=$*

if [[ -z $branch ]]; then
	echo "Branch name is required. $ERROR_MESSAGE" 1>&2
	exit 1
fi

echo "Target branch is: $branch"

function build {
	git checkout $branch || { echo Abort: Could not switch to $branch; exit 1; }
	git pull || { echo Abort: Failed to pull changes.; exit 1; }

	mvn clean install || { echo Abort: Build Failed; exit 1; }
}

###################################
# Build API and Impl if necessary
###################################

if [ "$FLG_U" ]; then
	echo "Building API bundles..."
	cd ../api
	build
	
	echo "Building IMPL bundles..."
	cd ../impl
	build
fi


########################
# Build distribution
########################
echo "Building GUI Distribution..."
cd ../gui-distribution
build

cd ./packaging
mvn clean install || { echo Abort: Installer Build Failed; exit 1; }
./sign-dmg.sh 'Developer ID Application' || { echo Abort: Could not sign DMG for Mac.; exit 1; }

cd ./target/install4j

echo "Uploading new installers..."
deployDir="build-$(date "+%Y-%m-%d-%H-%M-%S-%Z")"

server="grenache:~/public_html/data/cy3latest/$deployDir/"

# Create directory
ssh grenache mkdir /cellar/users/kono/public_html/data/cy3latest/$deployDir

scp ./signed/* $server || { echo Abort: Could not upload signed file.; exit 1; }
scp ./*.exe $server || { echo Abort: Could not upload Windows installers.; exit 1; }
scp ./*.zip $server || { echo Abort: Could not upload zipped file.; exit 1; }
scp ./*.sh $server || { echo Abort: Could not upload UNIX installer.; exit 1; }


if [ "$FLG_A" ]; then
	# API Bundles
	cd ../../../../app-developer/
	build
	cd ./target
	scp -r ./API $server || { echo Abort: Could not upload JavaDoc files.; exit 1; }
fi

deployURL="http://chianti.ucsd.edu/~kono/data/cy3latest/$deployDir"
echo "Finished.  Visit $deployURL"
open $deployURL
