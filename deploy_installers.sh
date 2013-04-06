#!/bin/sh
###############################################################################
#
# Build & Upload Installers 
#   Currently this script works only on Kei's Mac.
#
###############################################################################
if (($# != 1)); then
	printf "%b" "Usage: deploy_installers.sh branch_name\n" >&2
	exit 1
fi

branch=$1

function checkoutAndBuild {
	# 1. Clone remote if not available
	dist=$(ls -alh | grep $targetProject)
	echo "Directory is: $dist"

	if [[ -z $dist ]]; then
		echo "$targetProject does not exist.  Cloning..."
		git clone https://github.com/cytoscape/$targetProject.git
	fi

	# 2. Build it
	cd $targetProject 

	# Switch branch
	git checkout $branch || { echo Abort: Could not switch to $branch; exit 1; }
	git pull || { echo Abort: Failed to pull changes.; exit 1; }

	mvn clean install || { echo Abort: Build Failed; exit 1; }
}

targetProject="cytoscape-gui-distribution"
checkoutAndBuild

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


# API Bundles
targetProject="cytoscape-app-developer"
cd ../../../../
checkoutAndBuild
cd targetProject
ls
cd ./target
scp -r ./API $server || { echo Abort: Could not upload JavaDoc files.; exit 1; }

deployURL="http://chianti.ucsd.edu/~kono/data/cy3latest/$deployDir"

echo "Finished.  Visit $deployURL"
open $deployURL
