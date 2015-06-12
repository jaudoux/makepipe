#!/bin/bash
# Init & Create the makepipe in local directory

# Check we have a git bare repository at $MAKEPIPE_REPOSITORY from env
if [ -z $MAKEPIPE_REPOSITORY ]
then
	MAKEPIPE_REPOSITORY=/data/share/repository/makepipe
	# define a local
fi
# current patch
pwd=$(pwd)

###########
# Functions
#----------
# do init makepipe
init() {
	# Where to take the repository for makepipe
	repository="$1"
	# init git repository for this projects
	git init
	git submodule add $repository makepipe.module
	mkdir bricks
	cd bricks
	for f in ../makepipe.module/bricks/*; do ln -s $f;done
	cd ..
	ln -s makepipe.module/makepipe
	git add .
	git commit -m 'Initial release'
	version=$(git --git-dir=makepipe.module/.git tag | grep makepipe | tail -1)
	if [ $version ]
	then
		git tag $version
	fi
	echo "Makepipe installed"
	echo "Git Repository initialized"
	echo "Configure your yml file, ckeck makepipe.module/pipeline.yml for a sample configuration"
	echo "For more information, read makepipe.module/README.md"
}

#----------
# do update makepipe
update() {
	# check if current project repository is clean
	if ! git status -s
	then
		echo "You have files to commit or files to index"
 		echo "Do 'git add & git commit'"
		echo "or git ci -a -m 'your message'"
		exit 1
	fi
	# update current git module repository for this projects
	echo "You should do:"
	echo "git submodule update"
	#version=$(git $makepiperep tag | grep makepipe | tail -1)
	echo "git commit -a -m \"Update to version $version\""
	echo "git tag $version"
	echo "Makepipe updated"
	echo "Do make to update project"
}
############

# Test if makepipe is already here
if [ -x "$pwd/makepipe" ]
then
	makepipe="./makepipe"
	# if we don't have an argument, default is 'build'
	case $1 in
		"init")
			echo "You have already initialized this repository,"
		   	echo -e "you should do \nmakepipe update"
			;;
		"update") 
			# we do update
			if [ "$2" ]
			then
				MAKEPIPE_REPOSITORY=$2
			fi
			update $MAKEPIPE_REPOSITORY
			;;
		"help")
			case "$2" in
				"init")
					echo -e "Do \nmakepipe init [makepipe repository]"
					;;
				"update")
					echo -e "Do \nmakepipe update"
					;;
				*)
					$makepipe help $2
					;;
			esac
			;;
		*)
			$makepipe $*
			#makepipe build *.yml > Makefile
			;;
	esac
else
	# we have nothing yet, do init()
	if [ "$1" == "init" ]
	then
		shift
	fi
	if [ $1 ]
	then
		MAKEPIPE_REPOSITORY=$1
	fi
	init $MAKEPIPE_REPOSITORY

fi

