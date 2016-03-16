#!/bin/bash
# Init & Create the makepipe in local directory

# Where is located the $MAKEPIPE_REPOSITORY if not defined from env
if [ $GITREPOSITORY ]
then
	if [ -d $GITREPOSITORY/makepipe ]
	then
		MAKEPIPE_REPOSITORY=$GITREPOSITORY/makepipe
	fi
elif [ -d /data/share/repository/makepipe ]
then
	# define a local
	MAKEPIPE_REPOSITORY=/data/share/repository/makepipe
fi
if [ ! -d $MAKEPIPE_REPOSITORY ]
then
	echo "There no repository at $MAKEPIPE_REPOSITORY"
	echo "Check the environmental variable GITREPPOSITORY"
	echo "or do makepipe init makepipe_repository"
	exit 1
fi

# current patch
pwd=$(pwd)

###########
# Functions
#----------
# show help for all command
function usage() {
echo "
	Makepipe commands : 
	 * init: init a new project with makepipe
	 * update: update makepipe submodule from repository
	 * build yml: build Makefile from yml file
	 * brick: list all bricks available
	 * brick name: info about a brick
	 * list: list commands available
"
}

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
	echo <<END
	Makepipe installed
	Git Repository initialized
	Configure your yml file, ckeck makepipe.module/pipeline.yml for a sample configuration
	For more information, read makepipe.module/README.md
	#######
	The original bricks are followed by the submodule "makepipe.module" (via symlink).
	If you modify a brick, you should commit inside the submodule "makepipe.module" on
    the local branch.
	If you add a new brick, you can add to your current git repository (inside the bricks folder) 
	or put it in the bricks folder of the submodule "makepipe.module" and index it in the local
	branch (don't forget the symlink in your bricks folder).  	
END
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

if [ -x "$pwd/makepipe" ]
then
	makepipe="./makepipe"
else
	echo -e "\tWarning: makepipe not initialised\n"
fi


# Manage help commands
if [ "$1" == "help" ]
then
	case "$2" in
			"init")
				echo -e "Create git repository, add makepipe as submodule.\nmakepipe init [$MAKEPIPE_REPOSITORY]"
				;;
			"update")
				echo -e "Update local branch of submodule makepipe. Do \nmakepipe update [$MAKEPIPE_REPOSITORY]"
				;;
			"list")
				usage
				;;
			"build" | "brick")
				if [ $makepipe ]
				then
					$makepipe help $2
				else
					echo "To get the help on "$2", you have to init makepipe"
				fi
				;;
			*)	
				usage
				;;
	esac
	exit 0
fi

# Test if makepipe is already here
if [ $makepipe ]
then
	# if we don't have an argument, default is to send to local makepipe
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
		if [ $1 ]
		then
			MAKEPIPE_REPOSITORY=$1
		fi
		init $MAKEPIPE_REPOSITORY
	else
		usage
	fi

fi

