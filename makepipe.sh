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
  * help: display a command's help screen
  * build: yaml: generate a Makefile from a YAML pipeline
  * brick: list all bricks available
  * brick name: display a brick's help
  * init: init a new project with makepipe in current folder
  * update: update makepipe submodule from repository
"
}

# do init makepipe
init() {
	# Where to take the repository for makepipe
	repository="$1"
  # Check git not present in current folder
  if [ -d .git ]
  then
    echo "You have already a git repository in this directory"
    echo "You can add makepipe as submodule by doing:"
    echo "git submodule add $repository makepipe"
    exit 0
  fi
	# init git repository for this projects
	git init
	git submodule add $repository makepipe
	$(cd makepipe && git checkout -b local)
	git add .
	git commit -m 'Initial release on local branch'
	version=$(git --git-dir=makepipe/.git tag | grep makepipe | tail -1)
	if [ $version ]
	then
		git tag $version
	fi
	cat <<END
	[info] Makepipe installed
	[info] Git Repository initialized

# Configure your yml file, ckeck makepipee/pipeline.yml as sample configuration file.
For more information, read makepipe/README.md
 * The original bricks are managed in the submodule "makepipe/bricks".
 * If you modify a brick, you should commit inside the submodule "makepipe" on
   the local branch.
 * If you add a new brick, you can add to your current git repository (inside the bricks folder)
 * or put it in the bricks folder of the submodule "makepipe" and index it in the local
branch.
END
}

#----------
upgrade0.03info() {
  cat <<END
----
Something went wrong & Do it by yourself:
# Files NOT in bricks or makepipe.module
 * You have files to commit or files to index/commit"
 * Do 'git add & git commit'
 * or git commit -am 'your message'

# Files IN bricks or makepipe.module
 * If new file in bricks/*, move them in makepipe.module/bricks/
 * cd makepipe.module
 * create a local branch in makepipe.module:
      git checkout -b local
      git commit -am 'your message'
      git chekckout master
      git pull origin master
      git rebase local
      git checkout local
  * cd ..
  * git commit -a -m "Makepipe ready for upgrade from v0.03"
  * makepipe update

# More info on submodule
 * see https://git-scm.com/book/en/v2/Git-Tools-Submodules
END
  exit 1
}

# do update makepipe
update() {
	# check if current project repository is clean
	gitstatus=$(git status -s | grep -v '?')
	if [ "$(git status -s | grep -v '?')" != "" ]
	then
    echo '<-------'
    echo $gitstatus
    echo '------->'
    upgradev0.03info
    exit 0
	fi
	# update current git module repository for this projects
  echo "[info] Doing update of makepipe"
  # update from v0.03
  if [ -d makepipe.module ]
  then
    [ -L makepipe ] && rm makepipe
    for f in bricks/*; do [ -L $f ] && rm $f; done
    rmdir bricks && mv makepipe.module makepipe
    if [ $? ]
    then
      upgradev0.03info
    fi
    git commit -a -m 'Updated Makepipe from version 0.03'
    echo "Updated from version 0.03"
  fi
  git submodule update --remote --rebase
  echo "[info] Makepipe updated"
}
############

if [ -d makepipe.module ]
then
  echo "You HAVE to upgrade to makepipe v0.04, by doing:"
  echo "makepipe update"
	makepipe="./makepipe.offline"
fi

if [ -x "$pwd/makepipe/makepipe" ]
then
	makepipe="./makepipe/makepipe"
else
	echo -e "[Warning] Makepipe not initialised\n"
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
		"list")
			usage
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

