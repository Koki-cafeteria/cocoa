if [ -z "${ROOTDIR}" ]; then
    echo 'ERROR ROOTDIR not defined'
    exit 1
fi

git pull origin main
source $ROOTDIR/projects/.git_pull.sh