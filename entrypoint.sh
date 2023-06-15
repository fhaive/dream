#!/bin/sh

set -e

# use specified user name or use `default` if not specified
USER="${USER:-dream}"

# use the specified UID for the user
UID="${UID:-1000}"

# use the specified GID for the user
GID="${GID:-${UID}}"


groupadd -g $GID -o $USER
useradd -m -u $UID -g $GID -o -s /bin/bash $USER

# exec and run the actual process specified in the CMD of the Dockerfile (which gets passed as ${*})
exec gosu "${USER}" "${@}"