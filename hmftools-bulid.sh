#!/bin/bash

set -e

if [ $# -ne 1 ]; then
  echo "Invalid arguments. Usage: hmftools-build semver-version"
  exit 1
fi

SEMVER_VERSION="$1"

SEMVER_REGEX="^([a-zA-Z0-9_]+)-([0-9]+)\.([0-9]+)\.([0-9]+)(-(alpha|beta)\.[0-9]+)?$"
if ! [[ "$SEMVER_VERSION" =~ $SEMVER_REGEX ]]; then
  echo "Invalid semver version: $SEMVER_VERSION"
  exit 1
fi

BUILD_MODULE=${BASH_REMATCH[1]}
SEMVER="${SEMVER_VERSION/${BASH_REMATCH[1]}-}"

echo $BUILD_MODULE
echo $SEMVER

mvn versions:set -DnewVersion=${SEMVER}