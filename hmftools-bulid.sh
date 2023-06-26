#!/bin/bash

set -e

# Verify input
if [ $# -ne 1 ]; then
  echo "Invalid arguments. Usage: hmftools-build semver-version"
  exit 1
fi

SEMVER_VERSION="$1"
# The regex pattern is <module_name>-<semver>
SEMVER_REGEX="^([a-zA-Z0-9_-]+)-([0-9]+)\.([0-9]+)\.([0-9]+)(-(alpha|beta)\.[0-9]+)?$"
if ! [[ "$SEMVER_VERSION" =~ $SEMVER_REGEX ]]; then
  echo "Invalid semver version: $SEMVER_VERSION"
  exit 1
fi

BUILD_MODULE=${BASH_REMATCH[1]}
SEMVER="${SEMVER_VERSION/${BASH_REMATCH[1]}-}"

# Check if this module exists
if ! [[ -d "$BUILD_MODULE" ]]; then
  echo "Invalid semver version: module $BUILD_MODULE does not exist"
  exit 1
fi

# Version the appropriate files
mvn -f "./${BUILD_MODULE}/pom.xml" versions:set -DnewVersion=${SEMVER}

# Step 2: compile, package, release
mvn -f "./${BUILD_MODULE}/pom.xml" deploy -B


echo "done!"
