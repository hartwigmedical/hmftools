TAG_REGEX="^([a-z-]+)-v?([0-9]+\.[0-9]+\.[0-9]+(-(alpha|beta)\.[0-9]+)?)$"

if ! [[ $TAG =~ $TAG_REGEX ]]; then
  echo "Invalid tag: $TAG, expecting tag to pass the following regex: $TAG_REGEX"
  exit 1
fi

MODULE=${BASH_REMATCH[1]}
VERSION=${BASH_REMATCH[2]}
JAR_PATH="${MODULE}/target/${TAG}-jar-with-dependencies.jar"

## TODO: authenticate with gh

gh release create $TAG $JAR_PATH --draft --title "$MODULE $VERSION" --notes "This is an auto generated release."