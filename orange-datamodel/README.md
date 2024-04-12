# ORANGE Datamodel

The ORANGE datamodel module contains all the datamodels that [ORANGE](../orange) will produce as part of the JSON output.

The goal of this module is to make it easy to parse the datamodel:

- The module has as few dependencies as possible, to make it easy for downstream tools to depend on it.
- The module should be stable and present a clean datamodel to depend on.

## Releasing / deploying orange-datamodel

First step is to make sure the `orange-datamodel.version` property in the main hmftools POM is set correctly. This is the version under
which the `orange-datamodel` artifact will be released.

Assuming correct version, the `orange-datamodel` can be released and deployed via
running ```mvn clean deploy -B -Pshade -pl orange-datamodel -am``` locally. This will deploy the `orange-datamodel` artifact in GCP artifact
registry under the version which is configured in the hmftools
main pom