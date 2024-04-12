# ORANGE Datamodel

The ORANGE datamodel module contains all the datamodels that [ORANGE](../orange) will produce as part of the JSON output.

The goal of this module is to make it easy to parse the datamodel:

- The module has as few dependencies as possible, to make it easy for downstream tools to depend on it.
- The module should be stable and present a clean datamodel to depend on.

## Releasing / deploying orange-datamodel

The way to make a new `orange-datamodel` release and deployment is to create a tag in the proper format (`orange-datamodel-vX.Y.Z`).

This triggers a GCP-cloud-build ultimately depositing the right artifacts in the right place for downstream clients to depend on (see
also [hmftools-build.py](../hmftools-build.py))