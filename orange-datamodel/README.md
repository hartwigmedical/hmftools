# ORANGE Datamodel

The ORANGE datamodel module contains all the datamodels that [ORANGE](../orange) will produce as part of the JSON output.

The goal of this module is to make it easy to parse the datamodel:
 - The module has as few dependencies as possible, to make it easy for downstream tools to depend on it.
 - The module should be stable and present a clean datamodel to depend on.