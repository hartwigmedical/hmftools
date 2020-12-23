# iClusion-Importer

This module queries the iClusion API for all trials available in their database along with molecular inclusion criteria and allowed tumor types.
The iClusion database holds all cancer trials which are currently recruiting in the Netherlands and can be accessed through an account provided by iClusion.
More information about iClusion can be found on their [website](https://iclusion.org). 

Some specific notes about the iClusion Importer:
 * The trials are written to a fairly specific TSV format where every line represents a trial and some specific deserialization is needed to understand the actual trials.
 * There is no filtering or curation done by the importer. The data from the TSV is an exact representation of the iClusion database.
 * The module provides logic to map every mutation to an Event Type as defined by [SERVE](../serve/README.md).
 * The module provides mapping logic to map every trial into a set of [SERVE](../serve/README.md) actionable events

## Running the iClusion Importer

In order to run the iClusion importer one needs to get an account from iClusion. 
Assuming this has been provided, the trials can be imported by running the iClusion importer jar with the following params:
  
Argument  | Description
---|---
vicc_json  | Required: Path to the VICC json file containing the VICC data.
iclusion_endpoint | Required: The URL where the iClusion API runs that we are going to query
iclusion_client_id | Required: Provided by iClusion
iclusion_client_secret | Required: Provided by iClusion
iclusion_username | Required: Provided by iClusion
iclusion_password | Required: Provided by iClusion
iclusion_trial_tsv | Required: The path to where the trials are going to be written to.
