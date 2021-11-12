# iClusion-Importer

This module queries the iClusion API for all trials available in their database along with molecular inclusion criteria and eligible tumor types.
The iClusion database holds all cancer trials which are currently recruiting in the Netherlands and can be accessed through an account provided by iClusion.
More information about iClusion can be found on their [website](https://iclusion.org). 

Some specific notes about the iClusion Importer:
 * The trials are written to a specific TSV format where every line represents a trial and specific deserialization is needed to understand the actual trials.
 * There is no filtering or curation done by the importer. The data from the TSV is an exact representation of the iClusion database.
 * The module provides logic to map every mutation to an Event Type as defined by [SERVE](../serve/README.md).
 * The module provides logic to map every trial into a set of [SERVE](../serve/README.md) actionable events

## Running the iClusion Importer

In order to run the iClusion importer one needs to get an account from iClusion. 
Assuming this has been provided, the trials can be imported by running the iClusion importer jar with the following params:
  
Argument  | Description
---|---
iclusion_endpoint | Required: The URL where the iClusion API runs that one wants to query.
iclusion_client_id | Required: Provided by iClusion.
iclusion_client_secret | Required: Provided by iClusion.
iclusion_username | Required: Provided by iClusion.
iclusion_password | Required: Provided by iClusion.
iclusion_trial_tsv | Required: The path to where the trials are going to be written to.

## Version History and Download Links
- Upcoming
  - Extend config file of events
- [1.3](https://github.com/hartwigmedical/hmftools/releases/tag/iclusion-importer-v1.3)
  - This release is a bugfix release to be able to pick up DOIDs which were dropped from the API a while ago and now added back in
- [1.2](https://github.com/hartwigmedical/hmftools/releases/tag/iclusion-importer-v1.2)
  - New datamodel for iclusion
- [1.1](https://github.com/hartwigmedical/hmftools/releases/tag/iclusion-importer-v1.1)
- [1.0](https://github.com/hartwigmedical/hmftools/releases/tag/iclusion-importer-v1.0)
  - Initial release 
