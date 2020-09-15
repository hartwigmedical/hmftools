# ID Generator

The ID generator maps all original sample IDs to anonymous HMF IDs. The benefits of this are:
 - There is no reference anymore to the original sample ID which is known to the group that contributed the sample.
 - It becomes immediately clear which samples belong to the same patient.
 
 ### ID anonymization
 
The ID anonymization is based on the patient mapping performed by [AMBER](../amber/README.md).
For every sampleId we create a hash using the anonymization password and link this to the patient derived from amber patient mapping.
The sampleId and the hash can subsequently be used to map towards the HMF sample ID.

The tool expects a properly populated amberPatient in the database that is connected to, and can be run as follows:
 ```bash
java -jar /path/to/id_generator_jar \
    -password ${anonymization_password} \
    -in "/path/to/sample_hashes.csv" \
    -out "/path/to/new_sample_hashes.csv" \
    -db_user ${db_user} -db_pass ${db_pass} -db_url ${db_url}
 ```

Do note:
 - Parameters "in" and "out" can point to the same file in which case the sample hashes will be overwritten.
 - A "new_password" parameter can optionally be provided to reset the anonymization password to a new value.  

The following checks are done by the algorithm using the contents of the existing amberAnonymous mappings:
 - If there are no current mappings, an error is thrown since we expect mappings to exist.
 - If any of the existing mappings don't exist anymore after new round of anonymization an error is thrown since we expect all existing mappings to remain unchanged.

Assuming the above checks succeed, the ID generator will do the following:
 - Write the new hashes to the "out" parameter provided to the command line
 - Repopulate the table "amberAnonymous" with the new sampleId -> hmfSampleId mappings.


 
   