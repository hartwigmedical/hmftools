# ID Generator

The ID generator maps all original sample IDs to anonymous HMF IDs. The benefits of this are:
 - There is no reference anymore to the original sample ID which is known to the group that contributed the sample.
 - It becomes immediately clear which samples belong to the same patient.
 
 ### ID anonymization
 
The ID anonymization is based on the patient mapping performed by [AMBER](../amber/README.md).
For every sampleId we create a hash using the anonymization password and link this to the patient derived from amber patient mapping.
The sampleId and the hash can subsequently be used to map towards the HMF sample ID.

It occasionally happens that samples are changed after ingestion into the database (eg because the entry is pulled from the HMF database).
To support these events there is the option to soft-delete an entry in amberAnonymous. This way we can retain the HMF ID for this sample 
while making clear the sample is no longer part of the current HMF database. 

The tool expects a properly populated amberPatient in the database that is connected to, and can be run as follows:
 ```bash
java -jar /path/to/id_generator_jar \
    -password ${anonymization_password} \
    -input_sample_file "/path/to/sample_hashes.csv" \
    -output_sample_file "/path/to/new_sample_hashes.csv" \
    -db_user ${db_user} -db_pass ${db_pass} -db_url ${db_url}
 ```

Do note:
 - Parameters "input_sample_file" and "output_sample_file" can point to the same file in which case the sample hashes will be overwritten.
 - A "new_password" parameter can optionally be provided to reset the anonymization password to a new value.  

The following checks are done by the algorithm:
 - Every sample that has a mapping in amberAnonymous and is not (soft-)deleted is expected to exist in amberPatient
 - Every sample that exists in amberPatient should not be (soft-)deleted in amberAnonymous
 - The amberAnonymous table should be completely in sync with the samples implied by the input file hashes 

Assuming the above checks succeed, the ID generator will do the following:
 - Generate a new hash for every sample in amberPatient that doesn't have a hash yet. 
 - Write the new hashes to the "out" parameter provided to the command line
 - Repopulate the table "amberAnonymous" with the new sampleId -> hmfSampleId mappings.


 
   