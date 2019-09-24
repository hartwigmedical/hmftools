# ID Generator

The ID generator translates all original sample IDs to anonymous HMF IDs. The main uses of this tool are as follows:
 - Generate hashes for a selection of sampleIDs such that a subsequent release of ID generator can anonymize these samples
 - Create anonymous ID for a sample using the internal hashes of the ID generator
 
 ### ID anonymization
 
To generate a set of anonymous IDs for a set of sample IDs, you need:
 - to use the same password that was used to create the initial hashes for these samples.
 - to use a release of id generator that holds the hashes for every sample.
 
 ```bash
java -jar id-generator.jar \
    -anonymize_ids \ 
    -password ${secret_password_that_was_used_for_hash_creation} \
    -sample_ids_file /path/to/sample_ids_file \
    -patient_mapping_file /path/to/patient_mapping_file
 ```

The patient mapping file can map 2 sample IDs onto each other if you know these samples belong to the same patient 
(but participated in multiple studies for example). 
   