package com.hartwig.hmftools.idgenerator

// CLI options
const val CREATE_SINGLE_HASH_MODE = "create_single_hash"
const val CREATE_IDS_MODE = "create_ids"
const val UPDATE_IDS_MODE = "update_ids"
const val UPDATE_IDS_WITH_AMBER_MODE = "update_ids_amber"
const val ANONYMIZE_IDS_MODE = "anonymize_ids"

const val PASSWORD = "password"
const val NEW_PASSWORD = "new_password"
const val SAMPLE_ID = "sample_id"
const val SAMPLE_IDS_FILE = "sample_ids_file"
const val PATIENT_MAPPING_FILE = "patient_mapping_file"
const val AMBER_PATIENT_FILE = "amber_patient_file"
const val OUTPUT_FILE = "out"
const val SAMPLE_MAPPING_OUTPUT_FILE = "mapping_out"
const val ANONYMIZE_OUT = "anonymize_out"

// Resource locations
const val SAMPLE_HASHES_CSV = "/sample_hashes.csv"

object Version {
    override fun toString(): String = this::class.java.getPackage().implementationVersion
}
