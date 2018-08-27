package com.hartwig.hmftools.knowledgebaseimporter.diseaseOntology

// MIVO: Use string instead of numeric value for DOIDs since DOIDs can have leading zeros (e.g. 0080182)
data class Doid(val value: String)
