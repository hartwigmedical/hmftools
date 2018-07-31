package com.hartwig.hmftools.knowledgebaseimporter.knowledgebases

interface KnowledgebaseSource<T : KnownRecord, R : ActionableRecord> {
    val knownKbRecords: List<T>
    val actionableKbRecords: List<R>
}
