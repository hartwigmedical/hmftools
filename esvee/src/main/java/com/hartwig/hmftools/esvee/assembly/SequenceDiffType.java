package com.hartwig.hmftools.esvee.assembly;

public enum SequenceDiffType
{
    MATCH, // matches consensus
    BASE, // read has an SNV at this point
    DELETE, // read has a single-base delete
    INSERT, // read has a single-base insert
    REPEAT, // read matches the consensus but for a repeat expansion or contraction
    MISMATCH, // the read does not match the consensus at this location by any of the above types
    INCOMPLETE; // insufficient following bases to determine cause of difference
}
