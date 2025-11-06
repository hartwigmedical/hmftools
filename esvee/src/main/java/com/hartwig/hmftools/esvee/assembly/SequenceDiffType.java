package com.hartwig.hmftools.esvee.assembly;

public enum SequenceDiffType
{
    UNSET,
    MATCH, // matches consensus
    BASE, // read has an SNV at this point and subsequent bases match
    DELETE, // read has a single-base delete
    INSERT, // read has a single-base insert
    REPEAT,  // read matches the consensus but for a repeat expansion or contraction
    OTHER; // currently used
    // INCOMPLETE; // insufficient following bases to determine cause of difference
}
