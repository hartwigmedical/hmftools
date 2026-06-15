package com.hartwig.hmftools.tars.liftback.rescue;

// Reference-sequence interface for the rescue ref-verify path; decouples production FASTA from tests.
// Coordinates are 1-based inclusive. Returns null for unknown chromosomes or out-of-range queries.
@FunctionalInterface
public interface RefSequenceSource
{
    byte[] getBases(String chromosome, int posStart, int posEnd);
}
