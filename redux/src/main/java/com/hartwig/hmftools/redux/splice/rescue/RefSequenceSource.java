package com.hartwig.hmftools.redux.splice.rescue;

// Minimal reference-sequence interface so the rescue resolver's ref-verify path can be tested
// without a real FASTA file. Production implementation wraps htsjdk's IndexedFastaSequenceFile;
// tests construct an in-memory map keyed by chromosome.
//
// Coordinates are 1-based inclusive on both ends. Returns null if the requested range falls
// outside the source's known sequence (e.g. chromosome unknown, or range past end of contig).
@FunctionalInterface
public interface RefSequenceSource
{
    byte[] getBases(String chromosome, int posStart, int posEnd);
}
