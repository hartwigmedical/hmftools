package com.hartwig.hmftools.tars.liftback.rescue;

// Splice-motif tier of a candidate junction, ordered weakest to strongest. Declaration order is meaningful:
// callers rank candidates with compareTo (a stronger motif outranks a weaker one).
public enum Tier
{
    NONE,            // neither donor nor acceptor matches a known splice motif
    SEMI_CANONICAL,  // GC-AG / AT-AC (and reverse-complement equivalents)
    CANONICAL,       // GT-AG (~99% of splice sites)
    ANNOTATED        // matches an annotated junction from the sidecar
}
