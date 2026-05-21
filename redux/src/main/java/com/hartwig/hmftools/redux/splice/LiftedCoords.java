package com.hartwig.hmftools.redux.splice;

// minimal lift result — the translated genomic coordinates of a single alt-contig alignment. Used by
// SaTagRewriter, where mapQ / NM pass through unchanged so the full LiftBackResult machinery isn't needed.
public record LiftedCoords(
        String chromosome,
        int position,
        String cigarString)
{
}
