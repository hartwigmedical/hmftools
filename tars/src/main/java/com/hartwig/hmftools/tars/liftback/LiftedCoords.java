package com.hartwig.hmftools.tars.liftback;

// Translated genomic coordinates for a single alt-contig alignment. Used by SaTagRewriter (mapQ/NM pass through).
public record LiftedCoords(
        String chromosome,
        int position,
        String cigarString)
{
}
