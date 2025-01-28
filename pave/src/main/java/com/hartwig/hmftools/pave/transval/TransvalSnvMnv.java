package com.hartwig.hmftools.pave.transval;

import java.util.List;
import java.util.Set;

import org.jetbrains.annotations.NotNull;

public class TransvalSnvMnv extends TransvalVariant
{
    @NotNull
    public final String AlternateNucleotides;
    @NotNull
    public final String ReferenceCodon;
    @NotNull
    public final List<String> AlternateCodons;

    public TransvalSnvMnv(
            @NotNull final String transcriptId,
            @NotNull final String chromosome,
            final int position,
            final boolean spansMultipleExons,
            @NotNull final String referenceNucleotides,
            @NotNull final Set<TransvalHotspot> hotspots,
            @NotNull final String alternateNucleotides,
            @NotNull final String referenceCodon,
            @NotNull final List<String> alternateCodons)
    {
        super(transcriptId, chromosome, position, spansMultipleExons, referenceNucleotides, hotspots);
        AlternateNucleotides = alternateNucleotides;
        ReferenceCodon = referenceCodon;
        AlternateCodons = alternateCodons;
    }
}
