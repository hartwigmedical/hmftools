package com.hartwig.hmftools.pavereverse.dna;

import java.util.Objects;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.pavereverse.BaseSequenceChange;
import com.hartwig.hmftools.pavereverse.util.Checks;

public class DnaVariant
{
    public final GeneData Gene;
    public final TranscriptData Transcript;
    public final int Position;
    public final String Ref;
    public final String Alt;

    public DnaVariant(final GeneData gene, final TranscriptData transcript, final int position, final String ref, final String alt)
    {
        Preconditions.checkArgument(Objects.equals(transcript.GeneId, gene.GeneId));
        Preconditions.checkArgument(Checks.isNucleotideSequence(ref));
        Preconditions.checkArgument(Checks.isNucleotideSequence(alt));
        Gene = gene;
        Transcript = transcript;
        Position = position;
        Ref = ref;
        Alt = alt;
    }

    public BaseSequenceChange toGenomicVariant()
    {
        return null;
    }
}
