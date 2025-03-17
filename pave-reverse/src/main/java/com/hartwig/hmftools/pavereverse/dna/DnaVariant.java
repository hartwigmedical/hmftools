package com.hartwig.hmftools.pavereverse.dna;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.pavereverse.BaseSequenceChange;
import com.hartwig.hmftools.pavereverse.protein.HgvsVariant;
import com.hartwig.hmftools.pavereverse.util.Checks;

public class DnaVariant extends HgvsVariant
{
    public final HgvsAddress Address;
    public final String Ref;
    public final String Alt;

    public DnaVariant(GeneData gene, TranscriptData transcript, HgvsAddress address, String ref, String alt)
    {
        super(gene, transcript);
        Preconditions.checkArgument(Checks.isNucleotideSequence(ref));
        Preconditions.checkArgument(Checks.isNucleotideSequence(alt));
        Address = address;
        Ref = ref;
        Alt = alt;
    }

    public BaseSequenceChange toGenomicVariant()
    {
        return new BaseSequenceChange(Ref, Alt, chromosome(), Address.toStrandLocation(geneTranscript()));
    }
}
