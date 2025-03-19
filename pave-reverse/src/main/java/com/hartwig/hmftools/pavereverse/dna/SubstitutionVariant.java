package com.hartwig.hmftools.pavereverse.dna;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.pavereverse.BaseSequenceChange;
import com.hartwig.hmftools.pavereverse.protein.HgvsVariant;
import com.hartwig.hmftools.pavereverse.util.Checks;

public class SubstitutionVariant extends DnaVariant
{
    public final String Ref;
    public final String Alt;

    public SubstitutionVariant(GeneData gene, TranscriptData transcript, HgvsAddress address, String ref, String alt)
    {
        super(gene, transcript, address);
        Preconditions.checkArgument(Checks.isNucleotideSequence(ref));
        Preconditions.checkArgument(Checks.isNucleotideSequence(alt));
        Ref = ref;
        Alt = alt;
    }

    public BaseSequenceChange toGenomicVariant(RefGenomeInterface genome)
    {
        return new BaseSequenceChange(Ref, Alt, chromosome(), Address.toStrandLocation(geneTranscript()));
    }
}
