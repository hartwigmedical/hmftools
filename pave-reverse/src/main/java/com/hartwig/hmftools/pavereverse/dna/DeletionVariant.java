package com.hartwig.hmftools.pavereverse.dna;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.pavereverse.BaseSequenceChange;

public class DeletionVariant extends DnaVariant
{
    private final int mDeletionLength;

    public DeletionVariant(GeneData gene, TranscriptData transcript, HgvsAddress address, String deletedBases)
    {
        super(gene, transcript, address);
        mDeletionLength = deletedBases.length();
    }

    public BaseSequenceChange toGenomicVariant(RefGenomeInterface genome)
    {
        int position = Address.toStrandLocation(geneTranscript());
        int neighbourPosition = position - 1;
        String ref = genome.getBaseString(chromosome(), neighbourPosition, position);

        return new BaseSequenceChange(ref, ref.substring(0, 1), chromosome(), neighbourPosition);
    }
}
