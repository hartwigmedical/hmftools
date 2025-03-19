package com.hartwig.hmftools.pavereverse.dna;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.pavereverse.BaseSequenceChange;
import com.hartwig.hmftools.pavereverse.protein.HgvsVariant;

public abstract class DnaVariant extends HgvsVariant
{
    public final HgvsAddress Address;

    public DnaVariant(GeneData gene, TranscriptData transcript, HgvsAddress address)
    {
        super(gene, transcript);
        Address = address;
    }

    public abstract BaseSequenceChange toGenomicVariant(RefGenomeInterface genome);
}
