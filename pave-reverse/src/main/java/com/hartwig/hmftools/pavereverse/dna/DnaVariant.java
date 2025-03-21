package com.hartwig.hmftools.pavereverse.dna;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.pavereverse.BaseSequenceChange;
import com.hartwig.hmftools.pavereverse.protein.HgvsVariant;

import org.apache.commons.lang3.tuple.Pair;

public abstract class DnaVariant extends HgvsVariant
{
    public final HgvsAddress AddressOfChangeStart;
    public final HgvsAddress AddressOfChangeEnd;

    public DnaVariant(GeneData gene, TranscriptData transcript, HgvsAddress addressOfChangeStart, HgvsAddress addressOfChangeEnd)
    {
        super(gene, transcript);
        AddressOfChangeStart = addressOfChangeStart;
        AddressOfChangeEnd = addressOfChangeEnd;
    }

    public DnaVariant(GeneData gene, TranscriptData transcript, HgvsAddress addressOfChangeStart)
    {
        this(gene, transcript, addressOfChangeStart, addressOfChangeStart);
    }

    public abstract BaseSequenceChange toGenomicVariant(RefGenomeInterface genome);

    Pair<Integer, Integer> getAbsoluteLocationsOfChange()
    {
        int start = AddressOfChangeStart.toStrandLocation(geneTranscript());
        int stop = AddressOfChangeEnd.toStrandLocation(geneTranscript());
        if(reverseStrand())
        {
            return Pair.of(stop, start);
        }
        return Pair.of(start, stop);
    }
}
