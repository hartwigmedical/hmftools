package com.hartwig.hmftools.pavereverse.dna;

import static com.hartwig.hmftools.pavereverse.ReversePaveConfig.RPV_LOGGER;

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
    private final int mAbsoluteLocationOfStart;
    private final int mAbsoluteLocationOfEnd;

    public DnaVariant(GeneData gene, TranscriptData transcript, HgvsAddress addressOfChangeStart, HgvsAddress addressOfChangeEnd)
    {
        super(gene, transcript);
        AddressOfChangeStart = addressOfChangeStart;
        AddressOfChangeEnd = addressOfChangeEnd;
        String warningForStart = AddressOfChangeStart.consistencyWarnings(geneTranscript());
        if(warningForStart != null)
        {
            RPV_LOGGER.warn(warningForStart);
        }
        String warningForEnd = AddressOfChangeStart.consistencyWarnings(geneTranscript());
        if(warningForEnd != null)
        {
            if(!warningForEnd.equals(warningForStart))
            {
                RPV_LOGGER.warn(warningForEnd);
            }
        }
        int start = AddressOfChangeStart.toStrandLocation(geneTranscript());
        int stop = AddressOfChangeEnd.toStrandLocation(geneTranscript());
        if(reverseStrand())
        {
            mAbsoluteLocationOfStart = stop;
            mAbsoluteLocationOfEnd = start;
        }
        else
        {
            mAbsoluteLocationOfStart = start;
            mAbsoluteLocationOfEnd = stop;
        }
        if(mAbsoluteLocationOfStart > mAbsoluteLocationOfEnd)
        {
            RPV_LOGGER.warn("Start is beyond end: " + mAbsoluteLocationOfStart + " > " + mAbsoluteLocationOfEnd);
        }
    }

    public abstract BaseSequenceChange toGenomicVariant(RefGenomeInterface genome);

    Pair<Integer, Integer> getAbsoluteLocationsOfChange()
    {
        return Pair.of(mAbsoluteLocationOfStart, mAbsoluteLocationOfEnd);
    }
}
