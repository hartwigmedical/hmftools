package com.hartwig.hmftools.lilac;

import static com.hartwig.hmftools.common.gene.TranscriptUtils.calcCodingBases;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.bed.ImmutableNamedBed;
import com.hartwig.hmftools.common.genome.bed.NamedBed;

import java.util.List;

public class LociPosition
{
    private final List<TranscriptData> mTranscripts;

    public LociPosition(final List<TranscriptData> transcripts)
    {
        mTranscripts = transcripts;
    }

    public final int calcNucelotideLocus(int position)
    {
        for(TranscriptData transData : mTranscripts)
        {
            if(!positionWithin(position, transData.CodingStart, transData.CodingEnd))
                continue;

            // locus is a zero-based index, so the first coding base has locus of 0
            return calcCodingBases(transData, position).CodingBases - 1;
        }
        return -1;
    }
}
