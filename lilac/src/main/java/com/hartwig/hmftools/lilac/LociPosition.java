package com.hartwig.hmftools.lilac;

import static com.hartwig.hmftools.common.fusion.TranscriptUtils.calcCodingBases;

import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
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
            if(position < transData.CodingStart || position > transData.CodingEnd)
                continue;

            // locus is a zero-based index, so the first coding base has locus of 0
            return calcCodingBases(transData, position).CodingBases - 1;
        }
        return -1;
    }

}
