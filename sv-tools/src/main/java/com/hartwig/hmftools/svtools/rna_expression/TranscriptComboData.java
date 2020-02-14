package com.hartwig.hmftools.svtools.rna_expression;

import static com.hartwig.hmftools.linx.analysis.SvUtilities.appendStrList;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.TC_LONG;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.TC_MAX;

import java.util.List;

public class TranscriptComboData
{
    private final List<String> mTranscripts;
    private final int[] mCounts;

    public TranscriptComboData(final List<String> transcripts)
    {
        mTranscripts = transcripts;
        mCounts = new int[TC_MAX];
    }

    public final List<String> getTranscripts() { return mTranscripts; }

    public String getTranscriptStr(int maxCount)
    {
        if(mTranscripts.size() > maxCount)
            return appendStrList(mTranscripts.subList(0, 10), ';');
        else
            return appendStrList(mTranscripts, ';');
    }

    public boolean matches(final List<String> transcripts)
    {
        if(mTranscripts.size() != transcripts.size())
            return false;

        for(final String trans : transcripts)
        {
            if(!mTranscripts.stream().anyMatch(x -> x.equals(trans)))
                return false;
        }

        return true;
    }

    public final int[] getCounts() { return mCounts; }


}
