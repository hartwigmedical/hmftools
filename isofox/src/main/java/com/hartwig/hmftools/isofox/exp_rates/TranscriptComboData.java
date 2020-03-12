package com.hartwig.hmftools.isofox.exp_rates;

import static com.hartwig.hmftools.linx.analysis.SvUtilities.appendStr;

import java.util.Arrays;
import java.util.List;

import com.google.common.collect.Lists;

// counts of fragments which could support a set of transcripts and/or unspliced genes

public class TranscriptComboData
{
    private final List<Integer> mTranscripts;
    private final List<String> mUnsplicedGenes;
    private int mFragmentCount;

    private final String mCombinedKey;

    public TranscriptComboData(final List<Integer> transcripts, final List<String> unsplicedGenes)
    {
        mTranscripts = transcripts;
        mUnsplicedGenes = unsplicedGenes;
        mFragmentCount = 0;

        mCombinedKey = formTranscriptIds();
    }

    public final List<Integer> transcriptIds() { return mTranscripts; }
    public final List<String> unsplicedGeneIds() { return mUnsplicedGenes; }
    public final String combinedKey() { return mCombinedKey; }

    public String impliedType()
    {
        if(mUnsplicedGenes.isEmpty())
            return "SPLICED/LONG";
        else if(mTranscripts.isEmpty())
            return "UNSPLICED";
        else
            return "SHORT";
    }

    public boolean matches(final List<Integer> transcripts, final List<String> unsplicedGenes)
    {
        if(mTranscripts.size() != transcripts.size() || mUnsplicedGenes.size() != unsplicedGenes.size())
            return false;

        for(int transId : transcripts)
        {
            if(!mTranscripts.stream().anyMatch(x -> x == transId))
                return false;
        }

        for(String geneId : unsplicedGenes)
        {
            if(!mUnsplicedGenes.stream().anyMatch(x -> x.equals(geneId)))
                return false;
        }

        return true;
    }

    public final int fragmentCount() { return mFragmentCount; }

    public void addCounts(int count)
    {
        mFragmentCount += count;
    }

    private String formTranscriptIds()
    {
        // convert into an order list of ints
        List<Integer> transIds = Lists.newArrayList();

        for (Integer transId : mTranscripts)
        {
            int index = 0;
            while (index < transIds.size())
            {
                if (transId < transIds.get(index))
                    break;

                ++index;
            }

            transIds.add(index, transId);
        }

        String combinedKey = "";
        for (Integer transId : transIds)
        {
            combinedKey = appendStr(combinedKey, String.valueOf(transId), '-');
        }

        for(String geneId : mUnsplicedGenes)
        {
            combinedKey = appendStr(combinedKey, geneId, '-');
        }

        return combinedKey;
    }

    public String toString()
    {
        return String.format("trans(%d) genes(%d)  key(%s) count(%d)",
                mTranscripts.size(), mUnsplicedGenes.size(), mCombinedKey, mFragmentCount);
    }


}
