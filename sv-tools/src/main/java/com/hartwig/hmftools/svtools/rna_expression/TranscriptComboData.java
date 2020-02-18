package com.hartwig.hmftools.svtools.rna_expression;

import static com.hartwig.hmftools.linx.analysis.SvUtilities.appendStr;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.appendStrList;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.TC_LONG;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.TC_MAX;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.TC_SHORT;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.TC_SPLICED;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.TC_UNSPLICED;

import java.util.Arrays;
import java.util.List;

import com.google.common.collect.Lists;

public class TranscriptComboData
{
    private final List<String> mTranscripts;
    private final int[] mCounts;

    private final String mTranscriptsKey;

    public TranscriptComboData(final List<String> transcripts)
    {
        mTranscripts = transcripts;
        mCounts = new int[TC_MAX];

        mTranscriptsKey = formTranscriptIds(transcripts);
    }

    public final List<String> getTranscripts() { return mTranscripts; }
    public final String getTranscriptsKey() { return mTranscriptsKey; }

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
    public final int getCount(int type) { return type < mCounts.length ? mCounts[type] : 0; }
    public int totalCount() { return Arrays.stream(mCounts).sum(); }

    public static TranscriptComboData findMatchingData(final List<String> transcripts, final List<TranscriptComboData> dataList)
    {
        return dataList.stream().filter(x -> x.matches(transcripts)).findFirst().orElse(null);
    }

    public static String formTranscriptIds(final List<String> transcripts)
    {
        // convert into an order list of ints
        List<Integer> transIds = Lists.newArrayList();

        try
        {
            for (String trans : transcripts)
            {
                // attempt conversion to int
                int transInt = Integer.parseInt(trans.replaceAll("[A-Z]*", ""));

                int index = 0;
                while (index < transIds.size())
                {
                    if (transInt < transIds.get(index))
                        break;

                    ++index;
                }

                transIds.add(index, transInt);
            }

            String transIdStr = "";
            for (Integer transId : transIds)
            {
                transIdStr = appendStr(transIdStr, String.valueOf(transId), '-');
            }

            return transIdStr;
        }
        catch (Exception e)
        {
            return appendStrList(transcripts, '-');
        }
    }

    public String toString()
    {
        return String.format("ids(%d: {}) counts(spl=%d short=%d long=%d unspl=%d)",
                mTranscripts.size(), mTranscriptsKey, mCounts[TC_SPLICED], mCounts[TC_SHORT],mCounts[TC_LONG], mCounts[TC_UNSPLICED]);
    }


}
