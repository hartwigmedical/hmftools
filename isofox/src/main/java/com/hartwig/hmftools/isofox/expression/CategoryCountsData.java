package com.hartwig.hmftools.isofox.expression;

import static com.hartwig.hmftools.common.utils.Strings.appendStrList;
import static com.hartwig.hmftools.isofox.expression.ExpectedRatesGenerator.FL_FREQUENCY;

import java.util.List;

import com.google.common.collect.Lists;

// counts of fragments which could support a set of transcripts and/or unspliced genes

public class CategoryCountsData
{
    private final List<Integer> mTranscripts;
    private final List<String> mUnsplicedGenes;
    private double mFragmentCount;
    private int[] mFragmentCountsByLength;
    private double[] mFragmentCountsByGcRatio;
    private int mLowMapQualityCount;

    private final String mCombinedKey;

    public CategoryCountsData(final List<Integer> transcripts, final List<String> unsplicedGenes)
    {
        mTranscripts = transcripts;
        mUnsplicedGenes = unsplicedGenes;
        mFragmentCount = 0;

        mCombinedKey = formTranscriptIds();

        mFragmentCountsByLength = null;
        mFragmentCountsByGcRatio = null;
        mLowMapQualityCount = 0;
    }

    public CategoryCountsData(final String categoryStr, int fragLengths)
    {
        mCombinedKey = categoryStr;
        mTranscripts = Lists.newArrayList();
        mUnsplicedGenes = Lists.newArrayList();
        mFragmentCount = 0;
        mFragmentCountsByLength = new int[fragLengths];
        mLowMapQualityCount = 0;

        parseCombinedKey();
    }

    public void initialiseLengthCounts(int fragLengths)
    {
        if(fragLengths > 0)
            mFragmentCountsByLength = new int[fragLengths];
    }

    public void initialiseGcRatioCounts(int gcRatioBuckets)
    {
        if(gcRatioBuckets > 0)
            mFragmentCountsByGcRatio = new double[gcRatioBuckets];
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

    public final double fragmentCount() { return mFragmentCount; }
    public final int[] fragmentCountsByLength() { return mFragmentCountsByLength; }
    public final double[] fragmentCountsByGcRatio() { return mFragmentCountsByGcRatio; }

    public void addCounts(int count)
    {
        mFragmentCount += count;
    }
    public void adjustCounts(double factor)
    {
        mFragmentCount *= factor;
    }

    public void addFragLengthCounts(int count, int lengthIndex)
    {
        mFragmentCount += count;
        mFragmentCountsByLength[lengthIndex] += count;
    }

    public void applyFrequencies(final List<int[]> lengthFrequencies)
    {
        mFragmentCount = 0;

        for(int i = 0; i < mFragmentCountsByLength.length; ++i)
        {
            mFragmentCountsByLength[i] *= lengthFrequencies.get(i)[FL_FREQUENCY];
            mFragmentCount += mFragmentCountsByLength[i];
        }
    }

    public void addGcRatioCounts(int count, final int[] gcRatioIndex, final double[] counts)
    {
        mFragmentCount += count;

        if(gcRatioIndex != null && counts != null)
        {
            for (int i = 0; i < gcRatioIndex.length; ++i)
            {
                if (gcRatioIndex[i] >= 0)
                    mFragmentCountsByGcRatio[gcRatioIndex[i]] += counts[i];
            }
        }
    }

    public void applyGcAdjustments(final double[] gcAdjustments)
    {
        mFragmentCount = 0;

        for(int i = 0; i < mFragmentCountsByGcRatio.length; ++i)
        {
            mFragmentCountsByGcRatio[i] *= gcAdjustments[i];
            mFragmentCount += mFragmentCountsByGcRatio[i];
        }
    }

    private static final char DELIM = '-';

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

        List<String> items = Lists.newArrayList();

        for (Integer transId : transIds)
        {
            items.add(String.valueOf(transId));
        }

        items.addAll(mUnsplicedGenes);

        return appendStrList(items, DELIM);
    }

    private static final String GENE_INDENTIFIER = "ENSG";

    private void parseCombinedKey()
    {
        String[] items = mCombinedKey.split(String.valueOf(DELIM));

        for(int i = 0; i < items.length; ++i)
        {
            if(items[i].contains(GENE_INDENTIFIER))
            {
                mUnsplicedGenes.add(items[i]);
            }
            else
            {
                mTranscripts.add(Integer.parseInt(items[i]));
            }
        }
    }

    public int lowMapQualityCount() { return mLowMapQualityCount; }
    public void addLowMapQualityCount() { mLowMapQualityCount++; }

    public String toString()
    {
        return String.format("trans(%d) genes(%d)  key(%s) count(%.0f)",
                mTranscripts.size(), mUnsplicedGenes.size(), mCombinedKey, mFragmentCount);
    }


}
