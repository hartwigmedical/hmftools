package com.hartwig.hmftools.isofox.expression;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.Strings.appendStrList;

import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

// counts of fragments which could support a set of transcripts and/or unspliced genes

public class CategoryCountsData
{
    private final List<Integer> mTranscripts;
    private final List<String> mUnsplicedGenes;
    private double mFragmentCount;
    private int mLowMapQualFragments;
    private double[] mFragmentCountsByGcRatio;

    // counts by length is only used for expected not actual counts, and is then adjusted by the observed fragment length distribution
    private double[] mFragmentCountsByLength;

    private final String mCombinedKey;

    public CategoryCountsData(final List<Integer> transcripts, final List<String> unsplicedGenes)
    {
        mTranscripts = transcripts;
        mUnsplicedGenes = unsplicedGenes;
        mFragmentCount = 0;
        mLowMapQualFragments = 0;

        mCombinedKey = formTranscriptIds();

        mFragmentCountsByLength = null;
        mFragmentCountsByGcRatio = null;
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
    public final int lowMapQualFragments() { return mLowMapQualFragments; }
    public final double[] fragmentCountsByGcRatio() { return mFragmentCountsByGcRatio; }

    public void addCounts(double count)
    {
        mFragmentCount += count;

        if(count < 1)
            ++mLowMapQualFragments;
    }

    public void adjustCounts(double factor)
    {
        mFragmentCount *= factor;
    }

    public void addGcRatioCounts(double count, final int[] gcRatioIndex, final double[] counts)
    {
        mFragmentCount += count;

        if(count < 1)
            ++mLowMapQualFragments;

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

    // methods for expected counts routine and fragment length distribution adjustments
    public CategoryCountsData(final String categoryStr, int fragLengths)
    {
        mCombinedKey = categoryStr;
        mTranscripts = Lists.newArrayList();
        mUnsplicedGenes = Lists.newArrayList();
        mFragmentCount = 0;
        mFragmentCountsByLength = new double[fragLengths];

        parseCombinedKey();
    }

    public void initialiseLengthCounts(int fragLengths)
    {
        if(fragLengths > 0)
            mFragmentCountsByLength = new double[fragLengths];
    }

    public final double[] fragmentCountsByLength() { return mFragmentCountsByLength; }

    public void addFragLengthCounts(int count, int lengthIndex)
    {
        mFragmentCount += count;
        mFragmentCountsByLength[lengthIndex] += count;
    }

    public void applyFrequencies(final List<Double> lengthFrequencyRates)
    {
        if(lengthFrequencyRates.size() != mFragmentCountsByLength.length)
            return;

        mFragmentCount = 0;

        for(int i = 0; i < mFragmentCountsByLength.length; ++i)
        {
            mFragmentCountsByLength[i] *= lengthFrequencyRates.get(i);
            mFragmentCount += mFragmentCountsByLength[i];
        }
    }

    public void mergeCounts(final CategoryCountsData other)
    {
        mFragmentCount = max(mFragmentCount, other.fragmentCount());

        if(mFragmentCountsByLength != null && other.fragmentCountsByLength() != null)
        {
            for(int i = 0; i < mFragmentCountsByLength.length; ++i)
            {
                mFragmentCountsByLength[i] = max(mFragmentCountsByLength[i], other.fragmentCountsByLength()[i]);
            }
        }
    }

    private static final String DELIM = "-";

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

        StringJoiner sj = new StringJoiner(DELIM);
        items.forEach(x -> sj.add(x));
        return sj.toString();
    }

    private static final String GENE_INDENTIFIER = "ENSG";

    public static boolean hasGeneIdentifier(final String transName) { return transName.startsWith(GENE_INDENTIFIER); }

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

    public String toString()
    {
        return String.format("trans(%d) genes(%d)  key(%s) count(%.1f)",
                mTranscripts.size(), mUnsplicedGenes.size(), mCombinedKey, mFragmentCount);
    }
}
