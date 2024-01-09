package com.hartwig.hmftools.esvee.old;

import static java.lang.Math.max;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.read.Read;

import org.apache.commons.lang3.tuple.Pair;

public class SupportedAssembly extends com.hartwig.hmftools.esvee.old.Assembly
{
    private final Map<String, List<ReadSupport>> mReadSupportMap;
    private final List<Read> mReads;
    private final List<ReadSupport> mReadSupportList;

    private final byte[] mBaseQuality;
    private boolean mBaseQualityStale;
    private byte mAverageBaseQuality;

    public SupportedAssembly(final String name, final String assembly)
    {
        super(name, assembly);
        mReadSupportMap = new HashMap<>();
        mReads = Lists.newArrayList();
        mReadSupportList = Lists.newArrayList();

        mBaseQualityStale = true;
        mBaseQuality = new byte[Assembly.length()];
    }

    @Override
    public byte[] getBaseQuality()
    {
        if(mBaseQualityStale)
            recalculateBaseQuality();

        return mBaseQuality;
    }

    @Override
    public byte getAverageBaseQuality()
    {
        if(mBaseQualityStale)
            recalculateBaseQuality();

        return mAverageBaseQuality;
    }

    public boolean tryAddSupport(final SupportChecker checker, final Read read, final int suggestedIndex)
    {
        if(containsSupport(read))
            return true;

        if(checker.WeakSupport.supportsAt(this, read, suggestedIndex))
        {
            addEvidenceAt(read, suggestedIndex);
            return true;
        }
        else
            return tryAddSupport(checker, read);
    }

    public boolean tryAddSupport(final SupportChecker checker, final Read read)
    {
        Integer index = checker.WeakSupport.supportIndex(this, read);

        if(index == null)
            return false;

        addEvidenceAt(read, index);
        return true;
    }

    public void addEvidenceAt(final Read read, final int supportIndex)
    {
        if(containsSupport(read, supportIndex))
            return;

        mBaseQualityStale = true;

        ReadSupport readSupport = new ReadSupport(read, supportIndex);
        mReadSupportList.add(readSupport);
        mReads.add(read);

        List<ReadSupport> support = mReadSupportMap.get(read.getName());

        if(support == null)
        {
            support = Lists.newArrayListWithCapacity(2);
            mReadSupportMap.put(read.getName(), support);
        }

        support.add(readSupport);

        // mSupport.compute(read.getName(), (k, existing) -> new SupportEntry(read, supportIndex, existing));
    }

    public List<Read> supportingReads() { return mReads; }
    public List<ReadSupport> readSupport() { return mReadSupportList; }
    public Map<String, List<ReadSupport>> readSupportMap() { return mReadSupportMap; }

    public Set<String> getSupportReadNames()
    {
        return mReadSupportMap.keySet();
    }

    public int readSupportCount() { return mReads.size(); }

    public int fragmentSupportCount() { return mReadSupportMap.size(); }

    public boolean containsSupport(final Read read)
    {
        return mReads.contains(read);
    }

    public boolean containsSupport(final Read read, final int index)
    {
        List<ReadSupport> support = mReadSupportMap.get(read.getName());
        return support != null && support.stream().anyMatch(x -> x.Read == read && x.Index == index);
    }

    public List<ReadSupport> getReadSupport(final String readId)
    {
        return mReadSupportMap.get(readId);
    }

    public int getSupportIndex(final Read read)
    {
        List<ReadSupport> support = mReadSupportMap.get(read.getName());
        if(support == null)
            return -1;

        for(ReadSupport readSupport : support)
        {
            if(readSupport.Read == read)
                return readSupport.Index;
        }

        return -1;
        // throw new IllegalStateException(String.format("Record %s does not support assembly %s", read, getName()));
    }

    public Pair<int[], int[]> computeBaseSupportAndContradiction()
    {
        final int[] baseQualitySupporting = new int[Assembly.length()];
        final int[] baseQualityContradicting = new int[Assembly.length()];

        for(ReadSupport readSupport : mReadSupportList)
        {
            Read read = readSupport.Read;

            int assemblyOffset = max(readSupport.Index, 0);
            int readOffset = readSupport.Index < 0 ? -readSupport.Index : 0;

            int remainingAssemblyLength = Assembly.length() - assemblyOffset;
            int remainingReadLength = read.getLength() - readOffset;
            int supportLength = Math.min(remainingAssemblyLength, remainingReadLength);

            for(int i = 0; i < supportLength; i++)
            {
                byte assemblyBase = AssemblyBases[assemblyOffset + i];
                byte evidenceBase = read.getBases()[readOffset + i];
                byte evidenceQuality = read.getBaseQuality()[readOffset + i];
                if(assemblyBase == evidenceBase)
                    baseQualitySupporting[assemblyOffset + i] += evidenceQuality;
                else
                    baseQualityContradicting[assemblyOffset + i] += evidenceQuality;
            }
        }

        return Pair.of(baseQualitySupporting, baseQualityContradicting);
    }

    public void recalculateBaseQuality()
    {
        if(!mBaseQualityStale)
            return;

        final byte[] maxBaseQualitySupporting = new byte[Assembly.length()];
        final int[] baseQualitySupporting = new int[Assembly.length()];
        final int[] baseQualityContradicting = new int[Assembly.length()];

        for(ReadSupport readSupport : mReadSupportList)
        {
            Read read = readSupport.Read;

            final int assemblyOffset = max(readSupport.Index, 0);
            final int readOffset = readSupport.Index < 0 ? -readSupport.Index : 0;

            final int remainingAssemblyLength = Assembly.length() - assemblyOffset;
            final int remainingReadLength = read.getLength() - readOffset;
            final int supportLength = Math.min(remainingAssemblyLength, remainingReadLength);
            for(int i = 0; i < supportLength; i++)
            {
                final byte assemblyBase = AssemblyBases[assemblyOffset + i];
                final byte evidenceBase = read.getBases()[readOffset + i];
                final byte evidenceQuality = read.getBaseQuality()[readOffset + i];
                if(assemblyBase == evidenceBase)
                {
                    maxBaseQualitySupporting[assemblyOffset + i] =
                            (byte) max(maxBaseQualitySupporting[assemblyOffset + i], evidenceQuality);
                    baseQualitySupporting[assemblyOffset + i] += evidenceQuality;
                }
                else
                {
                    baseQualityContradicting[assemblyOffset + i] += evidenceQuality;
                }
            }
        }

        int baseQualitySum = 0;
        for(int i = 0; i < Assembly.length(); i++)
        {
            final int maxSupportedQuality = maxBaseQualitySupporting[i];
            final int supportingQuality = baseQualitySupporting[i];
            final int contraryQuality = baseQualityContradicting[i];

            if(supportingQuality == 0)
            {
                mBaseQuality[i] = 0;
                continue;
            }

            final int quality = max(0, maxSupportedQuality * max(supportingQuality - contraryQuality, 0) / supportingQuality);
            mBaseQuality[i] = (byte) quality;
            baseQualitySum += quality;
        }
        mAverageBaseQuality = (byte) (baseQualitySum / getLength());
        mBaseQualityStale = false;
    }

    @Override
    public String toString()
    {
        return String.format("%s (%s support)", Name, readSupportCount());
    }
}

