package com.hartwig.hmftools.esvee.models;

import java.util.AbstractMap;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.function.Supplier;

import com.hartwig.hmftools.esvee.read.Read;
import com.hartwig.hmftools.esvee.util.SizedIterable;
import com.hartwig.hmftools.esvee.assembly.SupportChecker;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.Nullable;

public class SupportedAssembly extends com.hartwig.hmftools.esvee.common.Assembly
{
    private final Map<String,SupportEntry> mSupport;
    private int mSupportCount;

    private final byte[] mBaseQuality;
    private boolean mBaseQualityStale;
    private byte mAverageBaseQuality;

    public SupportedAssembly(final String name, final String assembly)
    {
        super(name, assembly);
        mSupport = new HashMap<>();
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

    public Pair<int[], int[]> computeBaseSupportAndContradiction()
    {
        final int[] baseQualitySupporting = new int[Assembly.length()];
        final int[] baseQualityContradicting = new int[Assembly.length()];

        for(Map.Entry<Read, Integer> entry : getSupport())
        {
            final Read read = entry.getKey();
            final int supportIndex = entry.getValue();

            final int assemblyOffset = Math.max(supportIndex, 0);
            final int readOffset = supportIndex < 0 ? -supportIndex : 0;

            final int remainingAssemblyLength = Assembly.length() - assemblyOffset;
            final int remainingReadLength = read.getLength() - readOffset;
            final int supportLength = Math.min(remainingAssemblyLength, remainingReadLength);
            for(int i = 0; i < supportLength; i++)
            {
                final byte assemblyBase = AssemblyBases[assemblyOffset + i];
                final byte evidenceBase = read.getBases()[readOffset + i];
                final byte evidenceQuality = read.getBaseQuality()[readOffset + i];
                if(assemblyBase == evidenceBase)
                    baseQualitySupporting[assemblyOffset + i] += evidenceQuality;
                else
                    baseQualityContradicting[assemblyOffset + i] += evidenceQuality;
            }
        }

        return Pair.of(baseQualitySupporting, baseQualityContradicting);
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
        @Nullable
        final Integer index = checker.WeakSupport.supportIndex(this, read);

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
        mSupport.compute(read.getName(), (k, existing) -> new SupportEntry(read, supportIndex, existing));
        mSupportCount++;
    }

    public Set<String> getSupportFragments() { return mSupport.keySet(); }
    
    public int supportCount() { return mSupportCount; }

    public boolean containsSupport(final Read read)
    {
        @Nullable
        SupportEntry existingEntry = mSupport.get(read.getName());
        while(existingEntry != null)
        {
            if(existingEntry.Read.equals(read))
                return true;
            existingEntry = existingEntry.Next;
        }

        return false;
    }

    public boolean containsSupport(final Read read, final int index)
    {
        @Nullable
        SupportEntry existingEntry = mSupport.get(read.getName());
        while(existingEntry != null)
        {
            if(existingEntry.Read.equals(read) && existingEntry.SupportIndex == index)
                return true;
            existingEntry = existingEntry.Next;
        }

        return false;
    }

    public SizedIterable<Read> getSupportRecords()
    {
        return SizedIterable.create(supportCount(), this::supportIterator, entry -> entry.Read);
    }

    public SizedIterable<Map.Entry<Read, Integer>> getSupport()
    {
        return SizedIterable.create(supportCount(), this::supportIterator,
                entry -> new AbstractMap.SimpleEntry<>(entry.Read, entry.SupportIndex));
    }

    public SizedIterable<Map.Entry<Read, Integer>> getSupport(final String fragmentName)
    {
        return new SizedIterable<>(-1, new Supplier<Iterator<Map.Entry<Read, Integer>>>()
        {
            @Override
            public Iterator<Map.Entry<Read, Integer>> get()
            {
                return new Iterator<>()
                {
                    @Nullable
                    SupportEntry supportEntry = mSupport.get(fragmentName);

                    @Override
                    public boolean hasNext()
                    {
                        return supportEntry != null;
                    }

                    @Override
                    public Map.Entry<Read, Integer> next()
                    {
                        assert supportEntry != null;
                        final Map.Entry<Read, Integer> result = Map.entry(supportEntry.Read, supportEntry.SupportIndex);
                        supportEntry = supportEntry.Next;
                        return result;
                    }
                };
            }
        });
    }

    public int getSupportIndex(final Read read)
    {
        SupportEntry entry = mSupport.get(read.getName());
        while(entry != null)
        {
            if(entry.Read.equals(read))
                return entry.SupportIndex;
            entry = entry.Next;
        }
        throw new IllegalStateException(String.format("Record %s does not support assembly %s", read, getName()));
    }

    private Iterator<SupportEntry> supportIterator()
    {
        return new Iterator<>()
        {
            final Iterator<SupportEntry> entryIterator = mSupport.values().iterator();
            @Nullable
            SupportEntry current = null;

            @Override
            public boolean hasNext()
            {
                if(current != null)
                    return true;

                return entryIterator.hasNext();
            }

            @Override
            public SupportEntry next()
            {
                if(current == null)
                    current = entryIterator.next();

                final SupportEntry toReturn = current;
                current = current.Next;
                return toReturn;
            }
        };
    }

    public void recalculateBaseQuality()
    {
        if(!mBaseQualityStale)
            return;

        final byte[] maxBaseQualitySupporting = new byte[Assembly.length()];
        final int[] baseQualitySupporting = new int[Assembly.length()];
        final int[] baseQualityContradicting = new int[Assembly.length()];

        for(Map.Entry<Read, Integer> entry : getSupport())
        {
            final Read read = entry.getKey();
            final int supportIndex = entry.getValue();

            final int assemblyOffset = Math.max(supportIndex, 0);
            final int readOffset = supportIndex < 0 ? -supportIndex : 0;

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
                    maxBaseQualitySupporting[assemblyOffset + i] = (byte) Math.max(maxBaseQualitySupporting[assemblyOffset + i], evidenceQuality);
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

            final int quality = Math.max(0, maxSupportedQuality * Math.max(supportingQuality - contraryQuality, 0) / supportingQuality);
            mBaseQuality[i] = (byte) quality;
            baseQualitySum += quality;
        }
        mAverageBaseQuality = (byte) (baseQualitySum / getLength());
        mBaseQualityStale = false;
    }

    @Override
    public String toString()
    {
        return String.format("%s (%s support)", Name, supportCount());
    }

    private static class SupportEntry
    {
        public final Read Read;
        public final int SupportIndex;
        @Nullable
        public final SupportEntry Next;

        private SupportEntry(final Read read, final int supportIndex, @Nullable final SupportEntry next)
        {
            Read = read;
            SupportIndex = supportIndex;
            Next = next;
        }

        @Override
        public String toString()
        {
            final String core = String.valueOf(SupportIndex);
            if(Next == null)
                return core;
            else
                return core + "&" + Next;
        }
    }
}
