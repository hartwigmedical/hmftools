package com.hartwig.hmftools.sage.read;

import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.samtools.CigarHandler;
import com.hartwig.hmftools.common.samtools.CigarTraversal;
import com.hartwig.hmftools.sage.common.IndexedBases;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

// purpose of this class is to handle splits (ie 'N's) and fill them in with wildcards to aid with matching
public class ExpandedBasesFactory
{
    private final int mMaxSkippedReferenceRegions;
    private final int mReferenceRegionReplacementLength;

    public ExpandedBasesFactory(final int maxSkippedReferenceRegions, final int referenceRegionReplacementLength)
    {
        mMaxSkippedReferenceRegions = maxSkippedReferenceRegions;
        mReferenceRegionReplacementLength = referenceRegionReplacementLength;
    }

    public IndexedBases expand(int position, int readIndex, final SAMRecord record)
    {
        final byte[] src = record.getReadBases();
        final AtomicInteger indexAdjustment = new AtomicInteger(0);
        final List<Integer> indexes = Lists.newArrayList();

        final CigarHandler handler = new CigarHandler()
        {
            @Override
            public void handleSkippedReference(@NotNull final SAMRecord record, @NotNull final CigarElement e, final int cigarIndex,
                    final int refPosition)
            {
                if(e.getLength() >= mMaxSkippedReferenceRegions && mReferenceRegionReplacementLength > 0)
                {
                    indexes.add(cigarIndex);
                    if(cigarIndex < readIndex)
                    {
                        indexAdjustment.addAndGet(mReferenceRegionReplacementLength);
                    }
                }
            }
        };

        CigarTraversal.traverseCigar(record, handler);

        if(indexes.isEmpty())
        {
            return new IndexedBases(position, readIndex, record.getReadBases());
        }

        byte[] dest = new byte[src.length + indexes.size() * mReferenceRegionReplacementLength];
        int srcPos = 0;
        int destPos = 0;
        for(Integer index : indexes)
        {
            int length = index + 1 - srcPos;//

            // Copy prior
            System.arraycopy(src, srcPos, dest, destPos, length);
            destPos += length;
            srcPos += length;

            // Create skipped reference substitute
            for(int j = 0; j < mReferenceRegionReplacementLength; j++)
            {
                dest[destPos++] = IndexedBases.MATCH_WILDCARD;
            }
        }

        // Copy remainder
        System.arraycopy(src, srcPos, dest, destPos, src.length - srcPos);

        return new IndexedBases(position, readIndex + indexAdjustment.get(), dest);
    }

}
