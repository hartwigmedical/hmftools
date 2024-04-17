package com.hartwig.hmftools.sage.common;

import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.CigarHandler;
import com.hartwig.hmftools.common.bam.CigarTraversal;
import com.hartwig.hmftools.sage.old.IndexedBases;
import com.hartwig.hmftools.sage.evidence.ReadIndexBases;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public final class SplitReadUtils
{
    public static final int MAX_SKIPPED_REFERENCE_REGIONS = 50;

    // purpose of this class is to handle splits (ie 'N's) and fill them in with wildcards to aid with matching
    public static ReadIndexBases expandSplitRead(int readIndex, final SAMRecord record)
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
                if(e.getLength() >= MAX_SKIPPED_REFERENCE_REGIONS)
                {
                    indexes.add(cigarIndex);
                    if(cigarIndex < readIndex)
                    {
                        indexAdjustment.addAndGet(MAX_SKIPPED_REFERENCE_REGIONS);
                    }
                }
            }
        };

        CigarTraversal.traverseCigar(record, handler);

        if(indexes.isEmpty())
        {
            return new ReadIndexBases(readIndex, record.getReadBases());
        }

        byte[] dest = new byte[src.length + indexes.size() * MAX_SKIPPED_REFERENCE_REGIONS];
        int srcPos = 0;
        int destPos = 0;
        for(Integer index : indexes)
        {
            int length = index + 1 - srcPos;

            // Copy prior
            System.arraycopy(src, srcPos, dest, destPos, length);
            destPos += length;
            srcPos += length;

            // Create skipped reference substitute
            for(int j = 0; j < MAX_SKIPPED_REFERENCE_REGIONS; j++)
            {
                dest[destPos++] = IndexedBases.MATCH_WILDCARD;
            }
        }

        // Copy remainder
        System.arraycopy(src, srcPos, dest, destPos, src.length - srcPos);

        return new ReadIndexBases(readIndex + indexAdjustment.get(), dest);
    }

}
