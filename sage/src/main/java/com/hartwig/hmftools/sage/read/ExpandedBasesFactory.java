package com.hartwig.hmftools.sage.read;

import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.samtools.CigarHandler;
import com.hartwig.hmftools.common.samtools.CigarTraversal;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class ExpandedBasesFactory {

    private final int maxSkippedReferenceRegions;
    private final int referenceRegionReplacementLength;

    public ExpandedBasesFactory(final int maxSkippedReferenceRegions, final int referenceRegionReplacementLength) {
        this.maxSkippedReferenceRegions = maxSkippedReferenceRegions;
        this.referenceRegionReplacementLength = referenceRegionReplacementLength;
    }

    @NotNull
    public IndexedBases expand(int position, int readIndex, SAMRecord record) {
        final byte[] src = record.getReadBases();
        final AtomicInteger indexAdjustment = new AtomicInteger(0);
        final List<Integer> indexes = Lists.newArrayList();

        final CigarHandler handler = new CigarHandler() {
            @Override
            public void handleSkippedReference(@NotNull final SAMRecord record, @NotNull final CigarElement e, final int cigarIndex,
                    final int refPosition) {
                if (e.getLength() >= maxSkippedReferenceRegions && referenceRegionReplacementLength > 0) {
                    indexes.add(cigarIndex);
                    if (cigarIndex < readIndex) {
                        indexAdjustment.addAndGet(referenceRegionReplacementLength);
                    }
                }
            }
        };

        CigarTraversal.traverseCigar(record, handler);

        if (indexes.isEmpty()) {
            return new IndexedBases(position, readIndex, record.getReadBases());
        }

        byte[] dest = new byte[src.length + indexes.size() * referenceRegionReplacementLength];
        int srcPos = 0;
        int destPos = 0;
        for (Integer index : indexes) {
            int length = index + 1 - srcPos;//

            // Copy prior
            System.arraycopy(src, srcPos, dest, destPos, length);
            destPos += length;
            srcPos += length;

            // Create skipped reference substitute
            for (int j = 0; j < referenceRegionReplacementLength; j++) {
                dest[destPos++] = IndexedBases.MATCH_WILDCARD;
            }
        }

        // Copy remainder
        System.arraycopy(src, srcPos, dest, destPos, src.length - srcPos);

        return new IndexedBases(position, readIndex + indexAdjustment.get(), dest);
    }

}
