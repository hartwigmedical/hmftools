package com.hartwig.hmftools.sage.quality;

import java.util.Collection;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.sage.read.IndexedBases;
import com.hartwig.hmftools.sage.ref.RefSequence;
import com.hartwig.hmftools.sage.sam.CigarHandler;
import com.hartwig.hmftools.sage.sam.CigarTraversal;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class QualityCounter implements CigarHandler {

    private static final CigarElement SINGLE = new CigarElement(1, CigarOperator.M);

    private final IndexedBases refGenome;
    private final GenomeRegion bounds;

    private final Set<Integer> indelPositions = Sets.newHashSet();
    private final Map<QualityRecord, QualityCount> qualityMap = Maps.newHashMap();

    public QualityCounter(final RefSequence refGenome, final GenomeRegion bounds) {
        this.refGenome = refGenome.alignment();
        this.bounds = bounds;
    }

    public void processRecord(@NotNull final SAMRecord record) {
        CigarTraversal.traverseCigar(record, this);
    }

    @NotNull
    public Collection<QualityCount> counts() {
        final Set<QualityCount> result = Sets.newHashSet();

        for (QualityCount value : qualityMap.values()) {
            if (isValid(value)) {
                result.add(value);
            }
        }

        return QualityGrouping.groupWithoutPosition(result);
    }

    private boolean isValid(QualityCount count) {

        if (indelPositions.contains(count.position())) {
            return false;
        }
        return count.ref() == count.alt() || count.count() <= 3;
    }

    @Override
    public void handleInsert(@NotNull final SAMRecord record, @NotNull final CigarElement e, final int readIndex, final int refPos) {
        // Need to add one because indel is actually AFTER this by convention
        indelPositions.add(refPos + 1);
        handleAlignment(record, SINGLE, readIndex, refPos);
    }

    @Override
    public void handleDelete(@NotNull final SAMRecord record, @NotNull final CigarElement e, final int readIndex, final int refPos) {
        indelPositions.add(refPos + 1);
        handleAlignment(record, SINGLE, readIndex, refPos);
    }

    @Override
    public void handleAlignment(@NotNull final SAMRecord record, @NotNull final CigarElement e, final int readIndex, final int refPos) {
        for (int i = 0; i < e.getLength(); i++) {
            int index = readIndex + i;
            int position = refPos + i;

            if (position > bounds.end()) {
                return;
            }

            if (position < bounds.start()) {
                continue;
            }

            byte ref = refGenome.base(position);
            byte alt = record.getReadBases()[index];
            byte quality = record.getBaseQualities()[index];
            boolean firstOfPairFlag = record.getFirstOfPairFlag();

            final QualityRecord key = ImmutableQualityRecord.builder()
                    .position(position)
                    .ref(ref)
                    .alt(alt)
                    .qual(quality)
                    .firstOfPair(firstOfPairFlag)
                    .build();
            qualityMap.computeIfAbsent(key, QualityCount::new).increment();
        }
    }

}
