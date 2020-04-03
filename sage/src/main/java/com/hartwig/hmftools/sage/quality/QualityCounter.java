package com.hartwig.hmftools.sage.quality;

import java.util.Collection;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

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

        final Set<QualityRecord> altsToRemove = QualityGrouping.groupByAlt(qualityMap.values())
                .stream()
                .filter(x -> x.ref() != x.alt())
                .filter(x -> x.count() > 3)
                .map(QualityCount::key)
                .collect(Collectors.toSet());

        final Set<QualityCount> result = Sets.newHashSet();
        for (QualityCount count : qualityMap.values()) {
            final QualityRecord altKey = QualityGrouping.alt(count);
            if (!indelPositions.contains(count.position()) && !altsToRemove.contains(altKey)) {
                result.add(count);
            }
        }

        return QualityGrouping.groupWithoutPosition(result);
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
    public void handleAlignment(@NotNull final SAMRecord r, @NotNull final CigarElement e, final int startReadIndex, final int refPos) {
        for (int i = 0; i < e.getLength(); i++) {
            int readIndex = startReadIndex + i;
            int position = refPos + i;

            if (position > bounds.end()) {
                return;
            }

            if (position < bounds.start()) {
                continue;
            }

            byte ref = refGenome.base(position);
            byte alt = r.getReadBases()[readIndex];
            byte quality = r.getBaseQualities()[readIndex];
            boolean firstOfPairFlag = r.getFirstOfPairFlag();

            final QualityRecord key = ImmutableQualityRecord.builder()
                    .ref(ref)
                    .alt(alt)
                    .qual(quality)
                    .position(position)
                    .firstOfPair(firstOfPairFlag)
                    .trinucleotideContext(refGenome.trinucleotideContext(position))
                    .build();
            qualityMap.computeIfAbsent(key, QualityCount::new).increment();
        }
    }

}
