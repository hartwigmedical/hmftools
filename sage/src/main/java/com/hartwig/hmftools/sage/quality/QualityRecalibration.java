package com.hartwig.hmftools.sage.quality;

import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutorService;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class QualityRecalibration {

    private final ExecutorService executorService;
    private final IndexedFastaSequenceFile refGenome;
    private final List<CompletableFuture<Collection<QualityCounter>>> counters = Lists.newArrayList();

    public QualityRecalibration(final ExecutorService executorService, final IndexedFastaSequenceFile refGenome) {
        this.executorService = executorService;
        this.refGenome = refGenome;
    }

    @NotNull
    public CompletableFuture<List<QualityRecalibrationRecord>> qualityRecalibrationRecords(@NotNull final String bamFile) {

        for (final SAMSequenceRecord sequenceRecord : refGenome.getSequenceDictionary().getSequences()) {
            if (HumanChromosome.contains(sequenceRecord.getSequenceName())) {
                final HumanChromosome chromosome = HumanChromosome.fromString(sequenceRecord.getSequenceName());
                if (chromosome.isAutosome()) {
                    int start = sequenceRecord.getSequenceLength() - 3_000_000;
                    int end = sequenceRecord.getSequenceLength() - 1_000_001;
                    submitAllRegions(bamFile, sequenceRecord.getSequenceName(), start, end);
                }
            }
        }

        final Map<QualityCounterKey, QualityCounter> map = Maps.newHashMap();
        CompletableFuture<Void> done = CompletableFuture.completedFuture(null);

        final Iterator<CompletableFuture<Collection<QualityCounter>>> counterIterator = counters.iterator();
        // Merge each region as it comes in
        while (counterIterator.hasNext()) {
            CompletableFuture<Collection<QualityCounter>> region = counterIterator.next();
            done = done.thenCombine(region, (aVoid, qualityCounters) -> {

                for (QualityCounter count : qualityCounters) {
                    final QualityCounterKey key = withoutPosition(count);
                    map.computeIfAbsent(key, QualityCounter::new).increment(count.count());
                }

                return null;
            });

            counterIterator.remove();
        }

        return done.thenApply(aVoid -> {
            final List<QualityCounter> sortedList = Lists.newArrayList(map.values());
            Collections.sort(sortedList);
            return QualityRecalibrationFactory.create(sortedList);
        });
    }

    public void addRegion(String bam, String contig, int start, int end) {
        final GenomeRegion bounds = GenomeRegions.create(contig, start, end);
        counters.add(CompletableFuture.supplyAsync(() -> new QualityCounterFactory(bam, refGenome).regionCount(bounds), executorService));
    }

    public void submitAllRegions(@NotNull final String bam, @NotNull final String contig, int minPosition, int maxPosition) {
        final int regionSliceSize = 100_000;
        for (int i = 0; ; i++) {
            int start = minPosition + i * regionSliceSize;
            int end = start + regionSliceSize - 1;

            if (end < minPosition) {
                continue;
            }

            addRegion(bam, contig, Math.max(start, minPosition), Math.min(end, maxPosition));

            if (end >= maxPosition) {
                break;
            }
        }
    }

    @NotNull
    private static QualityCounterKey withoutPosition(@NotNull final QualityCounter count) {
        return ImmutableQualityCounterKey.builder().from(count).position(0).build();
    }

}
