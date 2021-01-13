package com.hartwig.hmftools.sage.quality;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.sage.config.SageConfig;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class QualityRecalibration {

    private final ExecutorService executorService;
    private final IndexedFastaSequenceFile refGenome;
    private final SageConfig config;

    public QualityRecalibration(final SageConfig config, final ExecutorService executorService,
            final IndexedFastaSequenceFile refGenome) {
        this.executorService = executorService;
        this.refGenome = refGenome;
        this.config = config;
    }

    @NotNull
    public CompletableFuture<List<QualityRecalibrationRecord>> qualityRecalibrationRecords(@NotNull final String bamFile) {
        final Map<QualityCounterKey, QualityCounter> map = new ConcurrentHashMap<>(Maps.newHashMap());
        final List<CompletableFuture<Void>> doneList = Lists.newArrayList();

        for (final SAMSequenceRecord sequenceRecord : refGenome.getSequenceDictionary().getSequences()) {
            final String contig = sequenceRecord.getSequenceName();

            if (HumanChromosome.contains(contig) && HumanChromosome.fromString(contig).isAutosome()) {
                int start = sequenceRecord.getSequenceLength() - 1_000_000 - config.baseQualityRecalibrationConfig().sampleSize();
                int end = sequenceRecord.getSequenceLength() - 1_000_001;
                for (CompletableFuture<Collection<QualityCounter>> region : submitAllRegions(bamFile, contig, start, end)) {
                    final CompletableFuture<Void> done = region.thenAccept(counts -> {
                        for (QualityCounter count : counts) {
                            final QualityCounterKey key = withoutPosition(count);
                            map.computeIfAbsent(key, QualityCounter::new).increment(count.count());
                        }
                    });
                    doneList.add(done);
                }
            }
        }

        return CompletableFuture.allOf(doneList.toArray(new CompletableFuture[0])).thenApply(aVoid -> {
            final List<QualityCounter> sortedList = Lists.newArrayList(map.values());
            Collections.sort(sortedList);
            return QualityRecalibrationFactory.create(sortedList);
        });
    }

    public CompletableFuture<Collection<QualityCounter>> addRegion(String bam, String contig, int start, int end) {
        final GenomeRegion bounds = GenomeRegions.create(contig, start, end);
        return CompletableFuture.supplyAsync(() -> new QualityCounterFactory(config, bam, refGenome).regionCount(bounds),
                executorService);
    }

    public List<CompletableFuture<Collection<QualityCounter>>> submitAllRegions(@NotNull final String bam, @NotNull final String contig,
            int minPosition, int maxPosition) {
        final List<CompletableFuture<Collection<QualityCounter>>> result = Lists.newArrayList();

        final int regionSliceSize = 100_000;
        for (int i = 0; ; i++) {
            int start = minPosition + i * regionSliceSize;
            int end = start + regionSliceSize - 1;

            if (end < minPosition) {
                continue;
            }

            result.add(addRegion(bam, contig, Math.max(start, minPosition), Math.min(end, maxPosition)));

            if (end >= maxPosition) {
                break;
            }
        }

        return result;
    }

    @NotNull
    private static QualityCounterKey withoutPosition(@NotNull final QualityCounter count) {
        return ImmutableQualityCounterKey.builder().from(count).position(0).build();
    }
}
