package com.hartwig.hmftools.sage.quality;

import java.util.Collection;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class QualityRecalibration {

    private final ExecutorService executorService;
    private final IndexedFastaSequenceFile refGenome;
    private final List<Future<Collection<QualityCounter>>> counters = Lists.newArrayList();

    public QualityRecalibration(final ExecutorService executorService, final IndexedFastaSequenceFile refGenome) {
        this.executorService = executorService;
        this.refGenome = refGenome;
    }

    @NotNull
    public List<QualityRecalibrationRecord> qualityRecalibrationRecords(@NotNull final String bamFile)
            throws ExecutionException, InterruptedException {

        for (final SAMSequenceRecord sequenceRecord : refGenome.getSequenceDictionary().getSequences()) {
            if (HumanChromosome.contains(sequenceRecord.getSequenceName())) {
                final HumanChromosome chromosome = HumanChromosome.fromString(sequenceRecord.getSequenceName());
                if (chromosome.isAutosome()) {
                    int start = sequenceRecord.getSequenceLength() - 3_000_000;
                    int end = sequenceRecord.getSequenceLength() - 1_000_001;
                    addAllRegions(bamFile, sequenceRecord.getSequenceName(), start, end);
                }
            }
        }

        List<QualityCounter> allCounts = Lists.newArrayList();
        for (Future<Collection<QualityCounter>> counter : counters) {
            allCounts.addAll(counter.get());
            allCounts = QualityCounterGrouping.groupWithoutPosition(allCounts);
        }

        return QualityRecalibrationFactory.create(allCounts);
    }

    public void addRegion(String bam, String contig, int start, int end) {
        final GenomeRegion bounds = GenomeRegions.create(contig, start, end);
        final Future<Collection<QualityCounter>> future =
                executorService.submit(() -> new QualityCounterFactory(bam, refGenome).regionCount(bounds));
        counters.add(future);
    }

    public void addAllRegions(String bam, String contig, int minPosition, int maxPosition) {
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
}
