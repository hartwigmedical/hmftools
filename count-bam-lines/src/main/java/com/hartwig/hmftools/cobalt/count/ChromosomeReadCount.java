package com.hartwig.hmftools.cobalt.count;

import java.io.File;
import java.util.List;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cobalt.ImmutableReadCount;
import com.hartwig.hmftools.common.cobalt.ReadCount;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.window.Window;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

class ChromosomeReadCount implements Callable<ChromosomeReadCount> {

    private static final Logger LOGGER = LogManager.getLogger(ChromosomeReadCount.class);

    private final File inputFile;
    private final SamReaderFactory readerFactory;
    private final String chromosome;
    private final long chromosomeLength;
    private final List<ReadCount> result = Lists.newArrayList();
    private final int minMappingQuality;
    private final Window window;

    private long start;
    private int count;

    ChromosomeReadCount(final File inputFile, final SamReaderFactory readerFactory, @NotNull final String chromosome,
            final long chromosomeLength, final int windowSize, final int minMappingQuality) {
        this.inputFile = inputFile;
        this.readerFactory = readerFactory;
        this.chromosome = chromosome;
        this.chromosomeLength = chromosomeLength;
        this.minMappingQuality = minMappingQuality;
        this.window = new Window(windowSize);

        start = 1;
        count = -1;
    }

    @Override
    public ChromosomeReadCount call() throws Exception {
        LOGGER.info("Generating windows on chromosome {}", chromosome);

        try (final SamReader reader = readerFactory.open(inputFile)) {
            final SAMRecordIterator iterator = reader.query(chromosome, 0, 0, true);
            while (iterator.hasNext()) {
                addRecord(iterator.next());
            }
        }
        return this;
    }

    @NotNull
    Chromosome chromosome() {
        return HumanChromosome.fromString(chromosome);
    }

    @NotNull
    List<ReadCount> readCount() {
        addReadCount(start, count);

        long lastWindowPosition = lastWindowPosition();
        if (result.get(result.size() - 1).position() < lastWindowPosition) {
            addReadCount(lastWindowPosition, -1);
        }

        return result;
    }

    private void addRecord(@NotNull SAMRecord record) {
        if (isEligible(record)) {
            long window = windowPosition(record.getAlignmentStart());
            if (start != window) {
                addReadCount(start, count);
                start = window;
                count = 0;
            }
            count++;
        }
    }

    private void addReadCount(long position, int count) {
        result.add(ImmutableReadCount.builder().chromosome(chromosome).position(position).readCount(count).build());
    }

    private boolean isEligible(@NotNull SAMRecord record) {
        return record.getMappingQuality() >= minMappingQuality && !(record.getReadUnmappedFlag() || record.getDuplicateReadFlag()
                || record.isSecondaryOrSupplementary());
    }

    private long lastWindowPosition() {
        return windowPosition(chromosomeLength);
    }

    private long windowPosition(long position) {
        return window.start(position);
    }
}
