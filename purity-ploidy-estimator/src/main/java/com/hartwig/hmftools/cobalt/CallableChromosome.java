package com.hartwig.hmftools.cobalt;

import java.io.File;
import java.util.concurrent.Callable;

import com.hartwig.hmftools.purple.LoadSomaticVariants;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

class CallableChromosome implements Callable<ChromosomeCount> {

    private static final Logger LOGGER = LogManager.getLogger(LoadSomaticVariants.class);

    private final File inputFile;
    private final String chromosome;
    private final SamReaderFactory readerFactory;
    private final ChromosomeCount chromosomeCount;

    CallableChromosome(final String chromosome, final long length, final int windowSize, final SamReaderFactory readerFactory,
            final File inputFile) {
        this.chromosomeCount = new ChromosomeCount(chromosome, length, windowSize);
        this.chromosome = chromosome;
        this.readerFactory = readerFactory;
        this.inputFile = inputFile;
    }

    @Override
    public ChromosomeCount call() throws Exception {

        LOGGER.info("Generating windows on chromosome {}", chromosome);
        SamReader reader = readerFactory.open(inputFile);
        final SAMRecordIterator iterator = reader.query(chromosome, 0, 0, true);
        while (iterator.hasNext()) {
            chromosomeCount.addRecord(iterator.next());
        }

        return chromosomeCount;
    }

}
