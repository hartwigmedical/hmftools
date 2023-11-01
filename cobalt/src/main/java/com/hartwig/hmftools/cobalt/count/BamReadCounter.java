package com.hartwig.hmftools.cobalt.count;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import com.hartwig.hmftools.cobalt.Chromosome;
import com.hartwig.hmftools.cobalt.ChromosomePositionCodec;
import com.hartwig.hmftools.cobalt.CobaltColumns;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.apache.commons.lang3.Validate;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import tech.tablesaw.api.*;

public class BamReadCounter
{
    private final int mWindowSize;
    private final int mMinMappingQuality;

    Table mReferenceCounts = null;
    Table mTumorCounts = null;

    private final ExecutorService mExecutorService;
    private final SamReaderFactory mReaderFactory;
    private Collection<Chromosome> mChromosomes = null;

    private final ChromosomePositionCodec mChromosomePosCodec;

    public Table getReferenceCounts() { return mReferenceCounts; }
    public Table getTumorCounts() { return mTumorCounts; }

    public BamReadCounter(
            final int windowSize, final int minMappingQuality,
            final ExecutorService executorService, final SamReaderFactory readerFactory,
            ChromosomePositionCodec chromosomePosCodec)
    {
        mWindowSize = windowSize;
        mMinMappingQuality = minMappingQuality;
        mExecutorService = executorService;
        mReaderFactory = readerFactory;
        mChromosomePosCodec = chromosomePosCodec;
    }

    public void generateCounts(
            @Nullable final String referenceBam, @Nullable final String tumorBam)
            throws ExecutionException, InterruptedException, IOException
    {
        if (referenceBam == null && tumorBam == null)
        {
            CB_LOGGER.error("No bam file supplied");
            return;
        }

        mChromosomes = loadChromosomes(mReaderFactory, referenceBam, tumorBam);

        List<Future<ChromosomeReadCount>> tumorFutures = null;
        List<Future<ChromosomeReadCount>> referenceFutures = null;

        if (tumorBam != null)
        {
            CB_LOGGER.info("Calculating Read Count from {}", tumorBam);
            tumorFutures = createFutures(mReaderFactory, new File(tumorBam));
        }

        if (referenceBam != null)
        {
            CB_LOGGER.info("Calculating Read Count from {}", referenceBam);
            referenceFutures = createFutures(mReaderFactory, new File(referenceBam));
        }

        if (tumorFutures != null)
        {
            mTumorCounts = fromFutures(tumorFutures);
        }
        if (referenceFutures != null)
        {
            mReferenceCounts = fromFutures(referenceFutures);
        }

        CB_LOGGER.info("Read Count Complete");
    }

    private List<Future<ChromosomeReadCount>> createFutures(
            final SamReaderFactory readerFactory, final File file)
    {
        final List<Future<ChromosomeReadCount>> futures = new ArrayList<>();
        for(Chromosome chromosome : mChromosomes)
        {
            final ChromosomeReadCount callable = new ChromosomeReadCount(
                    file, readerFactory, chromosome, mWindowSize, mMinMappingQuality);

            futures.add(mExecutorService.submit(callable));
        }

        return futures;
    }

    private Table fromFutures(List<Future<ChromosomeReadCount>> futures)
            throws ExecutionException, InterruptedException
    {
        final Table readCounts = Table.create("readCounts",
                StringColumn.create(CobaltColumns.CHROMOSOME),
                IntColumn.create(CobaltColumns.POSITION),
                IntColumn.create(CobaltColumns.READ_COUNT));

        for (Future<ChromosomeReadCount> future : futures)
        {
            final ChromosomeReadCount readCount = future.get();
            final Chromosome chromosome = readCount.chromosome();

            for (ReadCount rc : readCount.readCount())
            {
                Row row = readCounts.appendRow();
                row.setString(CobaltColumns.CHROMOSOME, chromosome.contig);
                row.setInt(CobaltColumns.POSITION, rc.position());
                row.setInt(CobaltColumns.READ_COUNT, rc.readCount());
            }
        }

        mChromosomePosCodec.addEncodedChrPosColumn(readCounts, false);

        return readCounts;
    }

    private Collection<Chromosome> loadChromosomes(final SamReaderFactory readerFactory,
            @Nullable final String referenceBam,
            @Nullable final String tumorBam) throws IOException
    {
        Collection<Chromosome> chromosomes = new ArrayList<>();

        Validate.isTrue(referenceBam != null || tumorBam != null);

        try (SamReader reader = readerFactory.open(new File(referenceBam != null ? referenceBam : tumorBam)))
        {
            SAMSequenceDictionary dictionary = reader.getFileHeader().getSequenceDictionary();

            for(final SAMSequenceRecord samSequenceRecord : dictionary.getSequences())
            {
                String sequenceName = samSequenceRecord.getSequenceName();

                if(HumanChromosome.contains(sequenceName))
                {
                    chromosomes.add(new Chromosome(sequenceName, samSequenceRecord.getSequenceLength()));
                }
            }
        }

        return chromosomes;
    }
}
