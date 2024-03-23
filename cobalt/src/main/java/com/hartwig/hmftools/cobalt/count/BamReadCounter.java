package com.hartwig.hmftools.cobalt.count;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.CobaltConstants.PARTITION_SIZE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Objects;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import com.hartwig.hmftools.cobalt.Chromosome;
import com.hartwig.hmftools.cobalt.ChromosomePositionCodec;
import com.hartwig.hmftools.cobalt.CobaltColumns;
import com.hartwig.hmftools.cobalt.CobaltConfig;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.bam.BamSlicer;

import org.apache.commons.lang3.Validate;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import tech.tablesaw.api.*;

public class BamReadCounter
{
    private final int mMinMappingQuality;
    private final boolean mIncludeDuplicates;

    private Table mReferenceDepths = null;
    private Table mTumorDepths = null;

    private final ExecutorService mExecutorService;
    private final SamReaderFactory mReaderFactory;

    private Collection<Chromosome> mChromosomes = null;

    private final ReadDepthAccumulator mRefReadDepthAccumulator;
    private final ReadDepthAccumulator mTumorReadDepthAccumulator;

    private final ChromosomePositionCodec mChromosomePosCodec;

    public Table getReferenceDepths() { return mReferenceDepths; }
    public Table getTumorDepths() { return mTumorDepths; }

    public BamReadCounter(
            final int windowSize, final CobaltConfig config,
            final ExecutorService executorService, final SamReaderFactory readerFactory,
            final ChromosomePositionCodec chromosomePosCodec)
    {
        mMinMappingQuality = config.MinMappingQuality;
        mIncludeDuplicates = config.IncludeDuplicates;
        mExecutorService = executorService;
        mReaderFactory = readerFactory;
        mChromosomePosCodec = chromosomePosCodec;
        mRefReadDepthAccumulator = new ReadDepthAccumulator(windowSize);
        mTumorReadDepthAccumulator = new ReadDepthAccumulator(windowSize);
    }

    public void generateDepths(
            @Nullable final String referenceBam, @Nullable final String tumorBam)
            throws ExecutionException, InterruptedException, IOException
    {
        if(referenceBam == null && tumorBam == null)
        {
            CB_LOGGER.error("no bam file supplied");
            return;
        }

        mChromosomes = loadChromosomes(mReaderFactory, referenceBam, tumorBam);

        List<Future<?>> tasks = new ArrayList<>();
        List<SamReader> samReaders = Collections.synchronizedList(new ArrayList<>());

        if(tumorBam != null)
        {
            CB_LOGGER.info("calculating read depths from {}", tumorBam);
            tasks.addAll(createFutures(tumorBam, mTumorReadDepthAccumulator, samReaders));
        }

        if(referenceBam != null)
        {
            CB_LOGGER.info("calculating read depths from {}", referenceBam);
            tasks.addAll(createFutures(referenceBam, mRefReadDepthAccumulator, samReaders));
        }

        // wait for all tasks to complete
        for(Future<?> f : tasks)
        {
            f.get();
        }

        if(tumorBam != null)
        {
            mTumorDepths = generateDepths(mTumorReadDepthAccumulator);
        }

        if(referenceBam != null)
        {
            mReferenceDepths = generateDepths(mRefReadDepthAccumulator);
        }

        // close all sam readers
        for(SamReader samReader : samReaders)
        {
            samReader.close();
        }

        CB_LOGGER.info("read Depth Complete");
    }

    private List<Future<?>> createFutures(final String bamFilePath, final ReadDepthAccumulator readDepthCounter, List<SamReader> samReaderList)
    {
        // add all the chromosomes
        for(Chromosome chromosome : mChromosomes)
        {
            readDepthCounter.addChromosome(chromosome.contig, chromosome.length);
        }

        final File bamFile = new File(bamFilePath);

        // One bam reader per thread instead of one per task
        final ThreadLocal<SamReader> threadBamReader = ThreadLocal.withInitial(() -> {
            SamReader samReader = mReaderFactory.open(bamFile);
            samReaderList.add(samReader);
            return samReader;
        });

        final List<ChrBaseRegion> partitions = partitionGenome();
        final List<Future<?>> futures = new ArrayList<>();
        for(ChrBaseRegion baseRegion : partitions)
        {
            Runnable task = () -> sliceRegionTask(threadBamReader, baseRegion, readDepthCounter);
            futures.add(mExecutorService.submit(task));
        }

        return futures;
    }

    private void sliceRegionTask(ThreadLocal<SamReader> samReaderSupplier, ChrBaseRegion region, ReadDepthAccumulator readDepthAccumulator)
    {
        CB_LOGGER.debug("region({}) accumulating read depth", region);

        final SamReader reader = samReaderSupplier.get();

        BamSlicer bamSlicer = new BamSlicer(mMinMappingQuality, mIncludeDuplicates, false, false);

        bamSlicer.slice(reader, region, samRecord -> processRead(samRecord, region, readDepthAccumulator));

        CB_LOGGER.debug("region({}) complete", region);
    }

    private void processRead(final SAMRecord record, ChrBaseRegion region, ReadDepthAccumulator readDepthAccumulator)
    {
        if(mIncludeDuplicates)
        {
            // revert to only analysing the raw reads
            if(record.hasAttribute(CONSENSUS_READ_ATTRIBUTE))
                return;
        }
        else
        {
            if(record.getDuplicateReadFlag())
                return;
        }

        Validate.isTrue(record.getContig().equals(region.Chromosome));

        for(AlignmentBlock currentBlock : record.getAlignmentBlocks())
        {
            accumulateAlignmentBlock(region, readDepthAccumulator, currentBlock.getReadStart(), currentBlock.getReferenceStart(),
                    currentBlock.getLength(), record.getReadBases());
        }
    }

    static void accumulateAlignmentBlock(
            final ChrBaseRegion region, final ReadDepthAccumulator readDepthAccumulator,
            final int alignmentBlockReadStart, final int alignmentBlockReferenceStart, final int alignmentBlockLength,
            byte[] readBases)
    {
        // NOTE: we need to adjust start and end to avoid adding counts to regions that belongs to another task
        // as if we do that they will be double counted
        int genomeStart = Math.max(alignmentBlockReferenceStart, region.start());
        int length = Math.min(alignmentBlockReferenceStart + alignmentBlockLength, region.end() + 1) - genomeStart;

        if(length <= 0)
        {
            return;
        }

        // use 0 based index here such that we can use it with java string
        int readStartIndex = alignmentBlockReadStart - 1;
        readStartIndex += (genomeStart - alignmentBlockReferenceStart);
        readDepthAccumulator.addReadAlignmentToCounts(region.Chromosome, genomeStart, length, readBases, readStartIndex);
    }

    private Table generateDepths(ReadDepthAccumulator readDepthAccumulator)
    {
        final Table readDepthTable = Table.create("readDepths",
                StringColumn.create(CobaltColumns.CHROMOSOME),
                IntColumn.create(CobaltColumns.POSITION),
                DoubleColumn.create(CobaltColumns.READ_DEPTH),
                DoubleColumn.create(CobaltColumns.READ_GC_CONTENT));

        for (Chromosome chromosome : mChromosomes)
        {
            List<ReadDepth> readDepths = readDepthAccumulator.getChromosomeReadDepths(chromosome.contig);
            Objects.requireNonNull(readDepths);
            for (ReadDepth readDepth : readDepths)
            {
                Row row = readDepthTable.appendRow();
                row.setString(CobaltColumns.CHROMOSOME, chromosome.contig);
                row.setInt(CobaltColumns.POSITION, readDepth.StartPosition);
                row.setDouble(CobaltColumns.READ_DEPTH, readDepth.ReadDepth);
                row.setDouble(CobaltColumns.READ_GC_CONTENT, readDepth.ReadGcContent);
            }
        }

        mChromosomePosCodec.addEncodedChrPosColumn(readDepthTable, false);

        return readDepthTable;
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

    private List<ChrBaseRegion> partitionGenome()
    {
        List<ChrBaseRegion> partitions = new ArrayList<>();
        for(Chromosome chromosome : mChromosomes)
        {
            /*if(mConfig.SpecificChrRegions.excludeChromosome(chromosomeStr))
                continue; */

            for(int startPos = 1; startPos < chromosome.length; startPos += PARTITION_SIZE)
            {
                int endPos = Math.min(startPos + PARTITION_SIZE - 1, chromosome.length);
                partitions.add(new ChrBaseRegion(chromosome.contig, startPos, endPos));
            }
        }
        return partitions;
    }
}
