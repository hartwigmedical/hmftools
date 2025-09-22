package com.hartwig.hmftools.cobalt.count;

import static com.hartwig.hmftools.cobalt.CobaltConfig.CB_LOGGER;
import static com.hartwig.hmftools.cobalt.CobaltConstants.PARTITION_SIZE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.cobalt.ChromosomeData;
import com.hartwig.hmftools.cobalt.ChromosomePositionCodec;
import com.hartwig.hmftools.cobalt.CobaltColumns;
import com.hartwig.hmftools.cobalt.CobaltConfig;
import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.commons.lang3.Validate;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import tech.tablesaw.api.DoubleColumn;
import tech.tablesaw.api.IntColumn;
import tech.tablesaw.api.Row;
import tech.tablesaw.api.StringColumn;
import tech.tablesaw.api.Table;

public class BRC
{
    private final CobaltConfig mConfig;
    private final String mBamPath;
    private final SamReaderFactory mReaderFactory;
    private final List<Future<?>> tasks = new ArrayList<>();
    private final List<ChromosomeData> mChromosomes = Lists.newArrayList();
    private final ReadDepthAccumulator mReadDepthAccumulator;
    private final ChromosomePositionCodec mChromosomePosCodec;

    public BRC(
            final int windowSize, final CobaltConfig config,
            final ExecutorService executorService,
            final String bamPath,
            final ChromosomePositionCodec chromosomePosCodec) throws IOException
    {
        mConfig = config;
        mBamPath = bamPath;
        mReaderFactory = config.readerFactory();
        mChromosomePosCodec = chromosomePosCodec;
        loadChromosomes();
        mReadDepthAccumulator = new ReadDepthAccumulator(windowSize);
        for(ChromosomeData chromosome : mChromosomes)
        {
            mReadDepthAccumulator.addChromosome(chromosome.Name, chromosome.Length);
        }

        CB_LOGGER.info("calculating read depths from {}", mBamPath);
        for(ChrBaseRegion baseRegion : partitionGenome())
        {
            Runnable task = () -> sliceRegionTask(baseRegion);
            tasks.add(executorService.submit(task));
        }
    }

    public List<ChromosomeData> chromosomes()
    {
        return mChromosomes;
    }

    public Table generateDepths() throws ExecutionException, InterruptedException
    {
        // wait for all tasks to complete
        for(Future<?> f : tasks)
        {
            f.get();
        }
        CB_LOGGER.info("read depth complete");
        return calcDepths();
    }

    private void sliceRegionTask(ChrBaseRegion region)
    {
        CB_LOGGER.debug("region({}) accumulating read depth", region);
        final File bamFile = new File(mBamPath);
        try(SamReader reader = mReaderFactory.open(bamFile))
        {
            BamSlicer bamSlicer = new BamSlicer(mConfig.MinMappingQuality, mConfig.IncludeDuplicates, false, false);
            bamSlicer.slice(reader, region, samRecord -> processRead(samRecord, region));
        }
        catch(IOException e)
        {
            CB_LOGGER.warn("bam reading failed", e);
        }
        CB_LOGGER.debug("region({}) complete", region);
    }

    private void processRead(final SAMRecord record, ChrBaseRegion region)
    {
        if(mConfig.IncludeDuplicates)
        {
            // revert to only analysing the raw reads
            if(record.hasAttribute(CONSENSUS_READ_ATTRIBUTE))
            {
                return;
            }
        }
        else
        {
            if(record.getDuplicateReadFlag())
            {
                return;
            }
        }

        Validate.isTrue(record.getContig().equals(region.Chromosome));

        for(AlignmentBlock currentBlock : record.getAlignmentBlocks())
        {
            accumulateAlignmentBlock(region, currentBlock, record.getReadBases());
        }
    }

    void accumulateAlignmentBlock(ChrBaseRegion region, AlignmentBlock block, byte[] readBases)
    {
        int alignmentBlockReadStart = block.getReadStart();
        int alignmentBlockReferenceStart = block.getReferenceStart();
        int alignmentBlockLength = block.getLength();
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
        mReadDepthAccumulator.addReadAlignmentToCounts(region.Chromosome, genomeStart, length, readBases, readStartIndex);
    }

    public ListMultimap<Chromosome, DepthReading> calculateReadDepths()throws ExecutionException, InterruptedException
    {
        for(Future<?> f : tasks)
        {
            f.get();
        }

        ListMultimap<Chromosome, DepthReading> result = ArrayListMultimap.create();
        for(ChromosomeData chromosome : mChromosomes)
        {
            List<DepthReading> readDepths = mReadDepthAccumulator.getChromosomeReadDepths(chromosome.Name);
            Objects.requireNonNull(readDepths);
            HumanChromosome humanChromosome = HumanChromosome.fromString(chromosome.Name);
            result.putAll(humanChromosome, readDepths);
        }
        return result;
    }
    private Table calcDepths()
    {
        final Table readDepthTable = Table.create("readDepths",
                StringColumn.create(CobaltColumns.CHROMOSOME),
                IntColumn.create(CobaltColumns.POSITION),
                DoubleColumn.create(CobaltColumns.READ_DEPTH),
                DoubleColumn.create(CobaltColumns.READ_GC_CONTENT));

        for(ChromosomeData chromosome : mChromosomes)
        {
            List<DepthReading> readDepths = mReadDepthAccumulator.getChromosomeReadDepths(chromosome.Name);
            Objects.requireNonNull(readDepths);
            for(DepthReading readDepth : readDepths)
            {
                Row row = readDepthTable.appendRow();
                row.setString(CobaltColumns.CHROMOSOME, chromosome.Name);
                row.setInt(CobaltColumns.POSITION, readDepth.StartPosition);
                row.setDouble(CobaltColumns.READ_DEPTH, readDepth.ReadDepth);
                row.setDouble(CobaltColumns.READ_GC_CONTENT, readDepth.ReadGcContent);
            }
        }
        mChromosomePosCodec.addEncodedChrPosColumn(readDepthTable, false);

        return readDepthTable;
    }

    private void loadChromosomes() throws IOException
    {
        try(SamReader reader = mReaderFactory.open(new File(mBamPath)))
        {
            SAMSequenceDictionary dictionary = reader.getFileHeader().getSequenceDictionary();
            for(final SAMSequenceRecord samSequenceRecord : dictionary.getSequences())
            {
                String sequenceName = samSequenceRecord.getSequenceName();
                if(!HumanChromosome.contains(sequenceName))
                {
                    continue;
                }
                if(mConfig.SpecificChrRegions.hasFilters() && mConfig.SpecificChrRegions.excludeChromosome(sequenceName))
                {
                    continue;
                }
                mChromosomes.add(new ChromosomeData(sequenceName, samSequenceRecord.getSequenceLength()));
            }
        }
    }

    private List<ChrBaseRegion> partitionGenome()
    {
        List<ChrBaseRegion> partitions = new ArrayList<>();
        for(ChromosomeData chromosome : mChromosomes)
        {
            for(int startPos = 1; startPos < chromosome.Length; startPos += PARTITION_SIZE)
            {
                int endPos = Math.min(startPos + PARTITION_SIZE - 1, chromosome.Length);
                partitions.add(new ChrBaseRegion(chromosome.Name, startPos, endPos));
            }
        }
        return partitions;
    }
}
