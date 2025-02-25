package com.hartwig.hmftools.bamtools.checker;

import static java.lang.Math.ceil;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome.formHumanChromosomeRegions;
import static com.hartwig.hmftools.common.region.PartitionUtils.buildPartitions;
import static com.hartwig.hmftools.common.region.PartitionUtils.partitionChromosome;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BAM_EXTENSION;

import java.io.File;
import java.util.List;
import java.util.NoSuchElementException;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.utils.TaskQueue;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class PartitionThread extends Thread
{
    private final SAMFileWriter mBamWriter;
    private final SamReader mSamReader;
    private final TaskQueue mPartitions;

    private final PartitionChecker mPartitionChecker;
    private final String mBamFilename;

    protected static final String UNSORTED_BAM_ID = "unsorted";
    protected static final String SORTED_BAM_ID = "sorted";

    public PartitionThread(
            final CheckConfig config, final FragmentCache fragmentCache, final TaskQueue partitions, final int threadId)
    {
        mPartitions = partitions;

        mSamReader = SamReaderFactory.makeDefault()
                .referenceSequence(new File(config.RefGenomeFile))
                .open(new File(config.BamFile));

        // create a BAM writer per thread
        mBamFilename = config.formFilename(format("%s_%02d", UNSORTED_BAM_ID, threadId), BAM_EXTENSION);

        SAMFileHeader fileHeader = mSamReader.getFileHeader().clone();
        fileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

        mBamWriter = new SAMFileWriterFactory().makeBAMWriter(fileHeader, false, new File(mBamFilename));

        mPartitionChecker = new PartitionChecker(config, fragmentCache, mSamReader, mBamWriter);
    }

    // public PartitionChecker partitionChecker() { return mPartitionChecker; }
    public String bamFilename() { return mBamFilename; }

    public void run()
    {
        while(true)
        {
            try
            {
                ChrBaseRegion region = (ChrBaseRegion)mPartitions.removeItem();

                mPartitionChecker.processPartition(region);
            }
            catch(NoSuchElementException e)
            {
                BT_LOGGER.trace("all tasks complete");
                break;
            }
            catch(Exception e)
            {
                e.printStackTrace();
                System.exit(1);
            }
        }
    }

    public void writeUnmappedReads()
    {
        SAMRecordIterator iterator = mSamReader.queryUnmapped();

        while(iterator.hasNext())
        {
            SAMRecord record = iterator.next();
            mBamWriter.addAlignment(record);
        }
    }

    public void writeIncompleteReads(final List<SAMRecord> reads)
    {
        reads.forEach(x -> mBamWriter.addAlignment(x));
    }

    public void close()
    {
        mBamWriter.close();
    }

    public static List<ChrBaseRegion> splitRegionsIntoPartitions(final CheckConfig config)
    {
        List<ChrBaseRegion> partitionRegions = Lists.newArrayList();

        if(!config.SpecificChrRegions.Regions.isEmpty())
        {
            // only split by thread count if can be done simply
            if(config.SpecificChrRegions.Regions.size() == 1)
            {
                if(config.Threads > 1)
                {
                    ChrBaseRegion specificRegion = config.SpecificChrRegions.Regions.get(0);
                    int intervalLength = (int)ceil(specificRegion.baseLength() / (double)config.Threads);
                    int regionStart = specificRegion.start();

                    for(int i = 0 ; i < config.Threads; ++i)
                    {
                        int regionEnd = min(regionStart + intervalLength - 1, specificRegion.end());
                        partitionRegions.add(new ChrBaseRegion(specificRegion.Chromosome, regionStart, regionEnd));
                        regionStart = regionEnd + 1;
                    }
                }
                else
                {
                    partitionRegions.add(config.SpecificChrRegions.Regions.get(0));
                }
            }
            else
            {
                for(int i = 0 ; i < config.Threads; ++i)
                {
                    partitionRegions.add(config.SpecificChrRegions.Regions.get(i));
                }
            }

            return partitionRegions;
        }

        SamReader samReader = SamReaderFactory.makeDefault()
                .referenceSequence(new File(config.RefGenomeFile))
                .open(new File(config.BamFile));

        SAMFileHeader fileHeader = samReader.getFileHeader().clone();

        for(final SAMSequenceRecord sequenceRecord : fileHeader.getSequenceDictionary().getSequences())
        {
            String chromosome = sequenceRecord.getSequenceName();

            if(config.SpecificChrRegions != null && config.SpecificChrRegions.excludeChromosome(chromosome))
                continue;

            partitionRegions.addAll(buildPartitions(chromosome, sequenceRecord.getEnd(), config.PartitionSize));
        }

        return partitionRegions;
    }
}
