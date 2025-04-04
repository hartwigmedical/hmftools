package com.hartwig.hmftools.bamtools.checker;

import static java.lang.Math.ceil;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.region.PartitionUtils.buildPartitions;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BAM_EXTENSION;

import java.io.File;
import java.util.List;
import java.util.NoSuchElementException;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.bamtools.common.PartitionTask;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.TaskQueue;

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

        if(config.writeBam())
        {
            // create a BAM writer per thread
            mBamFilename = config.formFilename(format("%s_%02d", UNSORTED_BAM_ID, threadId), BAM_EXTENSION);

            SAMFileHeader fileHeader = mSamReader.getFileHeader().clone();
            fileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

            mBamWriter = new SAMFileWriterFactory().makeBAMWriter(fileHeader, false, new File(mBamFilename));
        }
        else
        {
            mBamWriter = null;
            mBamFilename = "";
        }

        mPartitionChecker = new PartitionChecker(config, fragmentCache, mSamReader, mBamWriter);
    }

    public static List<ChrBaseRegion> splitRegionsIntoPartitions(final CheckConfig config)
    {
        return PartitionTask.splitRegionsIntoPartitions(
                config.BamFile, config.RefGenomeFile, config.Threads, config.SpecificChrRegions, config.PartitionSize);
    }

    public FragmentStats stats() { return mPartitionChecker.stats(); }
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
        if(mBamWriter != null)
            mBamWriter.close();
    }
}
