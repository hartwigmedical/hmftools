package com.hartwig.hmftools.redux.write;

import static com.hartwig.hmftools.common.bamops.BamMerger.buildCombinedHeader;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BAM_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BAM_INDEX_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.filenamePart;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.common.Constants.FILE_ID;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bamops.BamOperations;
import com.hartwig.hmftools.common.bamops.BamToolName;
import com.hartwig.hmftools.common.basequal.jitter.JitterAnalyser;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.redux.ReduxConfig;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;

public class FileWriterCache
{
    private final ReduxConfig mConfig;
    private final ReadDataWriter mReadDataWriter;

    private final List<BamWriter> mBamWriters;

    private final List<PartitionInfo> mPartitions;
    private final Queue<PartitionInfo> mCompletedPartitionsQueue;
    private int mCompletedPartitionCount;

    private final BamWriterSync mUnmappingWriter;
    private String mUnmappingSortedBamFilename;

    private String mFinalBamFilename;

    private final BamWriterSync mFullUnmappedWriter;

    private final JitterAnalyser mJitterAnalyser;

    private static final String SORTED_ID = "sorted";

    // specific BAM names
    private static final String UNMAPPING = "unmapping";
    private static final String UNMAPPING_SORTED = "unmapping_sorted";
    private static final String FULL_UNMAPPED = "full_unmapped";

    public FileWriterCache(final ReduxConfig config, @Nullable final JitterAnalyser jitterAnalyser)
    {
        mConfig = config;
        mJitterAnalyser = jitterAnalyser;

        mReadDataWriter = new ReadDataWriter(mConfig);

        mPartitions = Lists.newArrayList();
        mBamWriters = Lists.newArrayList();

        mCompletedPartitionsQueue = new LinkedBlockingQueue<>();
        mCompletedPartitionCount = 0;

        // create a shared BAM writer if either no multi-threading or using the sorted BAM writer
        if(!mConfig.WriteBam)
        {
            BamWriter bamWriter = new BamWriterNone("none", mConfig, mReadDataWriter, jitterAnalyser);
            mBamWriters.add(bamWriter);

            mUnmappingWriter = (BamWriterSync)bamWriter;
            mFullUnmappedWriter = (BamWriterSync)bamWriter;
            mUnmappingSortedBamFilename = "";
            mFinalBamFilename = "";
            return;
        }

        mFinalBamFilename = mConfig.OutputBam != null ? mConfig.OutputBam : formBamFilename(null, null);

        if(mConfig.UnmapRegions.enabled())
        {
            String unmappingFilename = formBamFilename(null, UNMAPPING);
            mUnmappingWriter = (BamWriterSync)createBamWriter(unmappingFilename, true);
            mUnmappingSortedBamFilename = formBamFilename(null, UNMAPPING_SORTED);

            String fullyUnmappedFilename = formBamFilename(null, FULL_UNMAPPED);
            mFullUnmappedWriter = (BamWriterSync)createBamWriter(fullyUnmappedFilename, true);
        }
        else
        {
            mUnmappingWriter = null;
            mFullUnmappedWriter = null;
            mUnmappingSortedBamFilename = "";
        }
    }

    public void addPartition(final List<ChrBaseRegion> regions)
    {
        int partitionIndex = mPartitions.size();

        String filename;

        if(partitionIndex == 0 && mConfig.ParallelConcatenation)
        {
            // use the final BAM name for the writer which writes the first partition's reads
            filename = mFinalBamFilename;
        }
        else
        {
            filename = formBamFilename(SORTED_ID, String.valueOf(partitionIndex));
        }

        BamWriter bamWriter = createBamWriter(filename, false);
        PartitionInfo partitionInfo = new PartitionInfo(partitionIndex, regions, bamWriter);
        mPartitions.add(partitionInfo);
    }

    public List<PartitionInfo> partitions() { return mPartitions; }
    public int partitionCount() { return mPartitions.size(); }
    public Queue<PartitionInfo> completedPartitionsQueue() { return mCompletedPartitionsQueue; }

    public synchronized void addCompletedPartition(final PartitionInfo partition)
    {
        ++mCompletedPartitionCount;
        RD_LOGGER.debug("completed {} partition readers", mCompletedPartitionCount);

        if(mConfig.ParallelConcatenation)
            mCompletedPartitionsQueue.add(partition);
        else
            partition.bamWriter().close();
    }

    public BamWriterSync getUnmappingBamWriter() { return mUnmappingWriter; }
    public BamWriterSync getFullUnmappedBamWriter() { return mFullUnmappedWriter; }

    public ReadDataWriter readDataWriter() { return mReadDataWriter; }
    public List<BamWriter> bamWriters() { return mBamWriters; }

    public String unmappedSortedBamFilename() { return mUnmappingSortedBamFilename; }

    public long totalWrittenReads()
    {
        return mBamWriters.stream().mapToLong(x -> x.nonConsensusWriteCount()).sum();
    }

    public long sortedWrittenReads()
    {
        return mBamWriters.stream().filter(x -> x.isSorted()).mapToLong(x -> x.nonConsensusWriteCount()).sum();
    }

    public long fullyUnmappedWrittenReads()
    {
        return mFullUnmappedWriter != null ? mFullUnmappedWriter.unsortedWriteCount() : 0;
    }

    public long sortedBamUnsortedWriteCount()
    {
        return mBamWriters.stream().filter(x -> x.isSorted()).mapToLong(x -> x.unsortedWriteCount()).sum();
    }

    public boolean prepareSortedUnmappingBam()
    {
        if(mUnmappingWriter == null)
            return true;

        mUnmappingWriter.close();

        if(mConfig.BamToolPath == null)
            return true;

        if(!BamOperations.sortBam(bamToolName(), bamToolPath(), mUnmappingWriter.filename(), mUnmappingSortedBamFilename, mConfig.Threads))
            return false;

        return BamOperations.indexBam(bamToolName(), bamToolPath(), mUnmappingSortedBamFilename, mConfig.Threads);
    }

    public BamToolName bamToolName() { return BamToolName.fromPath(mConfig.BamToolPath); }
    private String bamToolPath() { return mConfig.BamToolPath; }

    public boolean finaliseBams()
    {
        mReadDataWriter.close();

        if(mFullUnmappedWriter != null)
            mFullUnmappedWriter.close();

        // last thing to do is write fully unmapped read to the final BAM

        if(mConfig.WriteBam && mConfig.BamToolPath != null)
        {
            if(!mConfig.ParallelConcatenation)
                concatenateBams();

            RD_LOGGER.debug("indexing BAM: {}", mFinalBamFilename);

            if(!BamOperations.indexBam(bamToolName(), mConfig.BamToolPath, mFinalBamFilename, mConfig.Threads))
                return false;

            RD_LOGGER.debug("final BAM complete: {}", mFinalBamFilename);
        }

        deleteInterimBams();

        return true;
    }

    private boolean concatenateBams()
    {
        List<String> orderPartitionBams = mBamWriters.stream()
                .filter(x -> x.isSorted()).map(x -> x.filename()).collect(Collectors.toList());

        orderPartitionBams.add(mFullUnmappedWriter.filename());

        RD_LOGGER.debug("concatenating {} BAMs", orderPartitionBams.size());

        if(!BamOperations.concatenateBams(bamToolName(), mConfig.BamToolPath, mFinalBamFilename, orderPartitionBams, 1))
            return false;

        RD_LOGGER.debug("final concatenate complete: {}", mFinalBamFilename);

        return true;
    }

    private void deleteInterimBams()
    {
        if(mConfig.KeepInterimBams || !mConfig.WriteBam)
            return;

        List<String> interimBams = Lists.newArrayList();

        bamWriters().stream().filter(x -> !x.filename().equals(mFinalBamFilename)).forEach(x -> interimBams.add(x.filename()));
        interimBams.add(mUnmappingSortedBamFilename);

        try
        {
            for(String filename : interimBams)
            {
                Files.deleteIfExists(Paths.get(filename));
                Files.deleteIfExists(Paths.get(filename + BAM_INDEX_EXTENSION));
            }
        }
        catch(IOException e)
        {
            RD_LOGGER.error("error deleting interim bams: {}", e.toString());
        }
    }

    private BamWriter createBamWriter(final String filename, boolean synchronousUnsorted)
    {
        SAMFileWriter samFileWriter = null;

        if(mConfig.WriteBam)
        {
            RD_LOGGER.trace("writing temp BAM file: {}", filenamePart(filename));

            // no option to use library-based sorting
            samFileWriter = initialiseSamFileWriter(filename, !synchronousUnsorted);
        }

        // initiate the applicable type of BAM writer - synchronised or not
        BamWriter bamWriter;

        if(synchronousUnsorted)
        {
            bamWriter = new BamWriterSync(filename, mConfig, mReadDataWriter, samFileWriter, mJitterAnalyser);
        }
        else
        {
            bamWriter = new BamWriterNoSync(filename, mConfig, mReadDataWriter, samFileWriter, mJitterAnalyser);
        }

        mBamWriters.add(bamWriter);
        return bamWriter;
    }

    private String formBamFilename(@Nullable final String sorted, @Nullable final String multiId)
    {
        if(!mConfig.MultiBam && mConfig.Threads == 1 && mConfig.OutputBam != null && mConfig.BamToolPath == null)
            return mConfig.OutputBam; // no need to write a temporary BAM

        String filename = mConfig.OutputDir + mConfig.SampleId + "." + FILE_ID;

        if(mConfig.OutputId != null)
            filename += "." + mConfig.OutputId;

        if(multiId != null)
            filename += "." + multiId;

        if(sorted != null)
            filename += "." + sorted;

        filename += BAM_EXTENSION;

        return filename;
    }

    public SAMFileWriter initialiseSamFileWriter(final String filename, boolean isSorted)
    {
        SAMFileHeader fileHeader = buildCombinedHeader(mConfig.BamFiles, mConfig.RefGenomeFile);

        // note that while the sort order may be set to coordinate, the BAM writer is marked as presorted so
        // the BAM will not actually be sorted by the SAMTools library

        // even for the sorted BAM writer, specify unsorted to avoid the presorted last vs next alignment order check in SAMFileWriterImpl
        fileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

        /*
        if(isSorted)
            fileHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        else
            fileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        */

        boolean presorted = isSorted;

        // makeSAMWriter
        return new SAMFileWriterFactory().makeBAMWriter(fileHeader, presorted, new File(filename));
    }
}
