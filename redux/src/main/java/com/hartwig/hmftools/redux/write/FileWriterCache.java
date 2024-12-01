package com.hartwig.hmftools.redux.write;

import static com.hartwig.hmftools.common.bamops.BamMerger.buildCombinedHeader;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BAM_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.filenamePart;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.common.Constants.FILE_ID;

import java.io.File;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bamops.BamOperations;
import com.hartwig.hmftools.common.bamops.BamToolName;
import com.hartwig.hmftools.common.basequal.jitter.JitterAnalyser;
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
    private final BamWriter mSharedUnsortedWriter;

    private String mSortedBamFilename;

    private final JitterAnalyser mJitterAnalyser;

    private static final String SORTED_ID = "sorted";
    private static final String UNSORTED_ID = "unsorted";

    public FileWriterCache(final ReduxConfig config, @Nullable final JitterAnalyser jitterAnalyser)
    {
        mConfig = config;
        mJitterAnalyser = jitterAnalyser;

        mReadDataWriter = new ReadDataWriter(mConfig);

        mBamWriters = Lists.newArrayList();
        mSortedBamFilename = "";

        // create a shared BAM writer if either no multi-threading or using the sorted BAM writer
        if(!mConfig.WriteBam)
        {
            BamWriter bamWriter = new BamWriterNone("none", mConfig, mReadDataWriter, jitterAnalyser);
            mSharedUnsortedWriter = bamWriter;
            mBamWriters.add(bamWriter);
            return;
        }

        mSharedUnsortedWriter = createBamWriter(null, true);
    }

    public BamWriter getPartitionBamWriter(final String fileId)
    {
        if(!mConfig.MultiBam)
            return mSharedUnsortedWriter;

        return createBamWriter(fileId, false);
    }

    public BamWriter getUnsortedBamWriter() { return mSharedUnsortedWriter; }
    public ReadDataWriter readDataWriter() { return mReadDataWriter; }
    public List<BamWriter> bamWriters() { return mBamWriters; }

    public String sortedBamFilename() { return mSortedBamFilename; }

    public long totalWrittenReads()
    {
        return mBamWriters.stream().mapToLong(x -> x.nonConsensusWriteCount()).sum();
    }

    public long sortedBamUnsortedWriteCount()
    {
        return mBamWriters.stream().filter(x -> x.isSorted()).mapToLong(x -> x.unsortedWriteCount()).sum();
    }

    public void closeBams()
    {
        mBamWriters.forEach(x -> x.close());
    }

    public void close() { mReadDataWriter.close();  }

    private BamWriter createBamWriter(@Nullable final String multiId, boolean synchronousUnsorted)
    {
        SAMFileWriter samFileWriter = null;
        String filename = null;

        if(mConfig.WriteBam)
        {
            filename = formBamFilename(synchronousUnsorted ? UNSORTED_ID : SORTED_ID, multiId);

            if(multiId == null)
            {
                RD_LOGGER.trace("writing BAM file: {}", filenamePart(filename));

            }
            else
            {
                RD_LOGGER.trace("writing temp BAM file: {}", filenamePart(filename));
            }

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
            bamWriter = new BamWriterNoSync(
                    filename, mConfig, mReadDataWriter, samFileWriter, mJitterAnalyser, mSharedUnsortedWriter);
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
        if(isSorted)
            fileHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        else
            fileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

        boolean presorted = isSorted;
        return new SAMFileWriterFactory().makeBAMWriter(fileHeader, presorted, new File(filename));
    }

    private BamToolName bamToolName() { return BamToolName.fromPath(mConfig.BamToolPath); }
    private String bamToolPath() { return mConfig.BamToolPath; }

    public boolean sortUnsortedBam()
    {
        String unsortedBamFilename = mBamWriters.get(0).filename();
        mSortedBamFilename = unsortedBamFilename.replaceAll(UNSORTED_ID, SORTED_ID);

        if(!BamOperations.sortBam(bamToolName(), bamToolPath(), unsortedBamFilename, mSortedBamFilename, mConfig.Threads))
            return false;

        return BamOperations.indexBam(bamToolName(), bamToolPath(), mSortedBamFilename, mConfig.Threads);
    }
}
