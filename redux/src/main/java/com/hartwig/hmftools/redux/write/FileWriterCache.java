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
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bamops.BamOperations;
import com.hartwig.hmftools.common.bamops.BamToolName;
import com.hartwig.hmftools.common.basequal.jitter.JitterAnalyser;
import com.hartwig.hmftools.redux.BamReader;
import com.hartwig.hmftools.redux.ReduxConfig;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;

public class FileWriterCache
{
    private final ReduxConfig mConfig;
    private final ReadDataWriter mReadDataWriter;

    private final List<BamWriter> mBamWriters;

    // private final BamWriter mSharedUnsortedWriter;

    private final BamWriterSync mUnmappingWriter;
    private String mUnmappingSortedBamFilename;

    private final BamWriterSync mFullUnmappedWriter;
    private String mFullUnmappedSortedBamFilename;

    private final JitterAnalyser mJitterAnalyser;

    private static final String SORTED_ID = "sorted";

    // specific BAM names
    private static final String UNMAPPING = "unmapping";
    private static final String UNMAPPING_SORTED = "unmapping_sorted";
    private static final String FULL_UNMAPPED = "full_unmapped";
    private static final String FULL_UNMAPPED_SORTED = "full_unmapped_sorted";

    public FileWriterCache(final ReduxConfig config, @Nullable final JitterAnalyser jitterAnalyser)
    {
        mConfig = config;
        mJitterAnalyser = jitterAnalyser;

        mReadDataWriter = new ReadDataWriter(mConfig);

        mBamWriters = Lists.newArrayList();

        // create a shared BAM writer if either no multi-threading or using the sorted BAM writer
        if(!mConfig.WriteBam)
        {
            BamWriter bamWriter = new BamWriterNone("none", mConfig, mReadDataWriter, jitterAnalyser);
            mBamWriters.add(bamWriter);

            mUnmappingWriter = null;
            mFullUnmappedWriter = null;
            mUnmappingSortedBamFilename = "";
            mFullUnmappedSortedBamFilename = "";
            return;
        }

        if(mConfig.UnmapRegions.enabled())
        {
            String unmappingFilename = formBamFilename(null, UNMAPPING);
            mUnmappingWriter = (BamWriterSync)createBamWriter(unmappingFilename, true);
            mUnmappingSortedBamFilename = formBamFilename(null, UNMAPPING_SORTED);

            String fullyUnmappedFilename = formBamFilename(null, FULL_UNMAPPED);
            mFullUnmappedWriter = (BamWriterSync)createBamWriter(fullyUnmappedFilename, true);
            mFullUnmappedSortedBamFilename = formBamFilename(null, FULL_UNMAPPED_SORTED);
        }
        else
        {
            mUnmappingWriter = null;
            mFullUnmappedWriter = null;
            mUnmappingSortedBamFilename = "";
            mFullUnmappedSortedBamFilename = "";
        }
    }

    public BamWriter getPartitionBamWriter(final String fileId)
    {
        String filename = formBamFilename(SORTED_ID, fileId);
        return createBamWriter(filename, false);
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

    private BamToolName bamToolName() { return BamToolName.fromPath(mConfig.BamToolPath); }
    private String bamToolPath() { return mConfig.BamToolPath; }

    public boolean finaliseBams()
    {
        mReadDataWriter.close();

        mBamWriters.stream().filter(x -> x.isSorted()).forEach(x -> x.close());

        if(mConfig.BamToolPath == null)
            return true;

        // write any originally full-unmapped reads to the relevant BAM ahead of concatenation
        writeFullyUnmappedReads();

        if(!concatenateBams())
            return false;

        deleteInterimBams();

        return true;
    }

    private boolean writeFullyUnmappedReads()
    {
        BamReader bamReader = new BamReader(mConfig.BamFiles, mConfig.RefGenomeFile);

        AtomicLong unmappedCount = new AtomicLong();

        RD_LOGGER.debug("writing originally fully-unmapped reads");

        bamReader.queryUnmappedReads((final SAMRecord record) ->
        {
            mFullUnmappedWriter.writeRecordSync(record);
            unmappedCount.incrementAndGet();
        });

        if(unmappedCount.get() > 0)
        {
            RD_LOGGER.debug("wrote {} fully-unmapped reads", unmappedCount);
        }

        mFullUnmappedWriter.close();

        // sort and index - seems required
        if(!BamOperations.sortBam(bamToolName(), bamToolPath(), mFullUnmappedWriter.filename(), mFullUnmappedSortedBamFilename, mConfig.Threads))
            return false;

        return BamOperations.indexBam(bamToolName(), bamToolPath(), mFullUnmappedSortedBamFilename, mConfig.Threads);
    }

    private boolean concatenateBams()
    {
        List<String> orderPartitionBams = mBamWriters.stream()
                .filter(x -> x.isSorted()).map(x -> x.filename()).collect(Collectors.toList());

        orderPartitionBams.add(mFullUnmappedWriter.filename());

        String finalBamFilename = mConfig.OutputBam != null ? mConfig.OutputBam : formBamFilename(null, null);

        RD_LOGGER.debug("concatenating {} BAMs", orderPartitionBams.size());

        if(!BamOperations.concatenateBams(bamToolName(), mConfig.BamToolPath, finalBamFilename, orderPartitionBams, mConfig.Threads))
            return false;

        if(!BamOperations.indexBam(bamToolName(), mConfig.BamToolPath, finalBamFilename, mConfig.Threads))
            return false;

        return true;
    }

    private void deleteInterimBams()
    {
        if(mConfig.KeepInterimBams)
            return;

        List<String> interimBams = Lists.newArrayList();

        bamWriters().forEach(x -> interimBams.add(x.mFilename));
        interimBams.add(mUnmappingSortedBamFilename);
        interimBams.add(mFullUnmappedSortedBamFilename);

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
        if(isSorted)
            fileHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        else
            fileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

        boolean presorted = isSorted;
        return new SAMFileWriterFactory().makeBAMWriter(fileHeader, presorted, new File(filename));
    }
}
