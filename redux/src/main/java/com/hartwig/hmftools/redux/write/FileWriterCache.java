package com.hartwig.hmftools.redux.write;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BAM_INDEX_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.filenamePart;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.BamOperations;
import com.hartwig.hmftools.common.bam.BamToolName;
import com.hartwig.hmftools.common.basequal.jitter.JitterAnalyser;
import com.hartwig.hmftools.redux.ReduxConfig;
import com.hartwig.hmftools.redux.merge.BamMerger;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class FileWriterCache
{
    private final ReduxConfig mConfig;
    private final ReadDataWriter mReadDataWriter;

    private final List<BamWriter> mBamWriters;
    private final BamWriter mSharedUnsortedWriter;

    private final JitterAnalyser mJitterAnalyser;

    private static final String BAM_FILE_ID = "redux";
    private static final String SORTED_ID = "sorted";
    private static final String UNSORTED_ID = "unsorted";

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
            mSharedUnsortedWriter = bamWriter;
            mBamWriters.add(bamWriter);
            return;
        }

        mSharedUnsortedWriter = createBamWriter(null, true, false);
    }

    public BamWriter getPartitionBamWriter(final String fileId)
    {
        if(!mConfig.MultiBam)
            return mSharedUnsortedWriter;

        return createBamWriter(fileId, false, true);
    }

    public BamWriter getUnsortedBamWriter()
    {
        if(!mConfig.MultiBam)
            return mSharedUnsortedWriter;

        return mBamWriters.get(0);
    }

    public long totalWrittenReads()
    {
        return mBamWriters.stream().mapToLong(x -> x.nonConsensusWriteCount()).sum();
    }

    public void close()
    {
        mReadDataWriter.close();
        mBamWriters.forEach(x -> x.close());
    }

    private BamWriter createBamWriter(@Nullable final String multiId, boolean isSynchronous, boolean isSorted)
    {
        SAMFileWriter samFileWriter = null;
        String filename = null;

        if(mConfig.WriteBam)
        {
            filename = formBamFilename(isSorted ? SORTED_ID : UNSORTED_ID, multiId);

            if(multiId == null)
            {
                RD_LOGGER.debug("writing BAM file: {}", filenamePart(filename));

            }
            else
            {
                RD_LOGGER.debug("writing temp BAM file: {}", filenamePart(filename));
            }

            // no option to use library-based sorting
            samFileWriter = initialiseSamFileWriter(filename, isSorted);
        }

        // initiate the applicable type of BAM writer - synchronised or not
        BamWriter bamWriter;

        if(isSynchronous)
        {
            bamWriter = new BamWriterSync(filename, mConfig, mReadDataWriter, samFileWriter, mJitterAnalyser);
        }
        else
        {
            bamWriter = new BamWriterNoSync(
                    filename, mConfig, mReadDataWriter, samFileWriter, mJitterAnalyser, isSorted, mSharedUnsortedWriter);
        }

        mBamWriters.add(bamWriter);
        return bamWriter;
    }

    private String formBamFilename(@Nullable final String sorted, @Nullable final String multiId)
    {
        if(!mConfig.MultiBam && mConfig.Threads == 1 && mConfig.OutputBam != null && !runSortMergeIndex())
            return mConfig.OutputBam; // no need to write a temporary BAM

        String filename = mConfig.OutputDir + mConfig.SampleId + "." + BAM_FILE_ID;

        if(mConfig.OutputId != null)
            filename += "." + mConfig.OutputId;

        if(multiId != null)
            filename += "." + multiId;

        if(sorted != null)
            filename += "." + sorted;

        filename += ".bam";

        return filename;
    }

    private SAMFileWriter initialiseSamFileWriter(final String filename, boolean isSorted)
    {
        SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile))
                .open(new File(mConfig.BamFiles.get(0)));

        SAMFileHeader fileHeader = samReader.getFileHeader().clone();

        if(mConfig.BamFiles.size() > 1)
        {
            for(int i = 1; i < mConfig.BamFiles.size(); ++i)
            {
                SamReader nextReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile))
                        .open(new File(mConfig.BamFiles.get(i)));

                for(SAMReadGroupRecord readGroupRecord : nextReader.getFileHeader().getReadGroups())
                {
                    if(!fileHeader.getReadGroups().contains(readGroupRecord))
                        fileHeader.addReadGroup(readGroupRecord);
                }

                final SAMProgramRecord nextProgramRecord = nextReader.getFileHeader().getProgramRecords().get(0);
                String newProgramId = String.format("%s.%d", nextProgramRecord.getId(), i);

                fileHeader.addProgramRecord(new SAMProgramRecord(newProgramId, nextProgramRecord));
            }
        }

        // note that while the sort order may be set to coordinate, the BAM writer is marked as presorted so
        // the BAM will not actually be sorted by the SAMTools library
        if(isSorted)
            fileHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        else
            fileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

        boolean presorted = isSorted;
        return new SAMFileWriterFactory().makeBAMWriter(fileHeader, presorted, new File(filename));
    }

    public boolean runSortMergeIndex() { return mConfig.BamToolPath != null; }

    private BamToolName bamToolName() { return BamToolName.fromPath(mConfig.BamToolPath); }
    private String bamToolPath() { return mConfig.BamToolPath; }

    public boolean sortAndIndexBams()
    {
        if(!runSortMergeIndex())
            return true;

        String finalBamFilename = mConfig.OutputBam != null ? mConfig.OutputBam : formBamFilename(null, null);

        // MD_LOGGER.info("sorting, merging and indexing final BAM");

        String unsortedBamFilename = mBamWriters.get(0).filename();

        // collect up interim BAM files to delete and BAMs to merge
        List<String> interimBams = Lists.newArrayList();
        List<String> bamsToMerge = Lists.newArrayList();

        interimBams.add(unsortedBamFilename);

        String sortedBamFilename;

        if(mBamWriters.size() == 1)
        {
            sortedBamFilename = finalBamFilename;
        }
        else
        {
            sortedBamFilename = unsortedBamFilename.replaceAll(UNSORTED_ID, SORTED_ID);
            bamsToMerge.add(sortedBamFilename);
            interimBams.add(sortedBamFilename);

            for(BamWriter bamWriter : mBamWriters)
            {
                if(bamWriter.isSorted())
                {
                    interimBams.add(bamWriter.filename());
                    bamsToMerge.add(bamWriter.filename());
                }
            }
        }

        // sort the unsorted sync'ed BAM
        SortBamTask sortBamTask = new SortBamTask(unsortedBamFilename, sortedBamFilename, mConfig.Threads);
        sortBamTask.call();
        boolean sortingOk = sortBamTask.success();

        if(!sortingOk && mConfig.Threads > 1)
        {
            // try again with a single thread
            RD_LOGGER.debug("reattempting sort with single thread");
            sortBamTask = new SortBamTask(unsortedBamFilename, finalBamFilename, 1);
            sortBamTask.call();
            sortingOk = sortBamTask.success();
        }

        if(!sortingOk)
            return false;

        if(mBamWriters.size() > 1)
        {
            if(!mergeBams(finalBamFilename, bamsToMerge))
                return false;
        }

        if(!mConfig.KeepInterimBams)
            deleteInterimBams(interimBams);

        // no need for indexing since both merge methods now create an index
        // indexFinalBam(finalBamFilename);

        return true;
    }

    private boolean mergeBams(final String finalBamFilename, final List<String> sortedThreadBams)
    {
        if(bamToolName() == BamToolName.SAMBAMBA)
            return BamOperations.mergeBams(bamToolName(), bamToolPath(), finalBamFilename, sortedThreadBams, mConfig.Threads);

        // use internal BAM merge routine
        BamMerger bamMerger = new BamMerger(
                finalBamFilename, sortedThreadBams, mConfig.RefGenomeFile, bamToolPath(), mConfig.Threads, false);
        return bamMerger.merge();
    }

    private void deleteInterimBams(final List<String> interimBams)
    {
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

    private boolean indexFinalBam(String finalBamFilename)
    {
        // no need to index if Sambamba merge was used
        if(bamToolName() == BamToolName.SAMBAMBA && mBamWriters.size() > 1)
            return true;

        return BamOperations.indexBam(bamToolName(), bamToolPath(), finalBamFilename, mConfig.Threads);
    }

    private class SortBamTask implements Callable
    {
        private final String mBamfile;
        private final String mSortedBamfile;
        private final int mThreadCount;
        private boolean mSuccess;

        public SortBamTask(final String bamfile, final String sortedBamfile, final int threadCount)
        {
            mBamfile = bamfile;
            mThreadCount = threadCount;
            mSortedBamfile = sortedBamfile;
            mSuccess = true;
        }

        public boolean success() { return mSuccess; }

        @Override
        public Long call()
        {
            if(mSortedBamfile == null)
            {
                RD_LOGGER.error("invalid bam filename({})", mBamfile);
                mSuccess = false;
                return (long) 0;
            }

            // MD_LOGGER.debug("sorting unsorted bam({}) to sorted bam({})", mBamfile, mSortedBamfile);

            mSuccess = BamOperations.sortBam(bamToolName(), bamToolPath(), mBamfile, mSortedBamfile, mThreadCount);
            return (long)0;
        }
    }
}
