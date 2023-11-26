package com.hartwig.hmftools.markdups.write;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.filenamePart;
import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.markdups.MarkDupsConfig;

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
    private final MarkDupsConfig mConfig;
    private final ReadDataWriter mReadDataWriter;

    private final List<BamWriter> mBamWriters;
    private final BamWriter mSharedUnsortedWriter;

    private static final String BAM_FILE_ID = "mark_dups";
    private static final String SORTED_ID = "sorted";
    private static final String UNSORTED_ID = "unsorted";
    
    public FileWriterCache(final MarkDupsConfig config)
    {
        mConfig = config;

        mReadDataWriter = new ReadDataWriter(mConfig);

        mBamWriters = Lists.newArrayList();

        // create a shared BAM writer if either no multi-threading or using the sorted BAM writer
        if(!mConfig.WriteBam)
        {
            BamWriter bamWriter = new BamWriterNone("none", mConfig, mReadDataWriter);
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

    public BamWriter getFullyUnmappedReadsBamWriter()
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
                MD_LOGGER.debug("writing BAM file: {}", filenamePart(filename));

            }
            else
            {
                MD_LOGGER.debug("writing temp BAM file: {}", filenamePart(filename));
            }

            // no option to use library-based sorting
            samFileWriter = initialiseSamFileWriter(filename, isSorted);
        }

        // initiate the applicable type of BAM writer - synchronised or not
        BamWriter bamWriter;

        if(isSynchronous)
        {
            bamWriter = new BamWriterSync(filename, mConfig, mReadDataWriter, samFileWriter);
        }
        else
        {
            bamWriter = new BamWriterNoSync(filename, mConfig, mReadDataWriter, samFileWriter, isSorted, mSharedUnsortedWriter);
        }

        mBamWriters.add(bamWriter);
        return bamWriter;
    }

    private String formBamFilename(@Nullable final String sorted, @Nullable final String multiId)
    {
        if(!mConfig.MultiBam && mConfig.Threads == 1 && mConfig.OutputBam != null)
            return mConfig.OutputBam;

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

                nextReader.getFileHeader().getReadGroups().forEach(x -> fileHeader.addReadGroup(x));

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

    public boolean runSortMergeIndex() { return mConfig.SamToolsPath != null || mConfig.SambambaPath != null; }

    public boolean sortAndIndexBams()
    {
        if(!runSortMergeIndex())
            return true;

        String finalBamFilename = mConfig.OutputBam != null ? mConfig.OutputBam : formBamFilename(null, null);

        List<String> interimBams = Lists.newArrayList();
        List<String> sortedThreadBams = Lists.newArrayList();

        if(mConfig.SamToolsPath == null)
        {
            MD_LOGGER.error("samtools required for sort");
            return false;
        }

        MD_LOGGER.info("sorting, merging and indexing final BAM");

        boolean sortingOk = true;

        if(mBamWriters.size() == 1)
        {
            String unsortedBamFilename = mBamWriters.get(0).filename();
            SortBamTask sortBamTask = new SortBamTask(unsortedBamFilename, finalBamFilename, mConfig.Threads);

            MD_LOGGER.debug("sorting bam");
            sortingOk = sortBamTask.success();

            interimBams.add(unsortedBamFilename);
        }
        else
        {
            List<SortBamTask> sortTasks = Lists.newArrayList();

            int unsortedBamCount = (int)mBamWriters.stream().filter(x -> !x.isSorted()).count();
            int maxThreadCount = unsortedBamCount > 1 ? 1 : mConfig.Threads;

            for(BamWriter bamWriter : mBamWriters)
            {
                if(bamWriter.isSorted())
                {
                    interimBams.add(bamWriter.filename());
                    sortedThreadBams.add(bamWriter.filename());
                    continue;
                }

                String sortedBamFile = bamWriter.filename().replaceAll(UNSORTED_ID, SORTED_ID);

                sortTasks.add(new SortBamTask(bamWriter.filename(), sortedBamFile, maxThreadCount));

                interimBams.add(bamWriter.filename());
                interimBams.add(sortedBamFile);
                sortedThreadBams.add(sortedBamFile);
            }

            List<Callable> callableTasks = sortTasks.stream().collect(Collectors.toList());

            MD_LOGGER.debug("sorting {} bam file(s)", sortTasks.size());

            if(!TaskExecutor.executeTasks(callableTasks, mConfig.Threads))
            {
                sortingOk = false;
            }
            else
            {
                sortingOk = sortTasks.stream().allMatch(x -> x.success());
            }
        }

        if(!sortingOk)
            return false;

        MD_LOGGER.debug("sort complete");

        if(mBamWriters.size() > 1)
        {
            if(!mergeBams(finalBamFilename, sortedThreadBams))
                return false;
        }

        if(!mConfig.KeepInterimBams)
            deleteInterimBams(interimBams);

        if(!indexFinalBam(finalBamFilename))
            return false;

        return true;
    }

    private boolean mergeBams(final String finalBamFilename, final List<String> sortedThreadBams)
    {
        MD_LOGGER.debug("merging {} bams", mBamWriters.size());

        final String[] command = new String[5 + sortedThreadBams.size()];

        int index = 0;

        if(mConfig.SambambaPath != null)
        {
            command[index++] = mConfig.SambambaPath;
            command[index++] = "merge";
            command[index++] = "-t";
        }
        else
        {
            command[index++] = mConfig.SamToolsPath;
            command[index++] = "merge";
            command[index++] = "-@";
        }

        command[index++] = String.valueOf(mConfig.Threads);
        command[index++] = finalBamFilename;

        for(String threadBam : sortedThreadBams)
        {
            command[index++] = threadBam;
        }

        if(executeCommand(command, finalBamFilename))
        {
            MD_LOGGER.debug("merge complete");
            return true;
        }

        return false;
    }

    private void deleteInterimBams(final List<String> interimBams)
    {
        try
        {
            for(String filename : interimBams)
            {
                Files.deleteIfExists(Paths.get(filename));
            }
        }
        catch(IOException e)
        {
            MD_LOGGER.error("error deleting interim bams: {}", e.toString());
        }
    }

    private boolean indexFinalBam(String finalBamFilename)
    {
        // no need to index if Sambamba merge was used
        if(mConfig.SambambaPath != null && mBamWriters.size() > 1)
            return true;

        MD_LOGGER.debug("indexing final bam");

        final String[] command = new String[5];

        int index = 0;
        command[index++] = mConfig.SamToolsPath;
        command[index++] = "index";
        command[index++] = "-@";
        command[index++] = String.valueOf(mConfig.Threads);
        command[index++] = finalBamFilename;

        if(!executeCommand(command, finalBamFilename))
            return false;

        MD_LOGGER.debug("index complete");
        return true;
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
                MD_LOGGER.error("invalid bam filename({})", mBamfile);
                return (long)0;
            }

            MD_LOGGER.debug("sorting unsorted bam({}) to sorted bam({})", mBamfile, mSortedBamfile);

            // String sortArgs = format("sort -@ %s -m %dG -T tmp -O bam %s -o %s", Bash.allCpus(), SORT_MEMORY_PER_CORE, inputBam, outputBam);

            final String[] command = new String[11];

            int index = 0;
            command[index++] = mConfig.SamToolsPath;
            command[index++] = "sort";
            command[index++] = "-@";
            command[index++] = String.valueOf(mThreadCount);

            // default memory per thread according to samtools doco is 768MB
            // command[index++] = "-m";
            // command[index++] = "1G"; // could configure as a function of max heap used by MarkDups

            command[index++] = "-T";
            command[index++] = "tmp";
            command[index++] = "-O";
            command[index++] = "bam";
            command[index++] = mBamfile;
            command[index++] = "-o";
            command[index++] = mSortedBamfile;

            if(executeCommand(command, mSortedBamfile))
                return (long)0;

            mSuccess = false;
            return (long)1;
        }
    }

    private static boolean executeCommand(final String[] command, final String outputPrefix)
    {
        String redirectOutputFile = outputPrefix + ".out";
        String redirectErrFile = outputPrefix + ".err";

        try
        {
            int result = new ProcessBuilder(command)
                    .redirectOutput(new File(redirectOutputFile))
                    .redirectError(new File(redirectErrFile))
                    .start().waitFor();

            if(result != 0)
            {
                MD_LOGGER.error("error running command({}) for file({})", commandToStr(command), outputPrefix);
                return false;
            }

            // clean-up process log files
            Files.deleteIfExists(Paths.get(redirectOutputFile));
            Files.deleteIfExists(Paths.get(redirectErrFile));
        }
        catch(Exception e)
        {
            MD_LOGGER.error("error running command({}) for file({}): {}", commandToStr(command), outputPrefix, e.toString());
            return false;
        }

        return true;
    }

    private static String commandToStr(final String[] command)
    {
        return Arrays.stream(command).collect(Collectors.joining(" "));
    }
}
