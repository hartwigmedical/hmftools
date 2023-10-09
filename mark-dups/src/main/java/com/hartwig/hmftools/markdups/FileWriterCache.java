package com.hartwig.hmftools.markdups;

import static com.hartwig.hmftools.markdups.MarkDupsConfig.MD_LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.TaskExecutor;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class FileWriterCache
{
    private final MarkDupsConfig mConfig;
    private final ReadDataWriter mReadDataWriter;

    private final List<BamWriter> mBamWriters;

    private static final String BAM_FILE_ID = "mark_dups";
    private static final String SORTED_ID = "sorted";
    private static final String UNSORTED_ID = "unsorted";

    public FileWriterCache(final MarkDupsConfig config)
    {
        mConfig = config;

        mReadDataWriter = new ReadDataWriter(mConfig);

        mBamWriters = Lists.newArrayList();


        if(!mConfig.MultiBam)
            createBamWriter(null);
    }

    public BamWriter getBamWriter(final String fileId)
    {
        if(!mConfig.MultiBam)
            return mBamWriters.get(0);

        return createBamWriter(fileId);
    }

    public int totalWrittenReads()
    {
        return mBamWriters.stream().mapToInt(x -> x.recordWriteCount()).sum();
    }

    public void logUnwrittenReads()
    {
        mBamWriters.forEach(x -> x.logUnwrittenReads());
    }

    public void close()
    {
        mReadDataWriter.close();

        // closing a sorted BAM involves a final sort, so ensure this is also multi-threaded
        if(mConfig.SortedBam && mBamWriters.size() > 1)
        {
            List<SortBamCloseTask> closedSortedBamTasks = mBamWriters.stream().map(x -> new SortBamCloseTask(x)).collect(Collectors.toList());

            List<Callable> callableTasks = closedSortedBamTasks.stream().collect(Collectors.toList());

            MD_LOGGER.debug("closing {} sorted bam(s)", callableTasks.size());

            if(!TaskExecutor.executeTasks(callableTasks, mConfig.Threads))
                System.exit(1);

            MD_LOGGER.debug("sorted BAM close complete");
        }
        else
        {
            mBamWriters.forEach(x -> x.close());
        }
    }

    public BamWriter createBamWriter(@Nullable final String multiId)
    {
        SAMFileWriter samFileWriter = null;
        String filename = null;

        if(mConfig.WriteBam)
        {
            filename = formBamFilename(mConfig.SortedBam ? SORTED_ID : UNSORTED_ID, multiId);

            if(multiId == null)
            {
                MD_LOGGER.debug("writing BAM file: {}", filename);

            }
            else
            {
                MD_LOGGER.debug("writing tmp BAM file: {}", filename);
            }

            samFileWriter = initialiseSamFileWriter(filename);
        }

        BamWriter bamWriter = new BamWriter(filename, mConfig, mReadDataWriter, samFileWriter);
        mBamWriters.add(bamWriter);
        return bamWriter;
    }

    private SAMFileWriter initialiseSamFileWriter(final String filename)
    {
        SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile));

        SAMFileHeader fileHeader = samReader.getFileHeader().clone();

        if(mConfig.SortedBam)
            fileHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        else
            fileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

        return new SAMFileWriterFactory().makeBAMWriter(fileHeader, false, new File(filename));
    }

    private String formBamFilename(@Nullable final String sorted, @Nullable final String multiId)
    {
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

    public void sortAndIndexBams()
    {
        if(mConfig.SamToolsPath == null && mConfig.SambambaPath == null)
            return;

        String finalBamFilename = formBamFilename(null, null);

        List<String> interimBams = Lists.newArrayList();
        List<String> sortedThreadBams = Lists.newArrayList();

        if(!mConfig.SortedBam)
        {
            if(mConfig.SamToolsPath == null)
            {
                MD_LOGGER.error("samtools required for sort");
                return;
            }

            if(mBamWriters.size() == 1)
            {
                String unsortedBamFilename = mBamWriters.get(0).filename();
                SortBamTask sortBamTask = new SortBamTask(unsortedBamFilename, finalBamFilename, mConfig.Threads);

                MD_LOGGER.debug("sorting bam");
                sortBamTask.call();

                interimBams.add(unsortedBamFilename);
            }
            else
            {
                List<SortBamTask> sortTasks = Lists.newArrayList();

                for(BamWriter bamWriter : mBamWriters)
                {
                    String sortedBamFile = bamWriter.filename().replaceAll(UNSORTED_ID, SORTED_ID);

                    sortTasks.add(new SortBamTask(bamWriter.filename(), sortedBamFile, 1));

                    interimBams.add(bamWriter.filename());
                    interimBams.add(sortedBamFile);
                    sortedThreadBams.add(sortedBamFile);
                }

                List<Callable> callableTasks = sortTasks.stream().collect(Collectors.toList());

                MD_LOGGER.debug("sorting {} bam file(s)", sortTasks.size());

                if(!TaskExecutor.executeTasks(callableTasks, mConfig.Threads))
                    System.exit(1);
            }

            MD_LOGGER.debug("sort complete");
        }
        else
        {
            for(BamWriter bamWriter : mBamWriters)
            {
                interimBams.add(bamWriter.filename());
                sortedThreadBams.add(bamWriter.filename());
            }
        }

        if(mBamWriters.size() > 1)
            mergeBams(finalBamFilename, sortedThreadBams);

        if(!mConfig.KeepInterimBams)
            deleteInterimBams(interimBams);

        indexFinalBam(finalBamFilename);
    }

    private void mergeBams(final String finalBamFilename, final List<String> sortedThreadBams)
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

        executeCommand(command, finalBamFilename);

        MD_LOGGER.debug("merge complete");
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

    private void indexFinalBam(String finalBamFilename)
    {
        // no need to index if Sambamba merge was used
        if(mConfig.SambambaPath != null && mBamWriters.size() > 1)
            return;

        MD_LOGGER.debug("indexing final bam");

        final String[] command = new String[5];

        int index = 0;
        command[index++] = mConfig.SamToolsPath;
        command[index++] = "index";
        command[index++] = "-@";
        command[index++] = String.valueOf(mConfig.Threads);
        command[index++] = finalBamFilename;

        executeCommand(command, finalBamFilename);

        MD_LOGGER.debug("index complete");
    }

    private class SortBamTask implements Callable
    {
        private final String mBamfile;
        private final String mSortedBamfile;
        private final int mThreadCount;

        public SortBamTask(final String bamfile, final String sortedBamfile, final int threadCount)
        {
            mBamfile = bamfile;
            mThreadCount = threadCount;
            mSortedBamfile = sortedBamfile;
        }

        @Override
        public Long call()
        {
            if(mSortedBamfile == null)
            {
                MD_LOGGER.error("invalid bam filename({})", mBamfile);
                return (long)0;
            }

            // String sortArgs = format("sort -@ %s -m %dG -T tmp -O bam %s -o %s", Bash.allCpus(), SORT_MEMORY_PER_CORE, inputBam, outputBam);

            final String[] command = new String[9];

            int index = 0;
            command[index++] = mConfig.SamToolsPath;
            command[index++] = "sort";
            command[index++] = "-@";
            command[index++] = String.valueOf(mThreadCount);
            command[index++] = "-O";
            command[index++] = "bam";
            command[index++] = mBamfile;
            command[index++] = "-o";
            command[index++] = mSortedBamfile;

            executeCommand(command, mSortedBamfile);

            return (long)0;
        }
    }

    private class SortBamCloseTask implements Callable
    {
        private final BamWriter mBamWriter;

        public SortBamCloseTask(final BamWriter bamWriter)
        {
            mBamWriter = bamWriter;
        }

        @Override
        public Long call()
        {
            if(mBamWriter == null)
            {
                MD_LOGGER.error("invalid bam writer");
                return (long)0;
            }

            mBamWriter.close();

            return (long)0;
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
                MD_LOGGER.error("error running command({}:{}) for file({})", command[0], command[1], outputPrefix);
                return false;
            }

            // clean-up process log files
            Files.deleteIfExists(Paths.get(redirectOutputFile));
            Files.deleteIfExists(Paths.get(redirectErrFile));
        }
        catch(Exception e)
        {
            MD_LOGGER.error("error running command({}:{}) for file({}): {}", command[0], command[1], outputPrefix, e.toString());
            return false;
        }

        return true;
    }
}
