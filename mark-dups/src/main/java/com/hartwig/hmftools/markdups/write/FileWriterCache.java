package com.hartwig.hmftools.markdups.write;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.filenamePart;
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
import com.hartwig.hmftools.markdups.MarkDupsConfig;

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
        if(!mConfig.MultiBam || mConfig.UseSortCache)
        {
            String fileId = mConfig.UseSortCache ? "shared" : null;
            mSharedUnsortedWriter = createBamWriter(fileId, true, false);
        }
        else
        {
            mSharedUnsortedWriter = null;
        }
    }

    public BamWriter getPartitionBamWriter(final String fileId)
    {
        if(!mConfig.MultiBam && !mConfig.UseSortCache)
            return mSharedUnsortedWriter;

        return createBamWriter(fileId, false, mConfig.UseSortCache);
    }

    public int totalWrittenReads()
    {
        return mBamWriters.stream().mapToInt(x -> x.nonConsensusWriteCount()).sum();
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
        SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(new File(mConfig.RefGenomeFile)).open(new File(mConfig.BamFile));

        SAMFileHeader fileHeader = samReader.getFileHeader().clone();

        // note that while the sort order may be set to coordinate, the BAM writer is marked as presorted so
        // the BAM will not actually be sorted by the SAMTools library
        if(isSorted)
            fileHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        else
            fileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);

        boolean presorted = isSorted;
        return new SAMFileWriterFactory().makeBAMWriter(fileHeader, presorted, new File(filename));
    }

    public void sortAndIndexBams()
    {
        if(mConfig.SamToolsPath == null && mConfig.SambambaPath == null)
            return;

        String finalBamFilename = formBamFilename(null, null);

        List<String> interimBams = Lists.newArrayList();
        List<String> sortedThreadBams = Lists.newArrayList();

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
                System.exit(1);
        }

        MD_LOGGER.debug("sort complete");

        boolean finalBamOk = true;

        if(mBamWriters.size() > 1)
            finalBamOk = mergeBams(finalBamFilename, sortedThreadBams);

        if(!mConfig.KeepInterimBams && finalBamOk)
            deleteInterimBams(interimBams);

        if(finalBamOk)
            indexFinalBam(finalBamFilename);
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
