package com.hartwig.hmftools.neo.utils;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.DELIM;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_PEPTIDE;
import static com.hartwig.hmftools.neo.bind.BinderConfig.OUTPUT_ID;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.utils.TaskExecutor;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class RefGenomePeptideLocator
{
    private final int mFlankLength;
    private final int mThreads;

    private final Map<String,TranscriptAminoAcids> mTransAminoAcidMap;

    private final List<String> mPeptides;

    private BufferedWriter mWriter;

    private static final String FLANK_LENGTH = "flank_length";
    private static final String PEPTIDE_FILE = "peptide_file";
    private static final String THREADS = "threads";

    public RefGenomePeptideLocator(final CommandLine cmd)
    {
        String ensemblDataDir = cmd.getOptionValue(ENSEMBL_DATA_DIR);

        mPeptides = Lists.newArrayList();
        mTransAminoAcidMap = Maps.newHashMap();
        EnsemblDataLoader.loadTranscriptAminoAcidData(ensemblDataDir, mTransAminoAcidMap, Lists.newArrayList());

        mFlankLength = Integer.parseInt(cmd.getOptionValue(FLANK_LENGTH));
        mThreads = Integer.parseInt(cmd.getOptionValue(THREADS));

        loadPeptides(cmd.getOptionValue(PEPTIDE_FILE));

        mWriter = initialiseWriter(parseOutputDir(cmd), cmd.getOptionValue(OUTPUT_ID));
    }

    private void loadPeptides(final String filename)
    {
        if(filename == null)
            return;

        try
        {
            List<String> lines = Files.readAllLines(new File(filename).toPath());
            String header = lines.get(0);
            lines.remove(0);

            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DELIM);
            int peptideIndex = fieldsIndexMap.get(FLD_PEPTIDE);

            for(String line : lines)
            {
                String[] values = line.split(DELIM);
                mPeptides.add(values[peptideIndex]);
            }
        }
        catch (IOException e)
        {
            NE_LOGGER.warn("failed to load alleles file: {}", e.toString());
        }
    }

    public void run()
    {
        if(mPeptides.isEmpty())
        {
            NE_LOGGER.error("no peptides loaded");
            System.exit(1);
        }

        NE_LOGGER.info("searching for {} peptides", mPeptides.size());

        List<PeptideSearchTask> searchTasks = Lists.newArrayList();

        if(mThreads > 1)
        {
            List<List<String>> taskPeptideLists = Lists.newArrayList();

            int threads = min(mThreads, mPeptides.size());

            for(int i = 0; i < threads; ++i)
            {
                List<String> peptideList = Lists.newArrayList();
                taskPeptideLists.add(peptideList);
                searchTasks.add(new PeptideSearchTask(i, mTransAminoAcidMap, peptideList, mFlankLength, mWriter));
            }

            int taskIndex = 0;
            for (final String peptide : mPeptides)
            {
                taskPeptideLists.get(taskIndex).add(peptide);
                ++taskIndex;

                if(taskIndex >= threads)
                    taskIndex = 0;
            }

            final List<Callable> callableList = searchTasks.stream().collect(Collectors.toList());
            TaskExecutor.executeTasks(callableList, callableList.size());
        }
        else
        {
            PeptideSearchTask searchTask = new PeptideSearchTask(0, mTransAminoAcidMap, mPeptides, mFlankLength, mWriter);
            searchTasks.add(searchTask);
            searchTask.run();
        }

        int totalFound = searchTasks.stream().mapToInt(x -> x.peptidesFound()).sum();

        NE_LOGGER.info("peptide search complete, found {} of {}", totalFound, mPeptides.size());

        closeBufferedWriter(mWriter);
    }

    private BufferedWriter initialiseWriter(final String outputDir, final String outputId)
    {
        try
        {
            String outputFile = outputDir + "peptide_search";

            if(outputId != null)
                outputFile += "_" + outputId;

            outputFile += ".csv";

            BufferedWriter writer = createBufferedWriter(outputFile, false);

            writer.write("Peptide,GeneId,GeneName,TransId,AminoAcidPos,UpFlank,DownFlank");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to initialise CSV file output: {}", e.toString());
            return null;
        }
    }

    protected synchronized void writeData(final BufferedWriter writer, final TranscriptAminoAcids transAminoAcids, final String peptide,
            int aaPosition, final String upFlank, final String downFlank)
    {
        try
        {
            if(transAminoAcids != null)
            {
                mWriter.write(String.format("%s,%s,%s,%d,%s,%s",
                        peptide, transAminoAcids.GeneId, transAminoAcids.GeneName, aaPosition, upFlank, downFlank));
            }
            else
            {
                mWriter.write(String.format("%s,NONE,NONE,-1,,",peptide));
            }

            mWriter.newLine();
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write peptide search data: {}", e.toString());
        }
    }

    private class PeptideSearchTask implements Callable
    {
        private final int mTaskId;
        private final Map<String,TranscriptAminoAcids> mTransAminoAcidMap;
        private final List<String> mPeptides;
        private final int mFlankLength;
        private final BufferedWriter mWriter;
        private int mFound;

        public PeptideSearchTask(
                int taskId, final Map<String, TranscriptAminoAcids> transAminoAcidMap, final List<String> peptides,
                final int flankLength, final BufferedWriter writer)
        {
            mTaskId = taskId;
            mFlankLength = flankLength;
            mTransAminoAcidMap = transAminoAcidMap;
            mPeptides = peptides;
            mWriter = writer;
            mFound = 0;
        }

        @Override
        public Long call()
        {
            run();
            return (long)1;
        }

        public int peptidesFound() { return mFound; }

        public void run()
        {
            NE_LOGGER.info("{}: searching for {} peptides", mTaskId, mPeptides.size());

            for(int i = 0; i < mPeptides.size(); ++i)
            {
                findPeptide(mPeptides.get(i));

                if(i > 0 && (i % 1000) == 0)
                {
                    NE_LOGGER.info("{}: search count: {}", mTaskId, i);
                }
            }

            NE_LOGGER.info("{}: peptide search complete, found({} of {})", mTaskId, mFound, mPeptides.size());
        }

        private void findPeptide(final String peptide)
        {
            for(TranscriptAminoAcids transAminoAcids : mTransAminoAcidMap.values())
            {
                int aaIndex = transAminoAcids.AminoAcids.indexOf(peptide);
                if(aaIndex < 0)
                    continue;

                ++mFound;

                String upFlank = "";
                String downFlank = "";

                if(mFlankLength > 0)
                {
                    int upFlankBases = min(aaIndex, mFlankLength);

                    if(upFlankBases > 0)
                    {
                        upFlank = transAminoAcids.AminoAcids.substring(aaIndex - upFlankBases, aaIndex);
                    }

                    int downFlankStartPos = aaIndex + peptide.length();
                    int downFlankBases = min(transAminoAcids.AminoAcids.length() - downFlankStartPos - 1, mFlankLength);

                    if(downFlankBases > 0)
                    {
                        downFlank = transAminoAcids.AminoAcids.substring(downFlankStartPos, downFlankStartPos + downFlankBases);
                    }
                }

                writeData(mWriter, transAminoAcids, peptide, aaIndex, upFlank, downFlank);
                return;

                // NE_LOGGER.info("found {} random peptides from {} coding transcripts", totalPeptideCount, transCodingCount);
            }

            writeData(mWriter, null, peptide, -1, "", "");
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        options.addOption(PEPTIDE_FILE, true, "Peptides to search for");
        options.addOption(FLANK_LENGTH, true, "Number of amino acid flanks to retrieve");
        options.addOption(THREADS, true, "Threads (default none)");
        options.addOption(OUTPUT_ID, true, "Output file identifier");
        addEnsemblDir(options);
        addLoggingOptions(options);
        addOutputDir(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        RefGenomePeptideLocator refGenomePeptideLocator = new RefGenomePeptideLocator(cmd);
        refGenomePeptideLocator.run();
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
