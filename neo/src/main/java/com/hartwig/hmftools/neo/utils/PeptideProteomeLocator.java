package com.hartwig.hmftools.neo.utils;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.utils.ConfigUtils.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.loadDelimitedIdFile;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_PEPTIDE;
import static com.hartwig.hmftools.neo.bind.BindCommon.ITEM_DELIM;
import static com.hartwig.hmftools.neo.bind.TranscriptExpression.IMMUNE_EXPRESSION_FILE;
import static com.hartwig.hmftools.neo.bind.TranscriptExpression.IMMUNE_EXPRESSION_FILE_CFG;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.neo.bind.BindCommon;
import com.hartwig.hmftools.neo.bind.TranscriptExpression;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class PeptideProteomeLocator
{
    private final int mFlankLength;
    private final int mThreads;
    private final boolean mFindRepeats;

    private final Map<String,TranscriptAminoAcids> mTransAminoAcidMap;
    private final TranscriptExpression mTranscriptExpression;

    private final List<String> mPeptides;

    private BufferedWriter mWriter;

    private static final String FLANK_LENGTH = "flank_length";
    private static final String PEPTIDE_FILE = "peptide_file";
    private static final String FIND_REPEATS = "find_repeats";

    public PeptideProteomeLocator(final CommandLine cmd)
    {
        mPeptides = loadDelimitedIdFile(cmd.getOptionValue(PEPTIDE_FILE), FLD_PEPTIDE, CSV_DELIM);

        mTransAminoAcidMap = Maps.newHashMap();
        EnsemblDataLoader.loadTranscriptAminoAcidData(cmd.getOptionValue(ENSEMBL_DATA_DIR), mTransAminoAcidMap, Lists.newArrayList(), false);

        mTranscriptExpression = new TranscriptExpression(cmd.getOptionValue(IMMUNE_EXPRESSION_FILE));

        mFlankLength = Integer.parseInt(cmd.getOptionValue(FLANK_LENGTH));
        mThreads = parseThreads(cmd);
        mFindRepeats = cmd.hasOption(FIND_REPEATS);

        mWriter = initialiseWriter(parseOutputDir(cmd), cmd.getOptionValue(OUTPUT_ID));
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

                searchTasks.add(new PeptideSearchTask(
                        i, mTransAminoAcidMap, mTranscriptExpression, peptideList, mFindRepeats, mFlankLength, mWriter));
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
            PeptideSearchTask searchTask = new PeptideSearchTask(
                    0, mTransAminoAcidMap, mTranscriptExpression, mPeptides, mFindRepeats, mFlankLength, mWriter);

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
            String outputFile = BindCommon.formFilename(outputDir, "peptide_search", outputId);
            BufferedWriter writer = createBufferedWriter(outputFile, false);

            writer.write("Peptide,Genes,Transcripts,AminoAcidPos,UpFlank,DownFlank,MatchCount,TotalTPM");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to initialise CSV file output: {}", e.toString());
            return null;
        }
    }

    protected synchronized static void writeData(
            final BufferedWriter writer, final String geneNames, final String transNames, final String peptide,
            int aaPosition, final String upFlank, final String downFlank, int repeatCount, double tpmTotal)
    {
        try
        {
            if(geneNames != null)
            {
                writer.write(String.format("%s,%s,%s,%d,%s,%s,%d,%.4f",
                        peptide, geneNames, transNames, aaPosition, upFlank, downFlank, repeatCount, tpmTotal));
            }
            else
            {
                writer.write(String.format("%s,NONE,NONE,-1,-,-,0,0",peptide));
            }

            writer.newLine();
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
        private final TranscriptExpression mTranscriptExpression;
        private final List<String> mPeptides;
        private final int mFlankLength;
        private final BufferedWriter mWriter;
        private final boolean mFindRepeats;
        private int mFound;

        public PeptideSearchTask(
                int taskId, final Map<String,TranscriptAminoAcids> transAminoAcidMap, final TranscriptExpression transcriptExpression,
                final List<String> peptides, boolean findRepeats, final int flankLength, final BufferedWriter writer)
        {
            mTaskId = taskId;
            mFlankLength = flankLength;
            mTransAminoAcidMap = transAminoAcidMap;
            mTranscriptExpression = transcriptExpression;
            mPeptides = peptides;
            mFindRepeats = findRepeats;
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
            int matches = 0;

            int matchedAaIndex = -1;
            String upFlank = "";
            String downFlank = "";
            double tpmTotal = 0;

            Set<String> geneNames = Sets.newHashSet();
            List<String> transNames = Lists.newArrayList();

            for(TranscriptAminoAcids transAminoAcids : mTransAminoAcidMap.values())
            {
                int aaIndex = transAminoAcids.AminoAcids.indexOf(peptide);
                if(aaIndex < 0)
                    continue;

                ++matches;

                if(geneNames.isEmpty())
                {
                    ++mFound;
                    matchedAaIndex = aaIndex;

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
                }

                if(geneNames.size() < 10)
                    geneNames.add(transAminoAcids.GeneName);

                if(geneNames.size() == 1 || transNames.size() < 20)
                    transNames.add(transAminoAcids.TransName);

                if(mTranscriptExpression != null)
                {
                    Double tpm = mTranscriptExpression.getExpression(transAminoAcids.TransName);

                    if(tpm != null)
                        tpmTotal += tpm;
                }

                if(!mFindRepeats)
                    break;

                // NE_LOGGER.info("found {} random peptides from {} coding transcripts", totalPeptideCount, transCodingCount);
            }

            if(!geneNames.isEmpty())
            {
                StringJoiner sjGene = new StringJoiner(ITEM_DELIM);
                StringJoiner sjTrans = new StringJoiner(ITEM_DELIM);
                geneNames.forEach(x -> sjGene.add(x));
                transNames.forEach(x -> sjTrans.add(x));
                writeData(mWriter, sjGene.toString(), sjTrans.toString(), peptide, matchedAaIndex, upFlank, downFlank, matches, tpmTotal);
            }
            else
            {
                writeData(mWriter, null, null, peptide, -1, "", "", 0, 0);
            }
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        options.addOption(PEPTIDE_FILE, true, "Peptides to search for");
        options.addOption(FLANK_LENGTH, true, "Number of amino acid flanks to retrieve");
        options.addOption(FIND_REPEATS, false, "Look for repeated matches, default false");
        options.addOption(IMMUNE_EXPRESSION_FILE, true, IMMUNE_EXPRESSION_FILE_CFG);
        addEnsemblDir(options);
        addLoggingOptions(options);
        addOutputOptions(options);
        addThreadOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        PeptideProteomeLocator peptideProteomeLocator = new PeptideProteomeLocator(cmd);
        peptideProteomeLocator.run();
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
