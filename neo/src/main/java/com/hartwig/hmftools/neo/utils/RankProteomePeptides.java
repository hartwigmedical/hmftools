package com.hartwig.hmftools.neo.utils;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader.convertAminoAcidsToGeneMap;
import static com.hartwig.hmftools.common.utils.ConfigUtils.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.loadDelimitedIdFile;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.NeoCommon.OUTPUT_ID;
import static com.hartwig.hmftools.neo.NeoCommon.THREADS;
import static com.hartwig.hmftools.neo.bind.BindCommon.AMINO_ACID_21ST;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_ALLELE;
import static com.hartwig.hmftools.neo.bind.FlankCounts.FLANK_BASE_COUNT;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.neo.bind.BindData;
import com.hartwig.hmftools.neo.bind.BindScorer;
import com.hartwig.hmftools.neo.bind.ScoreConfig;
import com.hartwig.hmftools.neo.bind.TrainConfig;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class RankProteomePeptides
{
    private final Map<String,List<TranscriptAminoAcids>> mTransAminoAcidMap;

    private final List<String> mAlleles;
    private final BindScorer mScorer;
    private final double mRankCuttoff;
    private final int mThreads;
    private final BufferedWriter mWriter;

    private static final String ALLELE_FILE = "allele_file";
    private static final String RANK_CUTOFF = "rank_cutoff";

    protected static final List<Integer> RANKED_PROTEOME_PEPTIDE_LENGTHS = Lists.newArrayList(8, 9, 10, 11);

    public RankProteomePeptides(final CommandLine cmd)
    {
        mAlleles = loadDelimitedIdFile(cmd.getOptionValue(ALLELE_FILE), FLD_ALLELE, CSV_DELIM);

        Map<String,TranscriptAminoAcids> transAminoAcidMap = Maps.newHashMap();
        EnsemblDataLoader.loadTranscriptAminoAcidData(cmd.getOptionValue(ENSEMBL_DATA_DIR), transAminoAcidMap, Lists.newArrayList(), false);
        mTransAminoAcidMap = convertAminoAcidsToGeneMap(transAminoAcidMap);

        mRankCuttoff = Double.parseDouble(cmd.getOptionValue(RANK_CUTOFF, "0.01"));
        mScorer = new BindScorer(new ScoreConfig(cmd));

        mThreads = Integer.parseInt(cmd.getOptionValue(THREADS, "0"));

        mWriter = initialiseWriter(parseOutputDir(cmd), cmd.getOptionValue(OUTPUT_ID));
    }

    public void run()
    {
        if(mAlleles.isEmpty())
        {
            NE_LOGGER.error("no allele loaded");
            System.exit(1);
        }

        if(!mScorer.loadScoringData())
        {
            NE_LOGGER.error("failed to load scoring data");
            System.exit(1);
        }

        NE_LOGGER.info("searching for {} peptides", mAlleles.size());

        List<PeptideRankTask> searchTasks = Lists.newArrayList();

        if(mThreads > 1)
        {
            int threads = min(mThreads, mAlleles.size());

            for(int i = 0; i < threads; ++i)
            {
                searchTasks.add(new PeptideRankTask(i, mTransAminoAcidMap, mScorer, mRankCuttoff, mWriter));
            }

            int taskIndex = 0;
            for (final String allele : mAlleles)
            {
                searchTasks.get(taskIndex).getAlleles().add(allele);
                ++taskIndex;

                if(taskIndex >= threads)
                    taskIndex = 0;
            }

            final List<Callable> callableList = searchTasks.stream().collect(Collectors.toList());
            TaskExecutor.executeTasks(callableList, callableList.size());
        }
        else
        {
            PeptideRankTask searchTask = new PeptideRankTask(0, mTransAminoAcidMap, mScorer, mRankCuttoff, mWriter);
            searchTask.getAlleles().addAll(mAlleles);
            searchTasks.add(searchTask);
            searchTask.run();
        }

        closeBufferedWriter(mWriter);

        NE_LOGGER.info("allele ranked-peptide search complete");
    }

    private BufferedWriter initialiseWriter(final String outputDir, final String outputId)
    {
        try
        {
            String outputFile = outputDir + outputId + "_allele_ranked_ref_peptides.csv";

            BufferedWriter writer = createBufferedWriter(outputFile, false);

            writer.write("Allele,Peptide,LikelihoodRank,GeneName,TransName");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write results: {}", e.toString());
            return null;
        }
    }

    public synchronized static void writePeptides(final BufferedWriter writer, final String allele, final List<PeptideData> peptideDataList)
    {
        try
        {
            for(PeptideData peptideData : peptideDataList)
            {
                writer.write(String.format("%s,%s,%.6f,%s,%s",
                        allele, peptideData.Peptide, peptideData.LikelihoodRank,
                        peptideData.TransData.GeneName, peptideData.TransData.TransName));

                writer.newLine();
            }
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write results: {}", e.toString());
        }
    }

    private class PeptideRankTask implements Callable
    {
        private final int mTaskId;
        private final Map<String,List<TranscriptAminoAcids>> mTransAminoAcidMap;

        private final BindScorer mScorer;
        private final double mRankCuttoff;
        private final BufferedWriter mWriter;

        private final List<String> mAlleles;

        public PeptideRankTask(
                int taskId, final Map<String,List<TranscriptAminoAcids>> transAminoAcidMap,
                final BindScorer scorer, final double rankCuttoff, final BufferedWriter writer)
        {
            mTaskId = taskId;
            mTransAminoAcidMap = transAminoAcidMap;
            mRankCuttoff = rankCuttoff;
            mScorer = scorer;
            mWriter = writer;

            mAlleles = Lists.newArrayList();
        }

        @Override
        public Long call()
        {
            run();
            return (long)1;
        }

        public List<String> getAlleles() { return mAlleles; }

        public void run()
        {
            NE_LOGGER.info("{}: searching for ranked-peptides for {} alleles", mTaskId, mAlleles.size());

            for(int i = 0; i < mAlleles.size(); ++i)
            {
                findPeptides(mAlleles.get(i));
            }

            NE_LOGGER.info("{}: allele ranked peptide search complete", mTaskId, mAlleles.size());
        }

        private void findPeptides(final String allele)
        {
            List<PeptideData> results = Lists.newArrayList();

            for(int peptideLength : RANKED_PROTEOME_PEPTIDE_LENGTHS)
            {
                for(List<TranscriptAminoAcids> transAaList : mTransAminoAcidMap.values())
                {
                    Set<String> uniquePeptides = Sets.newHashSet();

                    for(TranscriptAminoAcids transAminoAcids : transAaList)
                    {
                        final String aminoAcids = transAminoAcids.AminoAcids;
                        int aaLength = aminoAcids.length();

                        for(int startIndex = 0; startIndex < aaLength - peptideLength; ++startIndex)
                        {
                            int endIndex = startIndex + peptideLength;
                            String aaPeptide = aminoAcids.substring(startIndex, endIndex);

                            if(uniquePeptides.contains(aaPeptide))
                                continue;

                            uniquePeptides.add(aaPeptide);

                            if(aaPeptide.contains(AMINO_ACID_21ST))
                                continue;

                            String upFlank = "";
                            String downFlank = "";

                            int upFlankBases = min(startIndex, FLANK_BASE_COUNT);

                            if(upFlankBases > 0)
                            {
                                upFlank = transAminoAcids.AminoAcids.substring(startIndex - upFlankBases, startIndex);
                            }

                            int downFlankBases = min(transAminoAcids.AminoAcids.length() - endIndex - 1, FLANK_BASE_COUNT);

                            if(downFlankBases > 0)
                            {
                                downFlank = transAminoAcids.AminoAcids.substring(endIndex, endIndex + downFlankBases);
                            }

                            BindData bindData = new BindData(allele, aaPeptide, "", upFlank, downFlank);

                            mScorer.calcScoreData(bindData);

                            if(bindData.likelihoodRank() < mRankCuttoff)
                            {
                                results.add(new PeptideData(aaPeptide, bindData.likelihoodRank(), transAminoAcids));
                            }
                        }
                    }
                }
            }

            NE_LOGGER.debug("{}: allele({}) search found {} ranked peptides", mTaskId, allele, results.size());

            writePeptides(mWriter, allele, results);
        }
    }

    private class PeptideData
    {
        public final String Peptide;
        public final double LikelihoodRank;
        public final TranscriptAminoAcids TransData;

        public PeptideData(final String peptide, double rank, final TranscriptAminoAcids transData)
        {
            Peptide = peptide;
            LikelihoodRank = rank;
            TransData = transData;
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        options.addOption(ALLELE_FILE, true, "Peptides to search for");
        options.addOption(RANK_CUTOFF, true, "Number of amino acid flanks to retrieve");
        options.addOption(THREADS, true, "Threads (default none)");
        options.addOption(OUTPUT_ID, true, "Output file identifier");
        TrainConfig.addCmdLineArgs(options);
        addEnsemblDir(options);
        addLoggingOptions(options);
        addOutputDir(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        RankProteomePeptides rankProteomePeptides = new RankProteomePeptides(cmd);
        rankProteomePeptides.run();
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
