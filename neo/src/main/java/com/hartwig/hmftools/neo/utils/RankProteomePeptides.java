package com.hartwig.hmftools.neo.utils;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader.convertAminoAcidsToGeneMap;
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
import static com.hartwig.hmftools.neo.bind.BindCommon.AMINO_ACID_21ST;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_ALLELE;
import static com.hartwig.hmftools.neo.bind.BindCommon.ITEM_DELIM;
import static com.hartwig.hmftools.neo.bind.FlankCounts.FLANK_AA_COUNT;

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
import com.hartwig.hmftools.neo.bind.BindData;
import com.hartwig.hmftools.neo.bind.BindScorer;
import com.hartwig.hmftools.neo.bind.ScoreConfig;

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

        mThreads = parseThreads(cmd);

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

        NE_LOGGER.info("searching for {} alleles", mAlleles.size());

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
            String outputFile = BindCommon.formFilename(outputDir, "ranked_prot_peptides", outputId);
            BufferedWriter writer = createBufferedWriter(outputFile, false);

            writer.write("Allele,Peptide,LikelihoodRank,Genes,TransNames");
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
                        allele, peptideData.Peptide, peptideData.LikelihoodRank, peptideData.geneNames(), peptideData.transNames()));

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
            Map<String,PeptideData> results = Maps.newHashMap();

            for(int peptideLength : RANKED_PROTEOME_PEPTIDE_LENGTHS)
            {
                for(List<TranscriptAminoAcids> transAaList : mTransAminoAcidMap.values())
                {
                    for(TranscriptAminoAcids transAminoAcids : transAaList)
                    {
                        final String aminoAcids = transAminoAcids.AminoAcids;
                        int aaLength = aminoAcids.length();

                        for(int startIndex = 0; startIndex < aaLength - peptideLength; ++startIndex)
                        {
                            int endIndex = startIndex + peptideLength;
                            String aaPeptide = aminoAcids.substring(startIndex, endIndex);

                            PeptideData peptideData = results.get(aaPeptide);

                            if(peptideData != null)
                            {
                                peptideData.GeneNames.add(transAminoAcids.GeneName);
                                peptideData.TransNames.add(transAminoAcids.TransName);
                                continue;
                            }

                            if(aaPeptide.contains(AMINO_ACID_21ST))
                                continue;

                            String upFlank = "";
                            String downFlank = "";

                            int upFlankLength = min(startIndex, FLANK_AA_COUNT);

                            if(upFlankLength > 0)
                                upFlank = transAminoAcids.AminoAcids.substring(startIndex - upFlankLength, startIndex);

                            int downFlankLength = min(transAminoAcids.AminoAcids.length() - endIndex - 1, FLANK_AA_COUNT);

                            if(downFlankLength > 0)
                                downFlank = transAminoAcids.AminoAcids.substring(endIndex, endIndex + downFlankLength);

                            BindData bindData = new BindData(allele, aaPeptide, "", upFlank, downFlank);

                            mScorer.calcScoreData(bindData);

                            if(bindData.likelihoodRank() < mRankCuttoff)
                            {
                                peptideData = new PeptideData(aaPeptide, bindData.likelihoodRank(), transAminoAcids);
                                peptideData.GeneNames.add(transAminoAcids.GeneName);
                                peptideData.TransNames.add(transAminoAcids.TransName);
                                results.put(aaPeptide, peptideData);
                            }
                        }
                    }
                }
            }

            NE_LOGGER.debug("{}: allele({}) search found {} ranked peptides", mTaskId, allele, results.size());

            writePeptides(mWriter, allele, results.values().stream().collect(Collectors.toList()));
        }
    }

    private class PeptideData
    {
        public final String Peptide;
        public final double LikelihoodRank;

        public final Set<String> GeneNames;
        public final Set<String> TransNames;

        public PeptideData(final String peptide, double rank, final TranscriptAminoAcids transData)
        {
            Peptide = peptide;
            LikelihoodRank = rank;
            GeneNames = Sets.newHashSet();
            TransNames = Sets.newHashSet();
        }

        public String geneNames()
        {
            if(GeneNames.isEmpty())
                return "";

            StringJoiner sj = new StringJoiner(ITEM_DELIM);
            GeneNames.forEach(x -> sj.add(x));
            return sj.toString();
        }

        public String transNames()
        {
            if(TransNames.isEmpty())
                return "";

            StringJoiner sj = new StringJoiner(ITEM_DELIM);
            TransNames.forEach(x -> sj.add(x));
            return sj.toString();
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        options.addOption(ALLELE_FILE, true, "Peptides to search for");
        options.addOption(RANK_CUTOFF, true, "Number of amino acid flanks to retrieve");
        ScoreConfig.addCmdLineArgs(options);
        addEnsemblDir(options);
        addLoggingOptions(options);
        addOutputOptions(options);
        addThreadOptions(options);

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
