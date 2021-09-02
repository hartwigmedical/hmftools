package com.hartwig.hmftools.neo.utils;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.codon.Codons.STOP_AMINO_ACID;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.utils.ConfigUtils.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.loadDelimitedIdFile;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_PEPTIDE;
import static com.hartwig.hmftools.neo.bind.BinderConfig.OUTPUT_ID;
import static com.hartwig.hmftools.neo.bind.BinderConfig.THREADS;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.neo.bind.BlosumMapping;
import com.hartwig.hmftools.neo.bind.RandomPeptideConfig;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

// routine for finding the most similar peptide in the ref genome
public class PeptideRefSimilarity
{
    private final List<String> mPeptides;

    private final Map<String,TranscriptAminoAcids> mTransAminoAcidMap;

    private final String mOutputDir;
    private final String mOutputId;
    private final int mThreads;

    private static final String PEPTIDES_FILE = "peptides_file";
    private static final String AMINO_ACID_21ST = "X";

    public PeptideRefSimilarity(final CommandLine cmd)
    {
        mPeptides = loadDelimitedIdFile(cmd.getOptionValue(PEPTIDES_FILE), FLD_PEPTIDE, CSV_DELIM);

        mTransAminoAcidMap = Maps.newHashMap();
        EnsemblDataLoader.loadTranscriptAminoAcidData(cmd.getOptionValue(ENSEMBL_DATA_DIR), mTransAminoAcidMap, Lists.newArrayList());

        mOutputDir = parseOutputDir(cmd);
        mOutputId = cmd.getOptionValue(OUTPUT_ID);

        mThreads = Integer.parseInt(cmd.getOptionValue(THREADS, "0"));
    }

    public void run()
    {
        NE_LOGGER.info("finding top similarities for {} peptides", mPeptides.size());

        List<PeptideSearchTask> searchTasks = Lists.newArrayList();

        if(mThreads > 1)
        {
            int threads = min(mThreads, mPeptides.size());

            for(int i = 0; i < threads; ++i)
            {
                searchTasks.add(new PeptideSearchTask(i, mTransAminoAcidMap));
            }

            int taskIndex = 0;
            for (final String peptide : mPeptides)
            {
                searchTasks.get(taskIndex).getPeptides().add(peptide);
                ++taskIndex;

                if(taskIndex >= threads)
                    taskIndex = 0;
            }

            final List<Callable> callableList = searchTasks.stream().collect(Collectors.toList());
            TaskExecutor.executeTasks(callableList, callableList.size());
        }
        else
        {
            PeptideSearchTask searchTask = new PeptideSearchTask(0, mTransAminoAcidMap);
            searchTask.getPeptides().addAll(mPeptides);
            searchTasks.add(searchTask);
            searchTask.run();
        }

        try
        {
            String filename = mOutputDir + "peptide_ref_gen_sim_" + mOutputId + ".csv";
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("Peptide,TopSimilarity,TopPeptide,GeneId,GeneName,TransName"); // PepTotal,TopTotal,DiffTotal,
            writer.newLine();

            for(PeptideSearchTask searchTask : searchTasks)
            {
                for(Map.Entry<String,PeptideSimiliarity> entry : searchTask.getResults().entrySet())
                {
                    final PeptideSimiliarity pepSim = entry.getValue();

                    writer.write(String.format("%s,%.1f,%s,%s,%s,%s",
                            entry.getKey(), pepSim.Similiarity, pepSim.Peptide,
                            // pepSim.SimValues[0], pepSim.SimValues[1], pepSim.SimValues[2],
                            pepSim.TransData.GeneId, pepSim.TransData.GeneName, pepSim.TransData.TransName));

                    writer.newLine();
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write peptide similarity results: {}", e.toString());
            return;
        }

        NE_LOGGER.info("similarity search complete");
    }

    private class PeptideSearchTask implements Callable
    {
        private final int mTaskId;
        private final Map<String,TranscriptAminoAcids> mTransAminoAcidMap;
        private final List<String> mPeptides;

        private final BlosumMapping mBlosumMapping;
        private final Map<String,PeptideSimiliarity> mResults;

        public PeptideSearchTask(
                int taskId, final Map<String, TranscriptAminoAcids> transAminoAcidMap)
        {
            mTaskId = taskId;
            mTransAminoAcidMap = transAminoAcidMap;
            mPeptides = Lists.newArrayList();
            mResults = Maps.newHashMap();

            mBlosumMapping = new BlosumMapping();
        }

        public List<String> getPeptides() { return mPeptides; }
        public Map<String,PeptideSimiliarity> getResults() { return mResults; }

        @Override
        public Long call()
        {
            run();
            return (long)1;
        }

        private void run()
        {
            NE_LOGGER.info("{}: searching for {} peptides", mTaskId, mPeptides.size());

            for(int i = 0; i < mPeptides.size(); ++i)
            {
                findPeptideSimilarity(mPeptides.get(i));

                if(i > 0 && (i % 10) == 0)
                {
                    NE_LOGGER.info("{}: search count: {}", mTaskId, i);
                }
            }

            NE_LOGGER.info("{}: peptide search complete", mTaskId);
        }

        private void findPeptideSimilarity(final String peptide)
        {
            String topPeptide = "";
            double topSimiliarity = 0;
            int[] topSimValues = new int[3];
            TranscriptAminoAcids topTrans = null;

            int peptideLength = peptide.length();

            for(TranscriptAminoAcids transAminoAcids : mTransAminoAcidMap.values())
            {
                final String aminoAcids = transAminoAcids.AminoAcids;

                if(aminoAcids.contains(peptide))
                {
                    mResults.put(peptide, new PeptideSimiliarity(peptide, 0, null, transAminoAcids));
                    break;
                }

                int aaLength = aminoAcids.length();

                for(int startIndex = 0; startIndex < aaLength - peptideLength; ++startIndex)
                {
                    String aaPeptide = aminoAcids.substring(startIndex, startIndex + peptideLength);

                    if(aaPeptide.contains(AMINO_ACID_21ST))
                        continue;


                    // int similarity = mBlosumMapping.calcSequenceBlosumDifference(aaPeptide, peptide);
                    // int[] simValues = new int[] {0, 0, 0};

                    double similarity = 0;

                    for(int i = 0; i < peptide.length(); ++i)
                    {
                        char aa1 = peptide.charAt(i);
                        char aa2 = aaPeptide.charAt(i);

                        int bs1 = mBlosumMapping.selfMapping(aa1);
                        int bs2 = mBlosumMapping.selfMapping(aa2);
                        int map = mBlosumMapping.map(aa1, aa2);

                        similarity += (bs1 + bs2) * 0.5 - map;

                        if(topTrans != null && similarity >= topSimiliarity) // skip if cannot be better
                            break;

                        /*
                        simValues[0] += mBlosumMapping.selfMapping(aa1) - diff;
                        simValues[1] += mBlosumMapping.selfMapping(aa2) - diff;
                        simValues[2] += abs(diff);
                        */
                    }

                    if(topTrans == null || similarity < topSimiliarity)
                    {
                        topPeptide = aaPeptide;
                        topSimiliarity = similarity;
                        topTrans = transAminoAcids;
                        // topSimValues = simValues;
                    }
                }
            }

            mResults.put(peptide, new PeptideSimiliarity(topPeptide, topSimiliarity, topSimValues, topTrans));
        }
    }

    private class PeptideSimiliarity
    {
        public final String Peptide;
        public final double Similiarity;
        public final int[] SimValues;
        public final TranscriptAminoAcids TransData;

        public PeptideSimiliarity(final String peptide, double similiarity, final int[] simValues, final TranscriptAminoAcids transData)
        {
            Peptide = peptide;
            Similiarity = similiarity;
            SimValues = simValues;
            TransData = transData;
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        addEnsemblDir(options);
        options.addOption(PEPTIDES_FILE, true, "Peptides file");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(OUTPUT_ID, true, "Output file id");
        options.addOption(THREADS, true, "Threads (default none)");
        options.addOption(LOG_DEBUG, false, "Log verbose");
        addLoggingOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        PeptideRefSimilarity peptideRefSimilarity = new PeptideRefSimilarity(cmd);
        peptideRefSimilarity.run();
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
