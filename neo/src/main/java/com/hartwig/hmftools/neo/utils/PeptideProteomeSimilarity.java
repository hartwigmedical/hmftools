package com.hartwig.hmftools.neo.utils;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader.convertAminoAcidsToGeneMap;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.neo.NeoCommon.APP_NAME;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.AMINO_ACID_21ST;
import static com.hartwig.hmftools.neo.bind.BindCommon.BIND_DELIM;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_ALLELE;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_PEPTIDE;
import static com.hartwig.hmftools.neo.bind.BindCommon.cleanAllele;
import static com.hartwig.hmftools.neo.utils.RankProteomePeptides.RANKED_PROTEOME_PEPTIDE_LENGTHS;
import static com.hartwig.hmftools.neo.utils.RankedProteomePeptides.PROTEOME_RANKS_FILE;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.aminoacid.BlosumMapping;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataLoader;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.neo.bind.BindCommon;
import com.hartwig.hmftools.neo.bind.BindData;
import com.hartwig.hmftools.neo.bind.BindScorer;
import com.hartwig.hmftools.neo.bind.ScoreConfig;

import org.jetbrains.annotations.NotNull;

public class PeptideProteomeSimilarity
{
    private final List<PeptideSimilarity> mPeptideSimilarities;

    private final Map<String, List<TranscriptAminoAcids>> mTransAminoAcidMap;
    private final RankedProteomePeptides mRankedProteomePeptides;
    private final BindScorer mScorer;

    private final String mOutputDir;
    private final String mOutputId;
    private final int mThreads;

    private static final String PEPTIDES_FILE = "peptides_file";

    public PeptideProteomeSimilarity(final ConfigBuilder configBuilder)
    {
        mPeptideSimilarities = Lists.newArrayList();
        loadPeptides(configBuilder.getValue(PEPTIDES_FILE));

        if(configBuilder.hasValue(ENSEMBL_DATA_DIR))
        {
            Map<String, TranscriptAminoAcids> transAminoAcidMap = Maps.newHashMap();
            EnsemblDataLoader.loadTranscriptAminoAcidData(configBuilder.getValue(ENSEMBL_DATA_DIR), transAminoAcidMap, Lists.newArrayList(), false);
            mTransAminoAcidMap = convertAminoAcidsToGeneMap(transAminoAcidMap);
        }
        else
        {
            mTransAminoAcidMap = null;
        }

        if(configBuilder.hasValue(PROTEOME_RANKS_FILE))
        {
            mRankedProteomePeptides = new RankedProteomePeptides(configBuilder.getValue(PROTEOME_RANKS_FILE));
        }
        else
        {
            mRankedProteomePeptides = null;
        }

        mScorer = new BindScorer(new ScoreConfig(configBuilder));

        mOutputDir = parseOutputDir(configBuilder);
        mOutputId = configBuilder.getValue(OUTPUT_ID);
        mThreads = parseThreads(configBuilder);
    }

    public void run()
    {
        if(mPeptideSimilarities.isEmpty())
        {
            NE_LOGGER.info("no allele peptides loaded");
            System.exit(1);
        }

        if(mRankedProteomePeptides != null && !mRankedProteomePeptides.isValid())
        {
            NE_LOGGER.info("ranked proteome peptides invalid");
            System.exit(1);
        }

        if(!mScorer.loadScoringData())
        {
            NE_LOGGER.error("failed to load scoring data");
            System.exit(1);
        }

        NE_LOGGER.info("finding top similarities for {} peptides", mPeptideSimilarities.size());

        List<PeptideSearchTask> searchTasks = Lists.newArrayList();

        if(mThreads > 1)
        {
            int threads = min(mThreads, mPeptideSimilarities.size());

            for(int i = 0; i < threads; ++i)
            {
                searchTasks.add(new PeptideSearchTask(i, mTransAminoAcidMap, mRankedProteomePeptides));
            }

            int taskIndex = 0;
            for(PeptideSimilarity peptideSim : mPeptideSimilarities)
            {
                searchTasks.get(taskIndex).getPeptides().add(peptideSim);
                ++taskIndex;

                if(taskIndex >= threads)
                    taskIndex = 0;
            }

            final List<Callable<Void>> callableList = searchTasks.stream().collect(Collectors.toList());
            TaskExecutor.executeTasks(callableList, callableList.size());
        }
        else
        {
            PeptideSearchTask searchTask = new PeptideSearchTask(0, mTransAminoAcidMap, mRankedProteomePeptides);
            searchTask.getPeptides().addAll(mPeptideSimilarities);
            searchTasks.add(searchTask);
            searchTask.run();
        }

        writeResults();

        NE_LOGGER.info("similarity search complete");
    }

    private void writeResults()
    {
        try
        {
            String filename = BindCommon.formFilename(mOutputDir, "peptide_prot_sim", mOutputId);
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("Allele,Peptide,PeptideLikelihoodRank");

            if(mTransAminoAcidMap != null)
            {
                writer.write(",TopPeptide,TopLikelihoodRank,TopDtoS,TopPosDiffs,GeneName,TransName");
            }

            if(mRankedProteomePeptides != null)
            {
                writer.write(",NearestBinderPeptide,NearestBinderLikelihoodRank,NearestBinderDtoS,NearestBinderPosDiffs");
            }

            boolean hasWildtypes = mPeptideSimilarities.stream().anyMatch(x -> !x.WildtypePeptide.isEmpty());

            if(hasWildtypes)
            {
                writer.write(",WildtypePeptide,WildtypeLikelihoodRank,WildtypeDtoS,WildtypePosDiffs");
            }

            writer.newLine();

            final BlosumMapping blosumMapping = new BlosumMapping();

            for(PeptideSimilarity peptideSim : mPeptideSimilarities)
            {
                double peptideRank = calcPeptideLikelihoodRank(peptideSim.Allele, peptideSim.Peptide);

                writer.write(String.format("%s,%s,%.6f",
                        peptideSim.Allele, peptideSim.Peptide, peptideRank));

                if(mTransAminoAcidMap != null)
                {
                    double topSimRank = calcPeptideLikelihoodRank(peptideSim.Allele, peptideSim.topPeptide());
                    String topSimPosDiff = peptideSim.positionDiffs(peptideSim.topPeptide());

                    writer.write(String.format(",%s,%.6f,%.1f,%s,%s,%s",
                            peptideSim.topPeptide(), topSimRank, peptideSim.topSimiliarity(), topSimPosDiff,
                            peptideSim.topTransData() != null ? peptideSim.topTransData().GeneName : "NONE",
                            peptideSim.topTransData() != null ? peptideSim.topTransData().TransName : "NONE"));
                }

                if(mRankedProteomePeptides != null)
                {
                    double nearestSim = blosumMapping.calcSequenceSimilarity(peptideSim.Peptide, peptideSim.nearestBinder());
                    String nearestPosDiff = peptideSim.positionDiffs(peptideSim.nearestBinder());

                    writer.write(String.format(",%s,%.6f,%.1f,%s",
                            peptideSim.nearestBinder(), peptideSim.nearestBinderLikelihoodRank(), nearestSim, nearestPosDiff));
                }

                if(hasWildtypes)
                {
                    double wildtypeSim = blosumMapping.calcSequenceSimilarity(peptideSim.Peptide, peptideSim.WildtypePeptide);
                    String wildtypePosDiff = peptideSim.positionDiffs(peptideSim.WildtypePeptide);
                    double wildtypeRank = calcPeptideLikelihoodRank(peptideSim.Allele, peptideSim.WildtypePeptide);

                    writer.write(String.format(",%s,%.6f,%.1f,%s",
                            peptideSim.WildtypePeptide, wildtypeRank, wildtypeSim, wildtypePosDiff));
                }

                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write peptide similarity results: {}", e.toString());
        }
    }

    private double calcPeptideLikelihoodRank(final String allele, final String peptide)
    {
        BindData bindData = new BindData(allele, peptide, "");
        mScorer.calcScoreData(bindData);
        return bindData.likelihoodRank();
    }

    private void loadPeptides(final String filename)
    {
        if(filename == null || !Files.exists(Paths.get(filename)))
        {
            NE_LOGGER.error("allele peptide file({}) not found", filename);
            return;
        }

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));
            String header = lines.get(0);
            lines.remove(0);

            Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, BIND_DELIM);

            // Peptide,Allele,Foreignness,WtPeptide,DisToSelf,DtsPeptide,Immunogenic,PRIME
            int alleleIndex = fieldsIndexMap.get(FLD_ALLELE);
            int peptideIndex = fieldsIndexMap.get(FLD_PEPTIDE);
            Integer wildtypeIndex = fieldsIndexMap.get("WtPeptide");

            for(String line : lines)
            {
                final String[] values = line.split(CSV_DELIM, -1);

                String allele = cleanAllele(values[alleleIndex]);
                String peptide = values[peptideIndex];

                if(!RANKED_PROTEOME_PEPTIDE_LENGTHS.contains(peptide.length()))
                    continue;

                String wtPeptide = wildtypeIndex != null ? values[wildtypeIndex] : "";

                mPeptideSimilarities.add(new PeptideSimilarity(peptide, allele, wtPeptide));
            }

            NE_LOGGER.info("loaded {} peptides from file({})", lines.size(), filename);
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to read allele peptides rank file: {}", e.toString());
        }
    }

    private class PeptideSearchTask implements Callable<Void>
    {
        private final int mTaskId;
        private final Map<String, List<TranscriptAminoAcids>> mTransAminoAcidMap;
        private final RankedProteomePeptides mRankedProteomePeptides;

        private final List<PeptideSimilarity> mPeptideSimilarities;

        private final BlosumMapping mBlosumMapping;

        public PeptideSearchTask(
                int taskId, final Map<String, List<TranscriptAminoAcids>> transAminoAcidMap,
                final RankedProteomePeptides rankedProteomePeptides)
        {
            mTaskId = taskId;
            mTransAminoAcidMap = transAminoAcidMap;
            mRankedProteomePeptides = rankedProteomePeptides;

            mPeptideSimilarities = Lists.newArrayList();

            mBlosumMapping = new BlosumMapping();
        }

        public List<PeptideSimilarity> getPeptides() { return mPeptideSimilarities; }

        @Override
        public Void call()
        {
            run();
            return null;
        }

        private void run()
        {
            NE_LOGGER.info("{}: searching for {} peptides", mTaskId, mPeptideSimilarities.size());

            for(int i = 0; i < mPeptideSimilarities.size(); ++i)
            {
                PeptideSimilarity peptideSim = mPeptideSimilarities.get(i);

                try
                {
                    findTopSimilarity(peptideSim);
                    findNearestBinder(peptideSim);
                }
                catch(Exception e)
                {
                    NE_LOGGER.error("{}: error processing peptide({}): {}",
                            mTaskId, peptideSim, e.toString());

                    e.printStackTrace();
                }

                if(i > 0 && (i % 10) == 0)
                {
                    NE_LOGGER.info("{}: search count: {}", mTaskId, i);
                }
            }

            NE_LOGGER.info("{}: peptide search complete", mTaskId);
        }

        private void findNearestBinder(final PeptideSimilarity peptideSim)
        {
            if(mRankedProteomePeptides == null)
                return;

            // search amongst the strong binders from the proteome - previously computed and cached per allele
            String peptide = peptideSim.Peptide;
            int peptideLength = peptide.length();

            String topPeptide = "";
            double topSimilarity = 0;
            double topRank = 0;

            Map<String, Double> allelePeptideRanks = mRankedProteomePeptides.getPeptideRanks(peptideSim.Allele, peptideLength);

            if(allelePeptideRanks == null)
            {
                NE_LOGGER.warn("peptide({}) missing ranking data", peptideSim);
                return;
            }

            for(Map.Entry<String, Double> entry : allelePeptideRanks.entrySet())
            {
                String otherPeptide = entry.getKey();

                if(otherPeptide.equals(peptide))
                {
                    peptideSim.setTopBinder(peptide, entry.getValue());
                    return;
                }

                double similarity = 0;
                boolean skip = false;

                for(int i = 0; i < peptide.length(); ++i)
                {
                    char aa1 = peptide.charAt(i);
                    char aa2 = otherPeptide.charAt(i);

                    int bs1 = mBlosumMapping.selfMapping(aa1);
                    int bs2 = mBlosumMapping.selfMapping(aa2);
                    int map = mBlosumMapping.map(aa1, aa2);

                    similarity += (bs1 + bs2) * 0.5 - map;

                    if(!topPeptide.isEmpty() && similarity >= topSimilarity) // early exit if cannot be better
                    {
                        skip = true;
                        break;
                    }
                }

                if(skip)
                    continue;

                if(topPeptide.isEmpty() || similarity < topSimilarity)
                {
                    topPeptide = otherPeptide;
                    topSimilarity = similarity;
                    topRank = entry.getValue();
                }
            }

            peptideSim.setTopBinder(topPeptide, topRank);
        }

        private void findTopSimilarity(final PeptideSimilarity peptideSim)
        {
            if(mTransAminoAcidMap == null)
                return;

            // search the proteome for the nearest peptide
            String topPeptide = "";
            double topSimiliarity = 0;
            TranscriptAminoAcids topTrans = null;

            String peptide = peptideSim.Peptide;

            int peptideLength = peptide.length();

            for(List<TranscriptAminoAcids> transAaList : mTransAminoAcidMap.values())
            {
                Set<String> uniquePeptides = Sets.newHashSet();

                for(TranscriptAminoAcids transAminoAcids : transAaList)
                {
                    final String aminoAcids = transAminoAcids.AminoAcids;

                    if(aminoAcids.contains(peptide))
                    {
                        peptideSim.setTopSimilarity(peptide, 0, transAminoAcids);
                        return;
                    }

                    int aaLength = aminoAcids.length();

                    for(int startIndex = 0; startIndex < aaLength - peptideLength; ++startIndex)
                    {
                        String aaPeptide = aminoAcids.substring(startIndex, startIndex + peptideLength);

                        if(uniquePeptides.contains(aaPeptide))
                            continue;

                        uniquePeptides.add(aaPeptide);

                        if(aaPeptide.contains(AMINO_ACID_21ST))
                            continue;

                        double similarity = 0;
                        boolean skip = false;

                        for(int i = 0; i < peptide.length(); ++i)
                        {
                            char aa1 = peptide.charAt(i);
                            char aa2 = aaPeptide.charAt(i);

                            int bs1 = mBlosumMapping.selfMapping(aa1);
                            int bs2 = mBlosumMapping.selfMapping(aa2);
                            int map = mBlosumMapping.map(aa1, aa2);

                            similarity += (bs1 + bs2) * 0.5 - map;

                            if(topTrans != null && similarity >= topSimiliarity) // skip if cannot be better
                            {
                                skip = true;
                                break;
                            }
                        }

                        if(skip)
                            continue;

                        if(topTrans == null || similarity < topSimiliarity)
                        {
                            topPeptide = aaPeptide;
                            topSimiliarity = similarity;
                            topTrans = transAminoAcids;
                        }
                    }
                }
            }

            peptideSim.setTopSimilarity(topPeptide, topSimiliarity, topTrans);
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        addEnsemblDir(configBuilder);
        configBuilder.addPath(PEPTIDES_FILE, true, "Peptides file");
        configBuilder.addPath(PROTEOME_RANKS_FILE, true, "Proteome ranks file");
        ScoreConfig.registerConfig(configBuilder);
        addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
        addOutputOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        PeptideProteomeSimilarity peptideProteomeSimilarity = new PeptideProteomeSimilarity(configBuilder);
        peptideProteomeSimilarity.run();
    }
}
