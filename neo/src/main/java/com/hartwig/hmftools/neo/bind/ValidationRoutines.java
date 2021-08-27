package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindConstants.INVALID_SCORE;
import static com.hartwig.hmftools.neo.bind.BindData.loadBindData;
import static com.hartwig.hmftools.neo.bind.HlaSequences.HLA_DEFINITIONS_FILE;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.MatrixUtils;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class ValidationRoutines
{
    private final BinderConfig mConfig;

    private final Map<String,Map<Integer,List<BindData>>> mAlleleTrainingData;
    private final Map<String,Map<Integer,List<BindData>>> mAlleleValidationData;
    private final List<String> mValidationAlleles;

    private final Map<String,Map<Integer,BindCountData>> mAlleleBindCounts; // counts data by peptide length>
    private final PosWeightModel mPosWeightModel;
    private final HlaSequences mHlaSequences;
    private final RandomPeptideDistribution mRandomDistribution;
    private final Map<String,Integer> mAlleleTotalCounts;

    private final BufferedWriter mPeptideWriter;

    public ValidationRoutines(final CommandLine cmd)
    {
        mConfig = new BinderConfig(cmd);
        mAlleleTrainingData = Maps.newHashMap();
        mAlleleValidationData = Maps.newHashMap();
        mAlleleBindCounts = Maps.newHashMap();
        mAlleleTotalCounts = Maps.newHashMap();
        mValidationAlleles = Lists.newArrayList();

        mHlaSequences = new HlaSequences();
        mHlaSequences.load(cmd.getOptionValue(HLA_DEFINITIONS_FILE));

        mPosWeightModel = new PosWeightModel(mConfig.Constants, mHlaSequences);

        mRandomDistribution = new RandomPeptideDistribution(mConfig.RandomPeptides);

        mPeptideWriter = initialisePeptideWriter();
    }

    public void run()
    {
        if(!loadData())
        {
            NE_LOGGER.error("failed to load data");
            return;
        }

        mAlleleValidationData.keySet().forEach(x -> mValidationAlleles.add(x));

        NE_LOGGER.info("running validation for {} alleles", mValidationAlleles.size());

        buildBindCounts();

        for(String allele : mValidationAlleles)
        {
            // if allele isn't in the training set then no need for evaluation
            if(!mAlleleTrainingData.containsKey(allele))
                continue;

            int trainingDataCount = mAlleleTrainingData.get(allele).values().stream().mapToInt(x -> x.size()).sum();
            NE_LOGGER.info("scoring without allele({}) and {} training data items", allele, trainingDataCount);
            testWithExcludedAllele(allele);
        }

        NE_LOGGER.info("validation routine complete");
    }

    private void buildBindCounts()
    {
        // translate binds to counts and weighted counts - ie the step prior to allele-blending
        NE_LOGGER.info("building bind counts", mValidationAlleles.size());

        for(Map.Entry<String, Map<Integer, List<BindData>>> alleleEntry : mAlleleTrainingData.entrySet())
        {
            final String allele = alleleEntry.getKey();

            final Map<Integer, List<BindData>> pepLenBindDataMap = alleleEntry.getValue();

            Map<Integer, BindCountData> pepLenCountsMap = Maps.newHashMap(); // counts data by peptide length
            mAlleleBindCounts.put(allele, pepLenCountsMap);

            for(Map.Entry<Integer, List<BindData>> pepLenEntry : pepLenBindDataMap.entrySet())
            {
                List<BindData> bindDataList = pepLenEntry.getValue();

                for(BindData bindData : bindDataList)
                {
                    int peptideLength = bindData.peptideLength();

                    BindCountData bindCounts = pepLenCountsMap.get(bindData.peptideLength());

                    if(bindCounts == null)
                    {
                        bindCounts = new BindCountData(allele, peptideLength);
                        pepLenCountsMap.put(peptideLength, bindCounts);
                    }

                    bindCounts.processBindData(bindData, false);
                }
            }
        }

        // fill in an gaps in alleles or peptide lengths if they are required
        if(!mValidationAlleles.isEmpty())
        {
            mValidationAlleles.stream()
                    .filter(x -> !mAlleleBindCounts.containsKey(x))
                    .forEach(x -> mAlleleBindCounts.put(x, Maps.newHashMap()));
        }

        for(Map.Entry<String, Map<Integer, BindCountData>> alleleEntry : mAlleleBindCounts.entrySet())
        {
            final String allele = alleleEntry.getKey();
            final Map<Integer, BindCountData> pepLenBindCountsMap = alleleEntry.getValue();

            if(!mConfig.RequiredPeptideLengths.isEmpty())
            {
                mConfig.RequiredPeptideLengths.stream()
                        .filter(x -> !pepLenBindCountsMap.containsKey(x))
                        .forEach(x -> pepLenBindCountsMap.put(x, new BindCountData(allele, x)));
            }
        }

        for(Map.Entry<String,Map<Integer,BindCountData>> alleleEntry : mAlleleBindCounts.entrySet())
        {
            final String allele = alleleEntry.getKey();

            final Map<Integer,BindCountData> pepLenBindCountsMap = alleleEntry.getValue();

            final List<BindCountData> pepLenBindCounts = pepLenBindCountsMap.values().stream().collect(Collectors.toList());

            int totalAlelleCount = 0;

            for(BindCountData bindCounts : pepLenBindCounts)
            {
                totalAlelleCount += bindCounts.totalBindCount();
                mPosWeightModel.buildWeightedCounts(bindCounts, pepLenBindCounts);
            }

            mAlleleTotalCounts.put(allele, totalAlelleCount);
        }
    }

    private void testWithExcludedAllele(final String targetAllele)
    {
        final Map<String,Map<Integer,BindScoreMatrix>> alleleBindMatrices = Maps.newHashMap();

        buildPositionWeightMatrixData(targetAllele, alleleBindMatrices);

        mRandomDistribution.buildDistribution(alleleBindMatrices);

        runScoring(targetAllele, alleleBindMatrices);
    }

    private void buildPositionWeightMatrixData(final String targetAllele, final Map<String,Map<Integer,BindScoreMatrix>> alleleBindMatrices)
    {
        Map<Integer,List<BindCountData>> countsByLength = Maps.newHashMap();
        mConfig.RequiredPeptideLengths.forEach(x -> countsByLength.put(x, Lists.newArrayList()));

        Map<Integer,BindCountData> targetAllelePepLenMap = mAlleleBindCounts.get(targetAllele);

        // now factor each allele's weighted counts into all the others
        for(Integer peptideLength : mConfig.RequiredPeptideLengths)
        {
            List<BindCountData> allBindCounts = Lists.newArrayList();
            countsByLength.put(peptideLength, allBindCounts);

            for(Map<Integer,BindCountData> pepLenMap : mAlleleBindCounts.values())
            {
                BindCountData bindCounts = pepLenMap.get(peptideLength);

                // exclude the target allele from contributing its counts
                if(bindCounts.Allele.equals(targetAllele))
                    continue;

                allBindCounts.add(bindCounts);
            }

            BindCountData bindCounts = targetAllelePepLenMap.get(peptideLength);
            
            // clear any previously set values
            MatrixUtils.clear(bindCounts.getFinalWeightedCounts());
            
            mPosWeightModel.buildFinalWeightedCounts(bindCounts, allBindCounts);
        }

        Map<Integer,BindCountData> pepLenBindCounts = mAlleleBindCounts.get(targetAllele);
        Map<Integer,BindScoreMatrix> peptideLengthMatrixMap = Maps.newHashMap();
        alleleBindMatrices.put(targetAllele, peptideLengthMatrixMap);

        for(BindCountData bindCounts : pepLenBindCounts.values())
        {
            BindScoreMatrix matrix = mPosWeightModel.createMatrix(bindCounts);
            peptideLengthMatrixMap.put(matrix.PeptideLength, matrix);
        }
    }

    private void runScoring(final String targetAllele, final Map<String,Map<Integer,BindScoreMatrix>> alleleBindMatrices)
    {
        final Map<Integer,List<BindData>> pepLenBindDataMap = mAlleleValidationData.get(targetAllele);

        Map<Integer,BindScoreMatrix> pepLenMatrixMap = alleleBindMatrices.get(targetAllele);

        if(pepLenMatrixMap == null)
        {
            NE_LOGGER.warn("allele({}) has no validation data", targetAllele);
            return;
        }

        for(Map.Entry<Integer,List<BindData>> pepLenEntry : pepLenBindDataMap.entrySet())
        {
            final List<BindData> bindDataList = pepLenEntry.getValue();
            for(BindData bindData : bindDataList)
            {
                BindScoreMatrix matrix = pepLenMatrixMap.get(bindData.peptideLength());

                double score = matrix.calcScore(bindData.Peptide);
                double rankPercentile = mRandomDistribution.getScoreRank(bindData.Allele, bindData.peptideLength(), score);
                bindData.setScoreData(score, rankPercentile, INVALID_SCORE, INVALID_SCORE);

                writePeptideResults(targetAllele, bindData);
            }
        }
    }

    private BufferedWriter initialisePeptideWriter()
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(mConfig.formOutputFilename("validation_scores"), false);
            writer.write("ExcludedAllele,Allele,Peptide,Score,Rank");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to initialise peptide scores file: {}", e.toString());
            return null;
        }
    }

    private void writePeptideResults(final String excludedAllele, final BindData bindData)
    {
        try
        {
            mPeptideWriter.write(String.format("%s,%s,%s,%.4f,%.6f",
                    excludedAllele, bindData.Allele, bindData.Peptide, bindData.score(), bindData.rankPercentile()));
            mPeptideWriter.newLine();
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write peptide scores file: {}", e.toString());
        }
    }

    private boolean loadData()
    {
        if(!loadBindData(
                mConfig.TrainingDataFile, true, Lists.newArrayList(), mConfig.RequiredPeptideLengths, mAlleleTrainingData))
        {
            return false;
        }

        if(!loadBindData(
                mConfig.ValidationDataFile, true, Lists.newArrayList(), mConfig.RequiredPeptideLengths, mAlleleValidationData))
        {
            return false;
        }

        return true;
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();

        BinderConfig.addCmdLineArgs(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        ValidationRoutines validation = new ValidationRoutines(cmd);
        validation.run();
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
