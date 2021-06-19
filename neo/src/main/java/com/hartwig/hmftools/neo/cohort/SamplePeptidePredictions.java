package com.hartwig.hmftools.neo.cohort;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;

import static com.hartwig.hmftools.common.sigs.DataUtils.convertList;
import static com.hartwig.hmftools.common.stats.Percentiles.PERCENTILE_COUNT;
import static com.hartwig.hmftools.common.stats.Percentiles.buildPercentiles;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.NeoCommon.LOG_DEBUG;
import static com.hartwig.hmftools.neo.NeoCommon.loadSampleIdsFile;
import static com.hartwig.hmftools.neo.cohort.AlleleCoverage.EXPECTED_ALLELE_COUNT;
import static com.hartwig.hmftools.neo.cohort.AlleleCoverage.getGeneStatus;
import static com.hartwig.hmftools.neo.cohort.DataLoader.loadAlleleCoverage;
import static com.hartwig.hmftools.neo.cohort.DataLoader.loadPredictionData;
import static com.hartwig.hmftools.neo.cohort.PredictionData.expandHomozygous;
import static com.hartwig.hmftools.neo.cohort.StatusResults.NORMAL;
import static com.hartwig.hmftools.neo.cohort.StatusResults.SIM_TUMOR;
import static com.hartwig.hmftools.neo.cohort.StatusResults.STATUS_MAX;
import static com.hartwig.hmftools.neo.cohort.StatusResults.TUMOR;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class SamplePeptidePredictions
{
    private final String mOutputDir;
    private final String mNeoDataDir;
    private final String mLilacDataDir;
    private final String mPredictionsDataDir;
    private final List<String> mSampleIds;

    private final Map<String,AllelePredictions> mAllelePredictions;

    private BufferedWriter mSummaryWriter;
    private BufferedWriter mPeptideWriter;

    private final double mAffinitySumFactor;
    private final int mAffinityLowCountThreshold;
    private final int mAffinityMediumCountThreshold;
    private final double mWritePeptideAffinityThreshold;

    public static final String SAMPLE_ID_FILE = "sample_id_file";
    public static final String NEO_DATA_DIR = "neo_data_dir";
    public static final String LILAC_DATA_DIR = "lilac_data_dir";
    public static final String PREDICTION_DATA_DIR = "prediction_data_dir";

    private static final String AFF_SUM_FACTOR = "sum_factor";
    private static final String AFF_LOW_THRESHOLD = "affinity_low_threshold";
    private static final String AFF_MED_THRESHOLD = "affinity_med_threshold";

    private static final String WRITE_PEPTIDES = "write_peptides";
    private static final String WRITE_PEPTIDE_AFF_THRESHOLD = "write_pep_aff_threshold";

    // constants
    private static final double PREDICTION_FACTOR = 2;

    public SamplePeptidePredictions(final CommandLine cmd)
    {
        mSampleIds = Lists.newArrayList();
        loadSampleIdsFile(cmd.getOptionValue(SAMPLE_ID_FILE), mSampleIds);

        mNeoDataDir = cmd.getOptionValue(NEO_DATA_DIR);
        mPredictionsDataDir = cmd.getOptionValue(PREDICTION_DATA_DIR);
        mLilacDataDir = cmd.getOptionValue(LILAC_DATA_DIR);

        mAffinityLowCountThreshold = Integer.parseInt(cmd.getOptionValue(AFF_LOW_THRESHOLD, "25"));
        mAffinityMediumCountThreshold = Integer.parseInt(cmd.getOptionValue(AFF_MED_THRESHOLD, "50"));
        mAffinitySumFactor = Double.parseDouble(cmd.getOptionValue(AFF_SUM_FACTOR, "2"));
        mWritePeptideAffinityThreshold = Double.parseDouble(cmd.getOptionValue(WRITE_PEPTIDE_AFF_THRESHOLD, "0"));

        mAllelePredictions = Maps.newHashMap();

        mOutputDir = parseOutputDir(cmd);
        mSummaryWriter = null;
        mPeptideWriter = null;

        initialiseSampleWriter();

        if(cmd.hasOption(WRITE_PEPTIDES))
            initialisePeptideWriter();
    }

    public void run()
    {
        if(mSampleIds.isEmpty())
            return;

        initialiseSampleWriter();

        NE_LOGGER.info("processing {} samples", mSampleIds.size());

        // check required inputs and config
        int processed = 0;

        for(final String sampleId : mSampleIds)
        {
            processSample(sampleId);
            ++processed;

            if(processed > 0 && (processed % 100) == 0)
            {
                NE_LOGGER.info("processed {} samples", processed);
            }
        }

        writeAllelePredictions();

        closeBufferedWriter(mSummaryWriter);
        closeBufferedWriter(mPeptideWriter);
    }

    private void processSample(final String sampleId)
    {
        List<AlleleCoverage> alleleCoverages = loadAlleleCoverage(sampleId, mLilacDataDir);

        List<Boolean> geneLostStatus = getGeneStatus(alleleCoverages);

        List<PredictionData> allPredictions = loadPredictionData(sampleId, mPredictionsDataDir);

        Map<String,PeptideScores> peptideScores = Maps.newHashMap();

        Map<String,List<PredictionData>> peptidePredictions = Maps.newHashMap();

        // organise into map by peptide, avoiding repeated peptides
        for(int i = 0; i < allPredictions.size(); ++i)
        {
            PredictionData predData = allPredictions.get(i);

            addAllelePrediction(predData);

            List<PredictionData> predictions = peptidePredictions.get(predData.Peptide);

            if(predictions == null)
            {
                peptidePredictions.put(predData.Peptide, Lists.newArrayList(predData));
            }
            else
            {
                if(predictions.stream().noneMatch(x -> x.Allele.equals(predData.Allele)))
                {
                    predictions.add(predData);
                }
            }
        }

        boolean allValid = true;
        for(Map.Entry<String,List<PredictionData>> entry : peptidePredictions.entrySet())
        {
            String peptide = entry.getKey();
            List<PredictionData> predictions = entry.getValue();

            expandHomozygous(predictions);

            if(predictions.size() != EXPECTED_ALLELE_COUNT)
            {
                NE_LOGGER.error("peptide({}) has incorrect allele prediction count({})", peptide, predictions.size());
                continue;
            }

            // process the 6 alleles using the coverage
            double maxAffinity = predictions.stream().mapToDouble(x -> x.Affinity).max().orElse(0);
            double minPresentation = predictions.stream().mapToDouble(x -> x.Presentation).min().orElse(0);

            PeptideScores scores = new PeptideScores(peptide, maxAffinity, minPresentation);

            for(int alleleIndex = 0; alleleIndex < alleleCoverages.size(); ++alleleIndex)
            {
                AlleleCoverage alleleCoverage = alleleCoverages.get(alleleIndex);

                PredictionData allelePrediction = predictions.stream()
                        .filter(x -> x.Allele.equals(alleleCoverage.Allele)).findFirst().orElse(null);

                if(allelePrediction == null)
                {
                    NE_LOGGER.error("peptide({}) missing allele({}) prediction", peptide, alleleCoverage.Allele);
                    allValid = false;
                    break;
                }

                scores.Affinity[NORMAL] = min(scores.Affinity[NORMAL], allelePrediction.Affinity);
                scores.Presentation[NORMAL] = max(scores.Presentation[NORMAL], allelePrediction.Presentation);

                if(!alleleCoverage.isLost())
                {
                    scores.Affinity[TUMOR] = min(scores.Affinity[TUMOR], allelePrediction.Affinity);
                    scores.Presentation[TUMOR] = max(scores.Presentation[TUMOR], allelePrediction.Presentation);
                }

                if(alleleCoverage.isLost() || !geneLostStatus.get(alleleIndex))
                {
                    scores.Affinity[SIM_TUMOR] = min(scores.Affinity[SIM_TUMOR], allelePrediction.Affinity);
                    scores.Presentation[SIM_TUMOR] = max(scores.Presentation[SIM_TUMOR], allelePrediction.Presentation);
                }
            }

            if(!allValid)
                break;

            peptideScores.put(peptide, scores);
        }

        if(!allValid)
        {
            NE_LOGGER.warn("sampleId({}) has consistent allele coverage vs prediction alleles", sampleId);
            return;
        }

        SampleSummary sampleSummary = new SampleSummary(peptideScores.size());

        for(Map.Entry<String,PeptideScores> entry : peptideScores.entrySet())
        {
            String peptide = entry.getKey();
            PeptideScores scores = entry.getValue();

            for(int i = 0; i < STATUS_MAX; ++i)
            {
                sampleSummary.Results[i].AffinityTotal += pow(1 / scores.Affinity[i], mAffinitySumFactor);

                if(scores.Affinity[i] <= mAffinityLowCountThreshold)
                    ++sampleSummary.Results[i].AffinityLowCount;

                if(scores.Affinity[i] <= mAffinityMediumCountThreshold)
                    ++sampleSummary.Results[i].AffinityMediumCount;

                sampleSummary.Results[i].PresentationTotal += pow(scores.Presentation[i], mAffinitySumFactor);

                if(scores.Presentation[i] >= 0.95)
                    ++sampleSummary.Results[i].PresentationCount;
            }

            writePeptideData(sampleId, peptide, scores);
        }

        writeSampleSummary(sampleId, sampleSummary);
    }

    private void addAllelePrediction(PredictionData predData)
    {
        AllelePredictions predictions = mAllelePredictions.get(predData.Allele);

        if(predictions == null)
        {
            predictions = new AllelePredictions(predData.Allele);
            mAllelePredictions.put(predData.Allele, predictions);
        }

        predictions.addPeptide(predData.Peptide, predData.Affinity);
    }

    private void initialiseSampleWriter()
    {
        try
        {
            final String outputFileName = mOutputDir + "NEO_SAMPLE_SUMMARY.csv";

            mSummaryWriter = createBufferedWriter(outputFileName, false);
            mSummaryWriter.write("SampleId,PeptideCount," + SampleSummary.header());
            mSummaryWriter.newLine();
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to create sample summary writer: {}", e.toString());
        }
    }

    private void writeSampleSummary(final String sampleId, final SampleSummary sampleSummary)
    {
        try
        {
            mSummaryWriter.write(String.format("%s,%d",
                    sampleId, sampleSummary.PeptideCount));

            for(int i = 0; i < STATUS_MAX; ++i)
            {
                final StatusResults results = sampleSummary.Results[i];

                mSummaryWriter.write(String.format(",%.4g,%d,%d,%4g,%d",
                        results.AffinityTotal, results.AffinityLowCount, results.AffinityMediumCount,
                        results.PresentationTotal, results.PresentationCount));
            }

            mSummaryWriter.newLine();
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to write neo-epitope data: {}", e.toString());
        }
    }

    private void writeAllelePredictions()
    {
        try
        {
            /*
            final String outputFileName = mOutputDir + "NEO_ALLELE_PERCENTILES.csv";

            BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write("Allele");

            for(int i = 0; i < PERCENTILE_COUNT; ++i)
            {
                writer.write(String.format(",Pct_%.2f", i * 0.01));
            }

            writer.newLine();

            for(Map.Entry<String,AllelePredictions> entry : mAllelePredictions.entrySet())
            {
                String allele = entry.getKey();

                writer.write(allele);

                AllelePredictions predictions = entry.getValue();
                int[] percentiles = predictions.buildPercentiles();

                for(int i = 0; i < percentiles.length; ++i)
                {
                    writer.write(String.format(",%d", percentiles[i]));
                }

                writer.newLine();
            }

            */

            final String outputFileName = mOutputDir + "NEO_ALLELE_FREQUENCIES.csv";

            BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write("Allele,Total");

            int maxBucket = 500;

            for(int i = 0; i < maxBucket; ++i)
            {
                writer.write(String.format(",Fq_%d", i));
            }

            writer.newLine();

            for(Map.Entry<String,AllelePredictions> entry : mAllelePredictions.entrySet())
            {
                String allele = entry.getKey();
                AllelePredictions predictions = entry.getValue();

                writer.write(String.format("%s,%d", allele, predictions.getTotal()));

                int[] affinityFrequencies = predictions.getFrequencies();

                for(int i = 0; i < maxBucket; ++i)
                {
                    writer.write(String.format(",%d", affinityFrequencies[i]));
                }

                writer.newLine();
            }

            writer.close();
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to write allele percentiles: {}", e.toString());
        }
    }

    private void initialisePeptideWriter()
    {
        try
        {
            final String outputFileName = mOutputDir + "NEO_PEPTIDE_SCORES.csv";

            mPeptideWriter = createBufferedWriter(outputFileName, false);
            mPeptideWriter.write("SampleId,Peptide");
            mPeptideWriter.write(",NormalAffinity,TumorAffinity,SimTumorAffinity");
            mPeptideWriter.write(",NormalPresentation,TumorPresentation,SimTumorPresentation");
            mPeptideWriter.newLine();
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to create peptide writer: {}", e.toString());
        }
    }

    private void writePeptideData(final String sampleId, final String peptide, final PeptideScores scores)
    {
        if(mPeptideWriter == null)
            return;

        if(mWritePeptideAffinityThreshold > 0 &&  Arrays.stream(scores.Affinity).noneMatch(x -> x < mWritePeptideAffinityThreshold))
            return;

        try
        {
            mPeptideWriter.write(String.format("%s,%s", sampleId, peptide));

            for(int i = 0; i < STATUS_MAX; ++i)
            {
                mPeptideWriter.write(String.format(",%.1f", scores.Affinity[i]));
            }

            for(int i = 0; i < STATUS_MAX; ++i)
            {
                mPeptideWriter.write(String.format(",%.4f", scores.Presentation[i]));
            }

            mPeptideWriter.newLine();
        }
        catch (IOException e)
        {
            NE_LOGGER.error("failed to write neo-epitope data: {}", e.toString());
        }
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(SAMPLE_ID_FILE, true, "SampleId file");
        options.addOption(NEO_DATA_DIR, true, "Directory for sample neo-epitope files");
        options.addOption(PREDICTION_DATA_DIR, true, "Directory for sample prediction result files");
        options.addOption(LILAC_DATA_DIR, true, "Directory for Lilac coverage files");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(LOG_DEBUG, false, "Log verbose");
        options.addOption(WRITE_PEPTIDES, false, "Write all peptide scores");

        options.addOption(AFF_SUM_FACTOR, true, "Affinity sum factor");
        options.addOption(AFF_LOW_THRESHOLD, true, "Affinity low count threshold");
        options.addOption(AFF_MED_THRESHOLD, true, "Affinity median count threshold");
        options.addOption(WRITE_PEPTIDE_AFF_THRESHOLD, true, "Only write peptides with affinity less than this if > 0");
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();

        SamplePeptidePredictions.addCmdLineArgs(options);

        final CommandLine cmd = createCommandLine(args, options);

        if(cmd.hasOption(LOG_DEBUG))
            Configurator.setRootLevel(Level.DEBUG);

        SamplePeptidePredictions samplePeptidePredictions = new SamplePeptidePredictions(cmd);
        samplePeptidePredictions.run();

        NE_LOGGER.info("cohort peptide predictions complete");
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
