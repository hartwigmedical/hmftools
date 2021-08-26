package com.hartwig.hmftools.lilac.cohort;

import static java.lang.Math.log10;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.DELIM;
import static com.hartwig.hmftools.lilac.cohort.CohortCommon.SAMPLE_FILES_DIR;
import static com.hartwig.hmftools.lilac.cohort.CohortCommon.SAMPLE_IDS_FILE;
import static com.hartwig.hmftools.lilac.cohort.CohortCommon.loadSampleIds;
import static com.hartwig.hmftools.lilac.cohort.CohortFrequency.ALLELE_FREQUENCY_FILE;
import static com.hartwig.hmftools.lilac.coverage.HlaComplexFile.parseCandidateCoverageData;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class CandidateScores
{
    private final CohortFrequency mCohortFrequency;
    private final SampleTruthSet mTruthSet;
    private final String mOutputDir;
    private final String mCandidateFilesDir;
    private final List<String> mSampleIds;

    private BufferedWriter mWriter;

    public CandidateScores(final CommandLine cmd)
    {
        mCohortFrequency = new CohortFrequency(cmd.getOptionValue(ALLELE_FREQUENCY_FILE));
        mTruthSet = new SampleTruthSet(cmd);
        mOutputDir = parseOutputDir(cmd);
        mCandidateFilesDir = checkAddDirSeparator(cmd.getOptionValue(SAMPLE_FILES_DIR));

        mSampleIds = loadSampleIds(cmd.getOptionValue(SAMPLE_IDS_FILE));
    }

    public void run()
    {
        if(mSampleIds.isEmpty())
        {
            LL_LOGGER.error("sample IDs loading failed");
            System.exit(1);
        }

        if(mCohortFrequency.getAlleleFrequencies().isEmpty() || mTruthSet.getSampleAlleleSet().isEmpty())
        {
            LL_LOGGER.error("reference data loading failed");
            System.exit(1);
        }

        if(mCandidateFilesDir == null || mOutputDir == null)
        {
            LL_LOGGER.error("invalid input/output paths");
            System.exit(1);
        }

        initialiseWriter();

        int sampleCount = 0;
        for(String sampleId : mSampleIds)
        {
            List<CandidateAlleles> candidates = loadSampleCandidates(sampleId);
            processSampleCandidates(sampleId, candidates);

            ++sampleCount;

            if((sampleCount % 100 == 0))
            {
                LL_LOGGER.info("processed {} samples", sampleCount);
            }
        }

        closeBufferedWriter(mWriter);
    }

    private void processSampleCandidates(final String sampleId, final List<CandidateAlleles> candidates)
    {
        final List<HlaAllele> truthSetAlleles = mTruthSet.getSampleAlleles(sampleId);
        double topScore = 0;

        for(CandidateAlleles candidate : candidates)
        {
            double cohortFrequencyTotal = 0;

            for(HlaAllele allele : candidate.Alleles)
            {
                double cohortFrequency = mCohortFrequency.getAlleleFrequency(allele);
                cohortFrequencyTotal += log10(max(cohortFrequency,0.0001));
            }

            candidate.setCohortFrequencyTotal(cohortFrequencyTotal);
            topScore = max(topScore, candidate.calcScore());
        }

        double maxTopScore = topScore;
        candidates.forEach(x -> writeCandidateData(sampleId, x, maxTopScore, truthSetAlleles));
    }

    private List<CandidateAlleles> loadSampleCandidates(final String sampleId)
    {
        final String filename = mCandidateFilesDir + sampleId + ".candidates.coverage.txt";
        List<CandidateAlleles> candidates = Lists.newArrayList();

        try
        {
            if(!Files.exists(Paths.get(filename)))
            {
                LL_LOGGER.warn("sample({}) no candidates file found", sampleId);
                return candidates;
            }

            final List<String> lines = Files.readAllLines(new File(filename).toPath());

            final Map<String,Integer> fieldsMap = createFieldsIndexMap(lines.get(0), DELIM);
            lines.remove(0);

            int totalCoverageIndex = fieldsMap.get("totalCoverage");
            int firstAlleleIndex = fieldsMap.get("allele1");

            for(String line : lines)
            {
                String[] items = line.split(DELIM);
                int totalCoverage = Integer.parseInt(items[totalCoverageIndex]);

                List<String> rawAlleleData = Lists.newArrayList();

                for(int index = firstAlleleIndex; index < min(firstAlleleIndex + 6, items.length); ++index)
                {
                    rawAlleleData.add(items[index]);
                }

                List<HlaAllele> allAlleles = parseCandidateCoverageData(rawAlleleData);

                if(allAlleles.size() != 6)
                {
                    LL_LOGGER.error("sample({}) invalid candidates: {}", sampleId, line);
                    continue;
                }

                candidates.add(new CandidateAlleles(candidates.size(), allAlleles, totalCoverage));
            }

            LL_LOGGER.debug("sample({}) load {} candidates", mSampleIds.size());
        }
        catch (IOException e)
        {
            LL_LOGGER.warn("failed to load sample candidates coverage file({}): {}", filename, e.toString());
        }

        return candidates;
    }

    private void initialiseWriter()
    {
        try
        {
            final String outputFileName = mOutputDir + "lilac_sample_candidate_scores.csv";
            mWriter = createBufferedWriter(outputFileName, false);

            mWriter.write("SampleId,");
            mWriter.write(CandidateAlleles.header());
            mWriter.newLine();

        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to write sample candidates file: {}", e.toString());
        }
    }

    private void writeCandidateData(
            final String sampleId, final CandidateAlleles candidate, double topScore, final List<HlaAllele> truthSetAlleles)
    {
        try
        {
            mWriter.write(String.format("%s,%s", sampleId, candidate.toCsv(topScore, truthSetAlleles)));
            mWriter.newLine();

        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to write sample candidates file: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        LL_LOGGER.info("processing cohort sample candidate files");

        Options options = new Options();
        options.addOption(SAMPLE_IDS_FILE, true, "Sample IDs");
        options.addOption(SAMPLE_FILES_DIR, true, "Path to candidate-coverage files");
        addOutputDir(options);
        SampleTruthSet.addCmdLineOptions(options);
        CohortFrequency.addCmdLineOptions(options);

        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args);

        CandidateScores candidateScores = new CandidateScores(cmd);
        candidateScores.run();

        LL_LOGGER.info("cohort processing complete");
    }

}
