package com.hartwig.hmftools.lilac.cohort;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.DELIM;
import static com.hartwig.hmftools.lilac.cohort.CohortCommon.SAMPLE_FILES_DIR;
import static com.hartwig.hmftools.lilac.cohort.CohortCommon.SAMPLE_IDS_FILE;
import static com.hartwig.hmftools.lilac.cohort.CohortCommon.loadSampleIds;
import static com.hartwig.hmftools.lilac.coverage.HlaComplexFile.parseCandidateCoverageData;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class SampleSummaries
{
    private final String mOutputDir;
    private final String mCandidateFilesDir;
    private final List<String> mSampleIds;

    private BufferedWriter mWriter;

    public SampleSummaries(final CommandLine cmd)
    {
        mOutputDir = parseOutputDir(cmd);
        mCandidateFilesDir = checkAddDirSeparator(cmd.getOptionValue(SAMPLE_FILES_DIR));

        mSampleIds = loadSampleIds(cmd.getOptionValue(SAMPLE_IDS_FILE));
    }

    public void run()
    {
        if(mSampleIds.isEmpty())
        {
            LL_LOGGER.error("sampleIds loading failed");
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
            processSampleCandidates(sampleId);

            ++sampleCount;

            if((sampleCount % 100 == 0))
            {
                LL_LOGGER.info("processed {} samples", sampleCount);
            }
        }

        closeBufferedWriter(mWriter);
    }

    private void processSampleCandidates(final String sampleId)
    {
    }

    /*
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
    */

    private void initialiseWriter()
    {
        try
        {
            final String outputFileName = mOutputDir + "LILAC_COHORT.csv";
            mWriter = createBufferedWriter(outputFileName, false);

            mWriter.write("SampleId,");
            // mWriter.write(CandidateAlleles.header());
            mWriter.newLine();

        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to write sample candidates file: {}", e.toString());
        }
    }

    private void writeSampleData(final String sampleId)
    {
        try
        {
            // mWriter.write(String.format("%s,%s", sampleId, candidate.toCsv(topScore, truthSetAlleles)));
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
        options.addOption(OUTPUT_DIR, true, "Path to output");
        SampleTruthSet.addCmdLineOptions(options);
        CohortFrequency.addCmdLineOptions(options);

        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args);

        SampleSummaries candidateScores = new SampleSummaries(cmd);
        candidateScores.run();

        LL_LOGGER.info("cohort processing complete");
    }

}
