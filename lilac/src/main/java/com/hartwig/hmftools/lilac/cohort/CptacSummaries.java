package com.hartwig.hmftools.lilac.cohort;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public class CptacSummaries
{
    private final String mOutputDir;
    private final String mSourceDir;
    private final List<String> mSampleIds;

    private BufferedWriter mWriter;

    private static final String CPTAC_DIR = "cptac_files_dir";
    private static final String CPTAC_SAMPLES = "samples_file";

    public CptacSummaries(final CommandLine cmd)
    {
        mOutputDir = parseOutputDir(cmd);
        mSourceDir = checkAddDirSeparator(cmd.getOptionValue(CPTAC_DIR));

        mSampleIds = Lists.newArrayList();

        try
        {
            Files.readAllLines(Paths.get(cmd.getOptionValue(CPTAC_SAMPLES))).stream()
                    .filter(x -> !x.equals("SampleId")).forEach(x -> mSampleIds.add(x));
        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to load samples: {}", e.toString());
            System.exit(1);
        }
    }

    public void run()
    {
        if(mSourceDir == null || mOutputDir == null || mSampleIds.isEmpty())
        {
            LL_LOGGER.error("invalid input/output paths or samples list");
            System.exit(1);
        }

        initialiseWriter();

        JsonParser jsonParser = new JsonParser();

        try
        {
            for(String sampleId : mSampleIds)
            {
                String filename = String.format("%s/report-%s-hla.json", mSourceDir, sampleId);
                StringBuilder sj = new StringBuilder();
                Files.readAllLines(new File(filename).toPath()).forEach(x -> sj.append(x));
                String jsonInput = sj.toString();

                final JsonElement mainElement = jsonParser.parse(jsonInput);
                final JsonObject mainObject = mainElement.getAsJsonObject();

                JsonArray alleles = mainObject.getAsJsonObject().get("hla").getAsJsonObject().getAsJsonArray("alleles");

                writeSampleData(sampleId, alleles);
            }
        }
        catch(Exception e)
        {
            LL_LOGGER.error("failed to parse json file or data: {}", e.toString());
        }

        closeBufferedWriter(mWriter);
    }

    private void writeSampleData(final String sampleId, final JsonArray alleleData)
    {
        try
        {
            mWriter.write(String.format("%s", sampleId));

            for(int i = 0; i < alleleData.size(); ++i)
            {
                String type = alleleData.get(i).getAsString();
                mWriter.write(String.format(",%s", type));
            }

            mWriter.newLine();
        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to write sample candidates file: {}", e.toString());
        }
    }

    private void initialiseWriter()
    {
        try
        {
            final String outputFileName = mOutputDir + "CPTAC_HLA_TYPES.csv";
            mWriter = createBufferedWriter(outputFileName, false);

            mWriter.write("SampleId,Xhla_A1,Xhla_A2,Xhla_B1,Xhla_B2,Xhla_C1,Xhla_C2");
            mWriter.newLine();

        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to initialise PCAWG types file: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        LL_LOGGER.info("processing CPTAC HLA types files");

        Options options = new Options();
        options.addOption(CPTAC_DIR, true, "Directory with CPTAC json file");
        options.addOption(CPTAC_SAMPLES, true, "File with CPTAC sampleIds");
        options.addOption(OUTPUT_DIR, true, "Path to output");

        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args);

        CptacSummaries cptacSummaries = new CptacSummaries(cmd);
        cptacSummaries.run();

        LL_LOGGER.info("cohort processing complete");
    }

}
