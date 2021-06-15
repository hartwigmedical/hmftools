package com.hartwig.hmftools.lilac.cohort;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.cohort.CohortCommon.SAMPLE_FILES_DIR;
import static com.hartwig.hmftools.lilac.cohort.CohortCommon.SAMPLE_IDS_FILE;
import static com.hartwig.hmftools.lilac.cohort.CohortCommon.loadSampleIds;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class PcawgSummaries
{
    private final String mOutputDir;
    private final String mTypesFile;

    private BufferedWriter mWriter;

    private static final String PCAWG_TYPES_FILE = "pcawg_types_file";

    public PcawgSummaries(final CommandLine cmd)
    {
        mOutputDir = parseOutputDir(cmd);
        mTypesFile = cmd.getOptionValue(PCAWG_TYPES_FILE);
    }

    public void run()
    {
        if(mTypesFile == null || mOutputDir == null)
        {
            LL_LOGGER.error("invalid input/output paths");
            System.exit(1);
        }

        initialiseWriter();

        try
        {
            //
            final List<String> lines = Files.readAllLines(new File(mTypesFile).toPath());

            JsonParser jsonParser = new JsonParser();

            for(String line : lines)
            {
                try
                {
                    final JsonElement mainElement = jsonParser.parse(line);
                    final JsonObject mainObject = mainElement.getAsJsonObject();

                    for(Map.Entry<String,JsonElement> entry : mainObject.entrySet())
                    {
                        String sampleId = entry.getKey();
                        JsonElement sampleData = entry.getValue();

                        LL_LOGGER.debug("sampleId({})", sampleId);

                        JsonObject polySolver = sampleData.getAsJsonObject().get("polysolver").getAsJsonObject();
                        JsonObject xHLA = sampleData.getAsJsonObject().get("xHLA").getAsJsonObject();
                        JsonArray lilac = sampleData.getAsJsonObject().get("LILAC").getAsJsonArray();

                        writeSampleData(sampleId, polySolver, xHLA, lilac);
                    }
                }
                catch(Exception je)
                {
                    LL_LOGGER.error("failed to parse json data: {}", je.toString());
                    break;
                }
            }
        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to load PCAWG types file({}): {}", mTypesFile, e.toString());
        }

        // initialiseWriter();
        // closeBufferedWriter(mWriter);
    }

    private void writeSampleData(final String sampleId, final JsonObject polysolver, final JsonObject xHLA, final JsonArray lilac)
    {
        try
        {
            mWriter.write(String.format("%s", sampleId));

            writeHlaTypes(polysolver.getAsJsonArray("normal"));
            writeHlaTypes(polysolver.getAsJsonArray("tumor"));
            writeHlaTypes(xHLA.getAsJsonArray("normal"));
            writeHlaTypes(xHLA.getAsJsonArray("tumor"));
            writeHlaTypes(lilac);
            mWriter.newLine();

        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to write sample candidates file: {}", e.toString());
        }
    }

    private void writeHlaTypes(final JsonArray hlaTypes) throws IOException
    {
        for(int i = 0; i < hlaTypes.size(); ++i)
        {
            String type = hlaTypes.get(i).getAsString();
            mWriter.write(String.format(",%s", type));
        }
    }

    private void initialiseWriter()
    {
        try
        {
            final String outputFileName = mOutputDir + "PCAWG_HLA_TYPES.csv";
            mWriter = createBufferedWriter(outputFileName, false);

            mWriter.write("SampleId,PsR_A1,PsR_A2,PsR_B1,PsR_B2,PsR_C1,PsR_C2");
            mWriter.write(",PsT_A1,PsT_A2,PsT_B1,PsT_B2,PsT_C1,PsT_C2");
            mWriter.write(",XhlaR_A1,XhlaR_A2,XhlaR_B1,XhlaR_B2,XhlaR_C1,XhlaR_C2");
            mWriter.write(",XhlaT_A1,XhlaT_A2,XhlaT_B1,XhlaT_B2,XhlaT_C1,XhlaT_C2");
            mWriter.write(",LiR_A1,LiR_A2,LiR_B1,LiR_B2,LiR_C1,LiR_C2");

            mWriter.newLine();

        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to initialise PCAWG types file: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        LL_LOGGER.info("processing PCAWG HLA types files");

        Options options = new Options();
        options.addOption(PCAWG_TYPES_FILE, true, "PCAWG json file");
        options.addOption(OUTPUT_DIR, true, "Path to output");

        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args);

        PcawgSummaries pcawgSummaries = new PcawgSummaries(cmd);
        pcawgSummaries.run();

        LL_LOGGER.info("cohort processing complete");
    }

}
