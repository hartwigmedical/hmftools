package com.hartwig.hmftools.bachelor;

import static com.hartwig.hmftools.bachelor.types.BachelorConfig.BACH_LOGGER;
import static com.hartwig.hmftools.bachelor.types.BachelorConfig.DB_PASS;
import static com.hartwig.hmftools.bachelor.types.BachelorConfig.DB_URL;
import static com.hartwig.hmftools.bachelor.types.BachelorConfig.DB_USER;
import static com.hartwig.hmftools.bachelor.types.BachelorConfig.SAMPLE;
import static com.hartwig.hmftools.bachelor.types.BachelorConfig.databaseAccess;

import java.util.List;

import com.hartwig.hmftools.bachelor.types.GermlineVariant;
import com.hartwig.hmftools.bachelor.types.GermlineVariantFile;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LoadGermlineVariants
{
    private static final String SAMPLE_DATA_DIR = "sample_data_dir";

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        final DatabaseAccess dbAccess = databaseAccess(cmd);

        if(dbAccess == null)
        {
            BACH_LOGGER.error("failed to create DB connection");
            return;
        }

        final String sampleId = cmd.getOptionValue(SAMPLE);
        final String dataPath = cmd.getOptionValue(SAMPLE_DATA_DIR);

        if(dataPath.isEmpty())
        {
            BACH_LOGGER.error("no input data path specified");
            return;
        }

        try
        {
            final List<GermlineVariant> germlineVariants = GermlineVariantFile.read(GermlineVariantFile.generateFilename(dataPath, sampleId));

            if(!germlineVariants.isEmpty())
            {
                BACH_LOGGER.info("sample({}) loading {} germline records", sampleId, germlineVariants.size());
                final GermlineVariantDAO germlineDAO = new GermlineVariantDAO(dbAccess.context());
                germlineDAO.write(sampleId, germlineVariants);
            }
        }
        catch(Exception e)
        {
            BACH_LOGGER.error("error loading and writing germline variants: {}", e.toString());
            return;
        }

        BACH_LOGGER.info("Complete");
    }


    @NotNull
    private static Options createBasicOptions()
    {
        final Options options = new Options();
        options.addOption(SAMPLE, true, "Name of the tumor sample");
        options.addOption(SAMPLE_DATA_DIR, true, "Path for germline variants TSV file");
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password");
        options.addOption(DB_URL, true, "Database url");

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
