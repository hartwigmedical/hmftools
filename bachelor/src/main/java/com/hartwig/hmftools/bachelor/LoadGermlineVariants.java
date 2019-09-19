package com.hartwig.hmftools.bachelor;

import static com.hartwig.hmftools.bachelor.types.BachelorConfig.SAMPLE_DATA_DIR;
import static com.hartwig.hmftools.bachelor.types.BachelorConfig.databaseAccess;

import java.sql.SQLException;
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
    private static final Logger LOGGER = LogManager.getLogger(LoadGermlineVariants.class);

    public static final String SAMPLE = "sample";
    public static final String DB_USER = "db_user";
    public static final String DB_PASS = "db_pass";
    public static final String DB_URL = "db_url";

    public static void main(@NotNull final String[] args) throws ParseException, SQLException
    {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        final DatabaseAccess dbAccess = databaseAccess(cmd);

        if(dbAccess == null)
        {
            LOGGER.error("failed to create DB connection");
            return;
        }

        final String sampleId = cmd.getOptionValue(SAMPLE);
        final String dataPath = cmd.getOptionValue(SAMPLE_DATA_DIR);

        if(dataPath.isEmpty())
        {
            LOGGER.error("no input data path specified");
            return;
        }

        try
        {
            List<GermlineVariant> germlineVariants = GermlineVariantFile.read(GermlineVariantFile.generateFilename(dataPath, sampleId));

            LOGGER.info("sample({}) loading {} germline records", sampleId, germlineVariants.size());
            final GermlineVariantDAO germlineDAO = new GermlineVariantDAO(dbAccess.context());
            germlineDAO.write(sampleId, germlineVariants);
        }
        catch(Exception e)
        {
            LOGGER.error("error loading and writing germline variants: {}", e.toString());
            return;
        }

        LOGGER.info("Complete");
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
