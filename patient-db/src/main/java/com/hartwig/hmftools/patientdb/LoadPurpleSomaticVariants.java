package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.IOException;
import java.sql.SQLException;

import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.patientdb.dao.BufferedWriter;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LoadPurpleSomaticVariants
{
    private static final Logger LOGGER = LogManager.getLogger(LoadPurpleSomaticVariants.class);

    private static final String SAMPLE = "sample";
    private static final String REFERENCE = "reference";
    private static final String RNA = "rna";

    private static final String SOMATIC_VCF = "somatic_vcf";

    public static void main(@NotNull String[] args) throws ParseException, SQLException
    {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        DatabaseAccess dbAccess = databaseAccess(cmd);

        String tumorSample = cmd.getOptionValue(SAMPLE);
        String referenceSample = cmd.getOptionValue(REFERENCE, null);
        String rnaSample = cmd.getOptionValue(RNA, null);

        String somaticVcf = cmd.getOptionValue(SOMATIC_VCF);

        LOGGER.info("Removing old data of sample {}", tumorSample);

        try
        {
            BufferedWriter<SomaticVariant> somaticWriter = dbAccess.somaticVariantWriter(tumorSample);

            LOGGER.info("Streaming data from {} to db", somaticVcf);
            new SomaticVariantFactory().fromVCFFile(tumorSample, referenceSample, rnaSample, somaticVcf, true, somaticWriter);
        }
        catch(IOException e)
        {
            LOGGER.error("failed to load somatic VCF({}): {}", somaticVcf, e.toString());
            System.exit(1);
        }

        LOGGER.info("Complete");
    }

    @NotNull
    private static Options createOptions()
    {
        Options options = new Options();
        options.addOption(SAMPLE, true, "Name of the tumor sample. This should correspond to the value used in PURPLE.");
        options.addOption(REFERENCE, true, "Optional name of the reference sample. This should correspond to the value used in PURPLE.");
        options.addOption(RNA, true, "Optional name of the rna sample. This should correspond to the value used in PURPLE.");
        options.addOption(SOMATIC_VCF, true, "Path to the PURPLE somatic variant VCF file.");
        addDatabaseCmdLineArgs(options);

        return options;
    }
}
