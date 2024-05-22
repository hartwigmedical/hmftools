package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.patientdb.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.patientdb.CommonUtils.LOGGER;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.common.peach.PeachGenotypeFile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class LoadPeachData
{
    private static final String PEACH_GENOTYPE_TSV = "peach_genotype_tsv";

    public static void main(@NotNull String[] args) throws ParseException, SQLException, IOException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addConfigItem(SAMPLE, SAMPLE_DESC);
        addDatabaseCmdLineArgs(configBuilder, true);
        configBuilder.addPath(PEACH_GENOTYPE_TSV, true, "Path towards the PEACH genotype TSV file");

        configBuilder.checkAndParseCommandLine(args);

        String sample = configBuilder.getValue(SAMPLE);

        String peachGenotypeTxt = configBuilder.getValue(PEACH_GENOTYPE_TSV);

        try(DatabaseAccess dbWriter = databaseAccess(configBuilder))
        {
            LOGGER.info("Reading PEACH genotypes from {}", peachGenotypeTxt);
            List<PeachGenotype> peachGenotype = PeachGenotypeFile.read(peachGenotypeTxt);
            LOGGER.info(" Read {} PEACH genotypes", peachGenotype.size());

            LOGGER.info("Writing PEACH into database for {}", sample);
            dbWriter.writePeach(sample, peachGenotype);

            LOGGER.info("Complete");
        }
        catch (Exception e)
        {
            LOGGER.error("Failed to load PEACH data", e);
            System.exit(1);
        }
    }
}
