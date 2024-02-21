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

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.virus.VirusBreakend;
import com.hartwig.hmftools.common.virus.VirusBreakendFile;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class LoadVirusBreakendData
{
    private static final String VIRUS_BREAKEND_TSV = "virus_breakend_tsv";

    public static void main(@NotNull String[] args) throws ParseException, SQLException, IOException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addConfigItem(SAMPLE, SAMPLE_DESC);
        addDatabaseCmdLineArgs(configBuilder, true);
        configBuilder.addPath(VIRUS_BREAKEND_TSV, true, "Path towards the virus breakend TSV file");

        configBuilder.checkAndParseCommandLine(args);

        String sample = configBuilder.getValue(SAMPLE);
        String virusBreakendTsv = configBuilder.getValue(VIRUS_BREAKEND_TSV);

        try (DatabaseAccess dbWriter = databaseAccess(configBuilder))
        {
            LOGGER.info("Reading virus breakend TSV {}", virusBreakendTsv);
            List<VirusBreakend> virusBreakends = VirusBreakendFile.read(virusBreakendTsv);
            LOGGER.info(" Read {} virus breakends", virusBreakends.size());

            LOGGER.info("Writing virus breakends into database for {}", sample);
            dbWriter.writeVirusBreakend(sample, virusBreakends);

            LOGGER.info("Complete");
        }
        catch (Exception e)
        {
            LOGGER.error("Failed to load virus breakends", e);
            System.exit(1);
        }
    }
}
