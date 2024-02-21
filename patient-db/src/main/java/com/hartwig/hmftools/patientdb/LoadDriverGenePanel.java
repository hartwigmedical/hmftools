package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION;
import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION_DESC;
import static com.hartwig.hmftools.patientdb.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.patientdb.CommonUtils.LOGGER;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class LoadDriverGenePanel
{
    public static void main(@NotNull final String[] args) throws ParseException, IOException, SQLException
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        DatabaseAccess.addDatabaseCmdLineArgs(configBuilder, true);
        configBuilder.addConfigItem(DRIVER_GENE_PANEL_OPTION, DRIVER_GENE_PANEL_OPTION_DESC);

        configBuilder.checkAndParseCommandLine(args);

        List<DriverGene> driverGenes = DriverGenePanelConfig.driverGenes(configBuilder);

        LOGGER.info("Loading {} driver genes to database", driverGenes.size());

        try (DatabaseAccess dbAccess = databaseAccess(configBuilder, true))
        {
            dbAccess.writeGenePanel(driverGenes);
            LOGGER.info("Complete");
        }
        catch (Exception e)
        {
            LOGGER.error("Failed to load driver genes", e);
            System.exit(1);
        }
    }
}
