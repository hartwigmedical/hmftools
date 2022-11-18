package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LoadDriverGenePanel
{

    private static final Logger LOGGER = LogManager.getLogger(LoadDriverGenePanel.class);

    public static void main(@NotNull final String[] args) throws ParseException, IOException, SQLException
    {
        Options options = new Options();

        DatabaseAccess.addDatabaseCmdLineArgs(options);
        DriverGenePanelConfig.addGenePanelOption(true, options);
        CommandLine cmd = new DefaultParser().parse(options, args);

        List<DriverGene> driverGenes = DriverGenePanelConfig.driverGenes(cmd);

        LOGGER.info("Loading {} driver genes to database", driverGenes.size());

        DatabaseAccess dbAccess = databaseAccess(cmd);
        dbAccess.writeGenePanel(driverGenes);
        LOGGER.info("Complete");
    }
}
