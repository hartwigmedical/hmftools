package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.common.cli.DriverGenePanelConfig;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelFactory;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LoadDriverGenePanel {

    private static final Logger LOGGER = LogManager.getLogger(LoadDriverGenePanel.class);

    private static final String ASSEMBLY = "assembly";

    public static void main(@NotNull final String[] args) throws ParseException, IOException, SQLException {
        Options options = new Options();
        // TODO: Rename ref genome versions
        options.addOption(ASSEMBLY, true, "Must be one of [hg19, hg38]");

        DatabaseAccess.addDatabaseCmdLineArgs(options);
        DriverGenePanelConfig.addGenePanelOption(true, options);
        CommandLine cmd = new DefaultParser().parse(options, args);

        String refGenomeVersionArg = cmd.getOptionValue(ASSEMBLY, "hg19");

        RefGenomeVersion refGenomeVersion = refGenomeVersionArg.equals("hg19") ? RefGenomeVersion.V37 : RefGenomeVersion.V38;
        List<DriverGene> driverGenes = DriverGenePanelConfig.driverGenes(cmd);
        DriverGenePanel panel = DriverGenePanelFactory.create(refGenomeVersion, driverGenes);

        DatabaseAccess dbAccess = databaseAccess(cmd);
        dbAccess.writeGenePanel(panel);
        LOGGER.info("Complete");
    }
}
