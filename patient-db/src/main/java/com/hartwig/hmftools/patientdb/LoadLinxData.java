package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.sv.linx.LinxBreakend;
import com.hartwig.hmftools.common.sv.linx.LinxCluster;
import com.hartwig.hmftools.common.sv.linx.LinxDriver;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.common.sv.linx.LinxLink;
import com.hartwig.hmftools.common.sv.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.sv.linx.LinxViralInsertion;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.patientdb.dao.StructuralVariantFusionDAO;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LoadLinxData {

    private static final Logger LOGGER = LogManager.getLogger(LoadLinxData.class);

    private static final String SAMPLE = "sample";
    private static final String LINX_DIR = "linx_dir";

    public static void main(@NotNull String[] args) throws ParseException, IOException
    {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);
        DatabaseAccess dbAccess = createDatabaseAccess(cmd);

        if (dbAccess == null)
        {
            LOGGER.error("Failed to create DB connection");
            System.exit(1);
        }

        String sampleId = cmd.getOptionValue(SAMPLE);
        String linxDir = cmd.getOptionValue(LINX_DIR);

        loadLinxData(dbAccess, sampleId, linxDir);

        LOGGER.info("Complete");
    }

    private static void loadLinxData(@NotNull DatabaseAccess dbAccess, @NotNull String sampleId, @NotNull String linxDir)
            throws IOException
    {
        List<LinxSvAnnotation> svAnnotations = LinxSvAnnotation.read(LinxSvAnnotation.generateFilename(linxDir, sampleId));
        LOGGER.info("Sample({}) loading {} SV annotation records", sampleId, svAnnotations.size());
        dbAccess.writeSvLinxData(sampleId, svAnnotations);

        List<LinxCluster> clusters = LinxCluster.read(LinxCluster.generateFilename(linxDir, sampleId));
        LOGGER.info("Sample({}) loading {} SV cluster records", sampleId, clusters.size());
        dbAccess.writeSvClusters(sampleId, clusters);

        List<LinxLink> links = LinxLink.read(LinxLink.generateFilename(linxDir, sampleId));
        LOGGER.info("Sample({}) loading {} SV links records", sampleId, links.size());
        dbAccess.writeSvLinks(sampleId, links);

        List<LinxBreakend> breakends = LinxBreakend.read(LinxBreakend.generateFilename(linxDir, sampleId));
        List<LinxFusion> fusions = LinxFusion.read(LinxFusion.generateFilename(linxDir, sampleId));
        LOGGER.info("Sample({}) loading {} breakends and {} fusion records", sampleId, breakends.size(), fusions.size());
        StructuralVariantFusionDAO annotationDAO = new StructuralVariantFusionDAO(dbAccess.context());
        annotationDAO.writeBreakendsAndFusions(sampleId, breakends, fusions);

        List<DriverCatalog> driverCatalog = DriverCatalogFile.read(LinxDriver.generateCatalogFilenameForReading(linxDir, sampleId));
        LOGGER.info("Sample({}) loading {} driver catalog records", sampleId, driverCatalog.size());
        dbAccess.writeLinxDriverCatalog(sampleId, driverCatalog);

        List<LinxDriver> drivers = LinxDriver.read(LinxDriver.generateFilename(linxDir, sampleId));
        LOGGER.info("Sample({}) loading {} SV driver records", sampleId, drivers.size());
        dbAccess.writeSvDrivers(sampleId, drivers);
    }

    @NotNull
    private static Options createOptions()
    {
        Options options = new Options();
        addDatabaseCmdLineArgs(options);
        options.addOption(SAMPLE, true, "Name of the tumor sample");
        options.addOption(LINX_DIR, true, "Directory to read LINX data from");

        return options;
    }
}
