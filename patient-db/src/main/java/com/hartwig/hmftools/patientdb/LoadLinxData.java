package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.common.drivercatalog.DriverType.DRIVERS_LINX_GERMLINE;
import static com.hartwig.hmftools.common.drivercatalog.DriverType.DRIVERS_LINX_SOMATIC;
import static com.hartwig.hmftools.patientdb.LoadPurpleData.hasMissingFiles;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.sv.linx.LinxBreakend;
import com.hartwig.hmftools.common.sv.linx.LinxCluster;
import com.hartwig.hmftools.common.sv.linx.LinxDriver;
import com.hartwig.hmftools.common.sv.linx.LinxFusion;
import com.hartwig.hmftools.common.sv.linx.LinxGermlineSv;
import com.hartwig.hmftools.common.sv.linx.LinxLink;
import com.hartwig.hmftools.common.sv.linx.LinxSvAnnotation;
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

    private static final String SOMATIC_ONLY = "somatic_only";
    private static final String GERMLINE_ONLY = "germline_only";

    public static void main(@NotNull String[] args) throws ParseException, IOException
    {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);
        DatabaseAccess dbAccess = createDatabaseAccess(cmd);

        if(dbAccess == null)
        {
            LOGGER.error("Failed to create DB connection");
            System.exit(1);
        }

        String sampleId = cmd.getOptionValue(SAMPLE);
        String linxDir = cmd.getOptionValue(LINX_DIR);

        boolean loadGermline = !cmd.hasOption(SOMATIC_ONLY);
        boolean loadSomatic = !cmd.hasOption(GERMLINE_ONLY);

        dbAccess.context().transaction(tr ->
        {
            if(loadSomatic)
                loadSomaticData(dbAccess, sampleId, linxDir);

            if(loadGermline)
                loadGermlineData(dbAccess, sampleId, linxDir);
        });

        LOGGER.info("Linx data loading complete");
    }

    private static void loadSomaticData(final DatabaseAccess dbAccess, final String sampleId, final String linxDir)
            throws IOException
    {
        LOGGER.info("sample({}) loading Linx somatic data", sampleId);

        final String svAnnotationFile = LinxSvAnnotation.generateFilename(linxDir, sampleId, false);
        final String svClusterFile = LinxCluster.generateFilename(linxDir, sampleId, false);
        final String svLinkFile = LinxLink.generateFilename(linxDir, sampleId, false);
        final String svBreakendFile = LinxBreakend.generateFilename(linxDir, sampleId);
        final String svFusionFile = LinxFusion.generateFilename(linxDir, sampleId);
        final String svDriverFile = LinxDriver.generateFilename(linxDir, sampleId);
        final String driverCatalogFile = LinxDriver.generateCatalogFilename(linxDir, sampleId, true);

        List<String> requiredFiles = Lists.newArrayList(
                svAnnotationFile, svClusterFile, svLinkFile, svBreakendFile, svFusionFile, svDriverFile, driverCatalogFile);

        if(requiredFiles.stream().noneMatch(x -> Files.exists(Paths.get(x))))
        {
            LOGGER.info("skipping somatic data - no files present");
            return;
        }

        if(hasMissingFiles(requiredFiles, "somatic"))
            System.exit(1);

        List<LinxSvAnnotation> svAnnotations = LinxSvAnnotation.read(svAnnotationFile);
        LOGGER.info("sample({}) loading {} SV annotation records", sampleId, svAnnotations.size());
        dbAccess.writeSvLinxData(sampleId, svAnnotations);

        List<LinxCluster> clusters = LinxCluster.read(svClusterFile);
        LOGGER.info("sample({}) loading {} SV cluster records", sampleId, clusters.size());
        dbAccess.writeSvClusters(sampleId, clusters);

        List<LinxLink> links = LinxLink.read(svLinkFile);
        LOGGER.info("sample({}) loading {} SV links records", sampleId, links.size());
        dbAccess.writeSvLinks(sampleId, links);

        List<LinxBreakend> breakends = LinxBreakend.read(svBreakendFile);
        List<LinxFusion> fusions = LinxFusion.read(svFusionFile);
        LOGGER.info("sample({}) loading {} breakends and {} fusion records", sampleId, breakends.size(), fusions.size());
        StructuralVariantFusionDAO annotationDAO = new StructuralVariantFusionDAO(dbAccess.context());
        annotationDAO.writeBreakendsAndFusions(sampleId, breakends, fusions);

        List<LinxDriver> drivers = LinxDriver.read(svDriverFile);
        LOGGER.info("sample({}) loading {} SV driver records", sampleId, drivers.size());
        dbAccess.writeSvDrivers(sampleId, drivers);

        List<DriverCatalog> driverCatalog = DriverCatalogFile.read(driverCatalogFile);
        LOGGER.info("sample({}) loading {} somatic driver catalog records", sampleId, driverCatalog.size());
        dbAccess.writeLinxDriverCatalog(sampleId, driverCatalog, DRIVERS_LINX_SOMATIC);
    }

    private static void loadGermlineData(final DatabaseAccess dbAccess, final String sampleId, final String linxDir)
            throws IOException
    {
        LOGGER.info("sample({}) loading Linx germline data", sampleId);

        final String germlineSvFile = LinxGermlineSv.generateFilename(linxDir, sampleId);
        final String driverCatalogFile = LinxDriver.generateCatalogFilename(linxDir, sampleId, false);

        List<String> requiredFiles = Lists.newArrayList(germlineSvFile, driverCatalogFile);

        if(requiredFiles.stream().noneMatch(x -> Files.exists(Paths.get(x))))
        {
            LOGGER.info("skipping germline data - no files present");
            return;
        }

        if(hasMissingFiles(requiredFiles, "germline"))
            System.exit(1);

        List<DriverCatalog> driverCatalog = DriverCatalogFile.read(driverCatalogFile);
        LOGGER.info("sample({}) loading {} germline driver catalog records", sampleId, driverCatalog.size());
        dbAccess.writeLinxDriverCatalog(sampleId, driverCatalog, DRIVERS_LINX_GERMLINE);

        List<LinxGermlineSv> germlineSVs = LinxGermlineSv.read(germlineSvFile);
        LOGGER.info("sample({}) loading {} germline SV records", sampleId, germlineSVs.size());
        dbAccess.writeGermlineSVs(sampleId, germlineSVs);
    }

    @NotNull
    private static Options createOptions()
    {
        Options options = new Options();
        addDatabaseCmdLineArgs(options);
        options.addOption(SAMPLE, true, "Name of the tumor sample");
        options.addOption(LINX_DIR, true, "Directory to read LINX data from");
        options.addOption(SOMATIC_ONLY, false, "Only load somatic data");
        options.addOption(GERMLINE_ONLY, false, "Only load germline data");

        return options;
    }
}
