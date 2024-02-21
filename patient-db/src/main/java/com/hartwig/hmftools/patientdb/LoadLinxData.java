package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.common.drivercatalog.DriverType.DRIVERS_LINX_GERMLINE;
import static com.hartwig.hmftools.common.drivercatalog.DriverType.DRIVERS_LINX_SOMATIC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_GERMLINE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_GERMLINE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.patientdb.LoadPurpleData.hasMissingFiles;
import static com.hartwig.hmftools.patientdb.CommonUtils.LOGGER;
import static com.hartwig.hmftools.patientdb.CommonUtils.logVersion;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxCluster;
import com.hartwig.hmftools.common.linx.LinxDriver;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.linx.LinxGermlineSv;
import com.hartwig.hmftools.common.linx.LinxLink;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.patientdb.dao.StructuralVariantFusionDAO;

import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class LoadLinxData
{
    private static final String SOMATIC_ONLY = "somatic_only";
    private static final String GERMLINE_ONLY = "germline_only";

    public static void main(@NotNull String[] args) throws ParseException, IOException
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        addConfig(configBuilder);

        if(!configBuilder.parseCommandLine(args))
        {
            configBuilder.logInvalidDetails();
            System.exit(1);
        }

        setLogLevel(configBuilder);
        logVersion();

        try (DatabaseAccess dbAccess = createDatabaseAccess(configBuilder))
        {
            if(dbAccess == null)
            {
                LOGGER.error("Failed to create DB connection");
                System.exit(1);
            }

            String sampleId = configBuilder.getValue(SAMPLE);
            String linxDir = configBuilder.getValue(LINX_DIR_CFG);
            String linxGermlineDir = configBuilder.hasValue(LINX_GERMLINE_DIR_CFG) ? configBuilder.getValue(LINX_GERMLINE_DIR_CFG) : linxDir;

            boolean loadGermline = !configBuilder.hasFlag(SOMATIC_ONLY);
            boolean loadSomatic = !configBuilder.hasFlag(GERMLINE_ONLY);

            dbAccess.context().transaction(tr ->
            {
                if(loadSomatic)
                    loadSomaticData(dbAccess, sampleId, linxDir);

                if(loadGermline)
                    loadGermlineData(dbAccess, sampleId, linxGermlineDir);
            });

            LOGGER.info("Linx data loading complete");
        }
        catch (Exception e)
        {
            LOGGER.error("Failed to load Linx data", e);
            System.exit(1);
        }
    }

    private static void loadSomaticData(final DatabaseAccess dbAccess, final String sampleId, final String linxDir)
            throws IOException
    {
        LOGGER.info("sample({}) loading Linx somatic data", sampleId);

        final String svAnnotationFile = LinxSvAnnotation.generateFilename(linxDir, sampleId, false);
        final String svClusterFile = LinxCluster.generateFilename(linxDir, sampleId, false);
        final String svLinkFile = LinxLink.generateFilename(linxDir, sampleId, false);
        final String svBreakendFile = LinxBreakend.generateFilename(linxDir, sampleId, false);
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
        final String germlineBreakendFile = LinxBreakend.generateFilename(linxDir, sampleId, true);
        final String driverCatalogFile = LinxDriver.generateCatalogFilename(linxDir, sampleId, false);

        List<String> requiredFiles = Lists.newArrayList(germlineSvFile, driverCatalogFile); // required after v5.31

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

        if(Files.exists(Paths.get(germlineBreakendFile)))
        {
            List<LinxBreakend> germlineBreakends = LinxBreakend.read(germlineBreakendFile);
            LOGGER.info("sample({}) loading {} germline breakend records", sampleId, germlineBreakends.size());
            dbAccess.writeGermlineBreakends(sampleId, germlineBreakends);
        }
    }

    private static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        configBuilder.addConfigItem(LINX_DIR_CFG, true, LINX_DIR_DESC);
        configBuilder.addConfigItem(LINX_GERMLINE_DIR_CFG, false, LINX_GERMLINE_DIR_DESC);
        configBuilder.addFlag(SOMATIC_ONLY, "Only load somatic data");
        configBuilder.addFlag(GERMLINE_ONLY, "Only load germline data");
        addDatabaseCmdLineArgs(configBuilder, true);
        ConfigUtils.addLoggingOptions(configBuilder);
    }
}
