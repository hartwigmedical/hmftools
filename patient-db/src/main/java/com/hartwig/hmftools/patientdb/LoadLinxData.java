package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.variant.structural.linx.LinxBreakend;
import com.hartwig.hmftools.common.variant.structural.linx.LinxCluster;
import com.hartwig.hmftools.common.variant.structural.linx.LinxDriver;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;
import com.hartwig.hmftools.common.variant.structural.linx.LinxLink;
import com.hartwig.hmftools.common.variant.structural.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.variant.structural.linx.LinxViralInsertion;
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
    private static final String SV_DATA_DIR = "sv_data_dir";
    private static final String VCF_FILE = "sv_vcf";

    public static void main(@NotNull String[] args) throws ParseException, IOException {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);
        DatabaseAccess dbAccess = createDatabaseAccess(cmd);

        if (dbAccess == null) {
            LOGGER.error("Failed to create DB connection");
            return;
        }

        String sampleId = cmd.getOptionValue(SAMPLE);
        String svDataPath = cmd.getOptionValue(SV_DATA_DIR);

        loadLinxData(dbAccess, sampleId, svDataPath);

        LOGGER.info("sample({}) data loading complete", sampleId);
    }

    private static void loadLinxData(@NotNull DatabaseAccess dbAccess, @NotNull String sampleId, @NotNull String svDataOutputDir)
            throws IOException {
        List<LinxSvAnnotation> linxSvData = LinxSvAnnotation.read(LinxSvAnnotation.generateFilename(svDataOutputDir, sampleId));
        LOGGER.info("Sample({}) loading {} SV annotation records", sampleId, linxSvData.size());
        dbAccess.writeSvLinxData(sampleId, linxSvData);

        List<LinxCluster> clusterData = LinxCluster.read(LinxCluster.generateFilename(svDataOutputDir, sampleId));
        LOGGER.info("Sample({}) loading {} SV cluster records", sampleId, clusterData.size());
        dbAccess.writeSvClusters(sampleId, clusterData);

        List<LinxLink> linksData = LinxLink.read(LinxLink.generateFilename(svDataOutputDir, sampleId));
        LOGGER.info("Sample({}) loading {} SV links records", sampleId, linksData.size());
        dbAccess.writeSvLinks(sampleId, linksData);

        String viralInsertFilename = LinxViralInsertion.generateFilename(svDataOutputDir, sampleId);
        if (Files.exists(Paths.get(viralInsertFilename))) {
            List<LinxViralInsertion> viralInserts = LinxViralInsertion.read(viralInsertFilename);

            LOGGER.info("Sample({}) loading {} SV viral inserts records", sampleId, viralInserts.size());

            dbAccess.writeSvViralInserts(sampleId, viralInserts);
        }

        String fusionsFilename = LinxFusion.generateFilename(svDataOutputDir, sampleId);
        String breakendsFilename = LinxBreakend.generateFilename(svDataOutputDir, sampleId);

        if (Files.exists(Paths.get(breakendsFilename))) {
            List<LinxBreakend> breakends = LinxBreakend.read(breakendsFilename);

            List<LinxFusion> fusions = Files.exists(Paths.get(fusionsFilename)) ? LinxFusion.read(fusionsFilename) : Lists.newArrayList();

            LOGGER.info("Sample({}) loading {} breakends and {} fusion records", sampleId, breakends.size(), fusions.size());

            StructuralVariantFusionDAO annotationDAO = new StructuralVariantFusionDAO(dbAccess.context());
            annotationDAO.writeBreakendsAndFusions(sampleId, breakends, fusions);
        }

        String driverCatalogFilename = LinxDriver.generateCatalogFilenameForReading(svDataOutputDir, sampleId);

        if (Files.exists(Paths.get(driverCatalogFilename))) {
            List<DriverCatalog> drivers = DriverCatalogFile.read(driverCatalogFilename);
            LOGGER.info("Sample({}) loading {} driver catalog records", sampleId, drivers.size());
            dbAccess.writeLinxDriverCatalog(sampleId, drivers);
        }

        String driversFilename = LinxDriver.generateFilename(svDataOutputDir, sampleId);

        if (Files.exists(Paths.get(driversFilename))) {
            List<LinxDriver> drivers = LinxDriver.read(driversFilename);
            LOGGER.info("Sample({}) loading {} SV driver records", sampleId, drivers.size());
            dbAccess.writeSvDrivers(sampleId, drivers);
        }
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();
        addDatabaseCmdLineArgs(options);
        options.addOption(SAMPLE, true, "Name of the tumor sample. This should correspond to the value used in PURPLE");
        options.addOption(VCF_FILE, true, "Path to the PURPLE structural variant VCF file");
        options.addOption(SV_DATA_DIR, true, "Directory to read or write SV data");

        return options;
    }
}
