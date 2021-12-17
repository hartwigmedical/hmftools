package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.common.drivercatalog.DriverType.DRIVERS_LINX_GERMLINE;
import static com.hartwig.hmftools.common.drivercatalog.DriverType.DRIVERS_LINX_SOMATIC;
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

    private static final String DATA_TYPE = "data_type";
    private static final String DATA_TYPE_BOTH = "both";
    private static final String DATA_TYPE_SOMATIC = "somatic";
    private static final String DATA_TYPE_GERMLINE = "germline";

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

        String dataType = cmd.getOptionValue(DATA_TYPE, DATA_TYPE_BOTH);

        loadLinxData(dbAccess, sampleId, linxDir, dataType);

        LOGGER.info("Complete");
    }

    private static void loadLinxData(final DatabaseAccess dbAccess, final String sampleId, final String linxDir, final String dataTypes)
            throws IOException
    {
        boolean loadSomaticData = dataTypes.equals(DATA_TYPE_BOTH) || dataTypes.equals(DATA_TYPE_SOMATIC);
        boolean loadGermlineData = dataTypes.equals(DATA_TYPE_BOTH) || dataTypes.equals(DATA_TYPE_GERMLINE);

        if (!(loadSomaticData || loadGermlineData))
            LOGGER.warn("No data will be loaded based on selected datatype '{}'", dataTypes);

        if(loadSomaticData)
        {
            LOGGER.info("Sample({}) loading Linx somatic data", sampleId);

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

            List<LinxDriver> drivers = LinxDriver.read(LinxDriver.generateFilename(linxDir, sampleId));
            LOGGER.info("Sample({}) loading {} SV driver records", sampleId, drivers.size());
            dbAccess.writeSvDrivers(sampleId, drivers);

            List<DriverCatalog> driverCatalog = DriverCatalogFile.read(LinxDriver.generateCatalogFilename(linxDir, sampleId, true));
            LOGGER.info("Sample({}) loading {} somatic driver catalog records", sampleId, driverCatalog.size());
            dbAccess.writeLinxDriverCatalog(sampleId, driverCatalog, DRIVERS_LINX_SOMATIC);
        }

        if(loadGermlineData)
        {
            LOGGER.info("Sample({}) loading Linx germline data", sampleId);

            List<DriverCatalog> driverCatalog = DriverCatalogFile.read(LinxDriver.generateCatalogFilename(linxDir, sampleId, false));
            LOGGER.info("Sample({}) loading {} germline driver catalog records", sampleId, driverCatalog.size());
            dbAccess.writeLinxDriverCatalog(sampleId, driverCatalog, DRIVERS_LINX_GERMLINE);

            List<LinxGermlineSv> germlineSVs = LinxGermlineSv.read(LinxGermlineSv.generateFilename(linxDir, sampleId));
            LOGGER.info("Sample({}) loading {} germline SV records", sampleId, germlineSVs.size());
            dbAccess.writeGermlineSVs(sampleId, germlineSVs);
        }
    }

    @NotNull
    private static Options createOptions()
    {
        Options options = new Options();
        addDatabaseCmdLineArgs(options);
        options.addOption(SAMPLE, true, "Name of the tumor sample");
        options.addOption(LINX_DIR, true, "Directory to read LINX data from");
        options.addOption(DATA_TYPE, true, "Data type(s) to load: 'somatic', 'germline' or 'both' (default)");

        return options;
    }
}
