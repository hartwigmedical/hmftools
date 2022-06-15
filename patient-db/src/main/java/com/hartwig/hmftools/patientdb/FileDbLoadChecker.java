package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.common.purple.PurpleCommon.PURPLE_SOMATIC_VCF_SUFFIX;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.patientdb.LoadPurpleData.hasMissingFiles;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.BufferedWriter;
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
import com.hartwig.hmftools.common.sv.linx.LinxLink;
import com.hartwig.hmftools.common.sv.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.utils.ConfigUtils;
import com.hartwig.hmftools.common.utils.FileWriterUtils;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jooq.Record;
import org.jooq.Record1;
import org.jooq.Result;

public class FileDbLoadChecker
{
    private final BufferedWriter mWriter;
    private final String mPurpleDir;
    private final String mLinxDir;
    private final List<String> mSampleIds;
    private final DatabaseAccess mDbAccess;

    private static final String SAMPLE = "sample";
    private static final String SAMPLE_ID_FILE = "sample_id_file";
    private static final String PURPLE_DIR = "purple_dir";
    private static final String LINX_DIR = "linx_dir";

    private static final Logger LOGGER = LogManager.getLogger(LoadPurpleSomaticVariants.class);

    public FileDbLoadChecker(final CommandLine cmd) throws Exception
    {
        String outputDir = parseOutputDir(cmd);

        mWriter = initialiseOutputFile(outputDir);

        mDbAccess = databaseAccess(cmd);

        mSampleIds = Lists.newArrayList();

        if(cmd.hasOption(SAMPLE))
        {
            mSampleIds.add(cmd.getOptionValue(SAMPLE));
        }
        else
        {
            mSampleIds.addAll(ConfigUtils.loadSampleIdsFile(cmd.getOptionValue(SAMPLE_ID_FILE)));
        }

        mPurpleDir = cmd.hasOption(PURPLE_DIR) ? checkAddDirSeparator(cmd.getOptionValue(PURPLE_DIR)) : null;
        mLinxDir = cmd.hasOption(LINX_DIR) ? checkAddDirSeparator(cmd.getOptionValue(LINX_DIR)) : null;
    }

    public void run()
    {
        if(mWriter == null || mDbAccess == null)
            System.exit(1);

        if(mSampleIds.isEmpty())
        {
            LOGGER.error("missing sample ID(s) config");
            System.exit(1);
        }

        LOGGER.info("processing {} samples", mSampleIds.size());

        int processed = 0;
        for(String sampleId : mSampleIds)
        {
            LOGGER.debug("processing sample({})", sampleId);

            ++processed;

            if(mPurpleDir != null)
            {
                String purpleDir = mPurpleDir.replaceAll("\\*", sampleId);
                checkPurpleFiles(sampleId, purpleDir);
            }

            if(mLinxDir != null)
            {
                String linxDir = mLinxDir.replaceAll("\\*", sampleId);
                checkLinxFiles(sampleId, linxDir);
            }

            if(processed > 0 && (processed % 100) == 0)
            {
                LOGGER.info("processed {} samples", processed);
            }
        }

        closeBufferedWriter(mWriter);

        LOGGER.info("File DB loading check complete");
    }

    private void checkPurpleFiles(final String sampleId, final String purpleDir)
    {
        if(!Files.exists(Paths.get(purpleDir)))
        {
            LOGGER.error("invalid Purple data directory({})", purpleDir);
            return;
        }

        final String somaticVcf = purpleDir + sampleId + PURPLE_SOMATIC_VCF_SUFFIX;

        if(!Files.exists(Paths.get(somaticVcf)))
        {
            LOGGER.warn("somatic VCF({}) missing");
            return;
        }

        SomaticVariantFactory somaticVariantFactory = new SomaticVariantFactory();

        try
        {
            List<SomaticVariant> somaticVariants = somaticVariantFactory.fromVCFFile(sampleId, somaticVcf);

            checkFileVsDatabase(sampleId, "SomaticVariant", somaticVariants.size(), "somaticVariant");
        }
        catch(Exception e)
        {
            LOGGER.error("failed to read sample({}) file: {}", sampleId, e.toString());
        }
    }

    private void checkLinxFiles(final String sampleId, final String linxDir)
    {
        if(!Files.exists(Paths.get(linxDir)))
        {
            LOGGER.error("invalid Linx data directory({})", linxDir);
            return;
        }

        try
        {
            final String svAnnotationFile = LinxSvAnnotation.generateFilename(linxDir, sampleId);
            final String svClusterFile = LinxCluster.generateFilename(linxDir, sampleId);
            final String svLinkFile = LinxLink.generateFilename(linxDir, sampleId);
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
            List<LinxCluster> clusters = LinxCluster.read(svClusterFile);
            List<LinxLink> links = LinxLink.read(svLinkFile);
            List<LinxBreakend> breakends = LinxBreakend.read(svBreakendFile);
            List<LinxFusion> fusions = LinxFusion.read(svFusionFile);
            List<LinxDriver> drivers = LinxDriver.read(svDriverFile);
            List<DriverCatalog> driverCatalog = DriverCatalogFile.read(driverCatalogFile);

            checkFileVsDatabase(sampleId, "SvAnnotation", svAnnotations.size(), "svAnnotation");
            checkFileVsDatabase(sampleId, "SvBreakend", breakends.size(), "svBreakend");
            checkFileVsDatabase(sampleId, "SvFusion", fusions.size(), "svFusion");
            checkFileVsDatabase(sampleId, "DriverCatalog", driverCatalog.size(), "driverCatalog");
            checkFileVsDatabase(sampleId, "SvDriver", drivers.size(), "svDriver");
            checkFileVsDatabase(sampleId, "SvCluster", clusters.size(), "svCluster");
            checkFileVsDatabase(sampleId, "SvLinx", links.size(), "svLink");
        }
        catch(Exception e)
        {
            LOGGER.error("failed to read sample({}) file: {}", sampleId, e.toString());
        }
    }

    private void checkFileVsDatabase(final String sampleId, final String fileType, int fileCount, final String tableName)
    {
        int dbCount = getDbVariantCount(sampleId, mDbAccess, tableName);

        if(dbCount == fileCount)
        {
            LOGGER.trace("sample({}) table({}) count({}) match", sampleId, tableName, fileCount);
            return;
        }

        LOGGER.trace("sample({}) table({}) file({}) DB({}) mismatch", sampleId, tableName, fileCount, dbCount);
        writeDiff(sampleId, fileType, fileCount, dbCount);
    }

    private BufferedWriter initialiseOutputFile(final String outputDir)
    {
        try
        {
            String outputFile = outputDir + "file_db_load_check_results.csv";

            LOGGER.info("writing output results: {}", outputFile);

            BufferedWriter writer = createBufferedWriter(outputFile, false);

            writer.write("SampleId,FileType,FileCount,DbCount");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            LOGGER.error("failed to initialise output file: {}", e.toString());
            return null;
        }
    }

    private void writeDiff(final String sampleId, final String fileType, int fileCount, int dbCount)
    {
        if(mWriter == null)
            return;

        try
        {
            mWriter.write(String.format("%s,%s,%d,%d", sampleId, fileType, fileCount, dbCount));
            mWriter.newLine();
        }
        catch(IOException e)
        {
            LOGGER.error("failed to write output file: {}", e.toString());
            System.exit(1);
        }
    }

    private static int getDbVariantCount(final String sampleId, final DatabaseAccess dbAccess, final String tableName)
    {
        // String queryStr = String.format("select count(*) from %s where sampleId = '%s'", tableName, sampleId);
        // dbAccess.context().query(queryStr).execute();

        Result<Record1<Integer>> result = dbAccess.context()
                .selectCount()
                .from(tableName)
                .where(String.format("sampleId = '%s'", sampleId))
                .fetch();

        for(Record record : result)
        {
            return Integer.parseInt(record.getValue(0).toString());
        }

        return 0;
    }

    public static void main(@NotNull String[] args)
    {
        Options options = createOptions();

        try
        {
            CommandLine cmd = new DefaultParser().parse(options, args);

            setLogLevel(cmd);

            FileDbLoadChecker fileDbLoadChecker = new FileDbLoadChecker(cmd);
            fileDbLoadChecker.run();
        }
        catch(Exception e)
        {
            LOGGER.error("Error: {}", e.toString());
            e.printStackTrace();
            System.exit(1);
        }
    }

    @NotNull
    private static Options createOptions()
    {
        Options options = new Options();
        options.addOption(SAMPLE, true, "Tumor sample ID");
        options.addOption(SAMPLE_ID_FILE, true, "CSV with SampleId");
        options.addOption(PURPLE_DIR, true, "Sample Purple directory");
        options.addOption(LINX_DIR, true, "Sample Linx directory");
        addDatabaseCmdLineArgs(options);
        ConfigUtils.addLoggingOptions(options);
        FileWriterUtils.addOutputOptions(options);
        return options;
    }
}
