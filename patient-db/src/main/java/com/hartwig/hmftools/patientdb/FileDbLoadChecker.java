package com.hartwig.hmftools.patientdb;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.drivercatalog.DriverType.DRIVERS_LINX_SOMATIC;
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
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberFile;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantFileLoader;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxCluster;
import com.hartwig.hmftools.common.linx.LinxDriver;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.linx.LinxLink;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.utils.ConfigUtils;
import com.hartwig.hmftools.common.utils.FileWriterUtils;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;
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
    private final int mThreads;

    private static final String SAMPLE = "sample";
    private static final String SAMPLE_ID_FILE = "sample_id_file";
    private static final String PURPLE_DIR = "purple_dir";
    private static final String LINX_DIR = "linx_dir";
    private static final String THREADS = "threads";

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
        mThreads = Integer.parseInt(cmd.getOptionValue(THREADS, "1"));
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

        List<SampleCheckerTask> sampleTasks = Lists.newArrayList();

        for(int i = 0; i < min(mSampleIds.size(), mThreads); ++i)
        {
            sampleTasks.add(new SampleCheckerTask(i));
        }

        int taskIndex = 0;
        for(String sampleId : mSampleIds)
        {
            if(taskIndex >= sampleTasks.size())
                taskIndex = 0;

            sampleTasks.get(taskIndex).getSampleIds().add(sampleId);

            ++taskIndex;
        }

        final List<Callable> callableList = sampleTasks.stream().collect(Collectors.toList());
        TaskExecutor.executeTasks(callableList, mThreads);

        closeBufferedWriter(mWriter);

        LOGGER.info("File DB loading check complete");
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

    private synchronized void writeDiff(final String sampleId, final String fileType, int fileCount, int dbCount)
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

    public class SampleCheckerTask implements Callable
    {
        private final int mTaskId;
        private final List<String> mSampleIds;

        public SampleCheckerTask(int taskId)
        {
            mTaskId = taskId;
            mSampleIds = Lists.newArrayList();
        }

        public List<String> getSampleIds() { return mSampleIds; }

        @Override
        public Long call()
        {
            for(int i = 0; i < mSampleIds.size(); ++i)
            {
                String sampleId = mSampleIds.get(i);

                processSample(sampleId);

                if(i > 0 && (i % 100) == 0)
                {
                    LOGGER.info("{}: processed {} samples", mTaskId, i);
                }
            }

            LOGGER.info("{}: tasks complete for {} samples", mTaskId, mSampleIds.size());
            return (long)0;
        }

        private void processSample(final String sampleId)
        {
            LOGGER.debug("processing sample({})", sampleId);

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
        }

        private void checkPurpleFiles(final String sampleId, final String purpleDir)
        {
            if(!Files.exists(Paths.get(purpleDir)))
            {
                LOGGER.error("invalid Purple data directory({})", purpleDir);
                return;
            }

            final String somaticVcf = PurpleCommon.purpleSomaticVcfFile(purpleDir, sampleId);

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

                final String geneCopyNumberFile = GeneCopyNumberFile.generateFilenameForReading(purpleDir, sampleId);
                final String copyNumberFile = PurpleCopyNumberFile.generateFilenameForReading(purpleDir, sampleId);
                final String svVcf = PurpleCommon.purpleSvFile(purpleDir, sampleId);

                List<GeneCopyNumber> geneCopyNumbers = GeneCopyNumberFile.read(geneCopyNumberFile);
                List<PurpleCopyNumber> copyNumbers = PurpleCopyNumberFile.read(copyNumberFile);
                List<StructuralVariant> svs = StructuralVariantFileLoader.fromFile(svVcf, new AlwaysPassFilter());

                checkFileVsDatabase(sampleId, "GeneCopyNumber", geneCopyNumbers.size(), "geneCopyNumber");
                checkFileVsDatabase(sampleId, "CopyNumber", copyNumbers.size(), "copyNumber");
                checkFileVsDatabase(sampleId, "StructuralVariant", svs.size(), "structuralVariant");
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
                List<LinxCluster> clusters = LinxCluster.read(svClusterFile);
                List<LinxLink> links = LinxLink.read(svLinkFile);
                List<LinxBreakend> breakends = LinxBreakend.read(svBreakendFile);
                List<LinxFusion> fusions = LinxFusion.read(svFusionFile);
                List<LinxDriver> drivers = LinxDriver.read(svDriverFile);

                checkFileVsDatabase(sampleId, "SvAnnotation", svAnnotations.size(), "svAnnotation");
                checkFileVsDatabase(sampleId, "SvBreakend", breakends.size(), "svBreakend");
                checkFileVsDatabase(sampleId, "SvFusion", fusions.size(), "svFusion");
                checkFileVsDatabase(sampleId, "SvDriver", drivers.size(), "svDriver");
                checkFileVsDatabase(sampleId, "SvCluster", clusters.size(), "svCluster");
                checkFileVsDatabase(sampleId, "SvLink", links.size(), "svLink");

                List<DriverCatalog> driverCatalog = DriverCatalogFile.read(driverCatalogFile);
                List<DriverCatalog> dbDrivers = mDbAccess.readDriverCatalog(sampleId);

                int fileDriverCount = (int)driverCatalog.stream().filter(x -> DRIVERS_LINX_SOMATIC.contains(x.driver())).count();
                int dbDriverCount = (int)dbDrivers.stream().filter(x -> DRIVERS_LINX_SOMATIC.contains(x.driver())).count();

                if(dbDriverCount != fileDriverCount)
                {
                    LOGGER.trace("sample({}) driverCatalog file({}) DB({}) mismatch", sampleId, fileDriverCount, dbDriverCount);
                    writeDiff(sampleId, "DriverCatalog", fileDriverCount, dbDriverCount);
                }
                else
                {
                    LOGGER.trace("sample({}) driverCatalog({}) count({}) match", sampleId, fileDriverCount);
                }
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

        private int getDbVariantCount(final String sampleId, final DatabaseAccess dbAccess, final String tableName)
        {
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
        options.addOption(THREADS, true, "Thread count");
        addDatabaseCmdLineArgs(options);
        ConfigUtils.addLoggingOptions(options);
        FileWriterUtils.addOutputOptions(options);
        return options;
    }
}
