package com.hartwig.hmftools.patientdb;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.driver.DriverType.DRIVERS_LINX_SOMATIC;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.patientdb.CommonUtils.LOGGER;
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
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverCatalogFile;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxCluster;
import com.hartwig.hmftools.common.linx.LinxDriver;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.linx.LinxLink;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantFileLoader;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.common.utils.file.FileWriterUtils;
import com.hartwig.hmftools.common.variant.SmallVariant;
import com.hartwig.hmftools.common.variant.SmallVariantFactory;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

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

    public FileDbLoadChecker(final ConfigBuilder configBuilder) throws Exception
    {
        String outputDir = parseOutputDir(configBuilder);

        mWriter = initialiseOutputFile(outputDir);

        mDbAccess = databaseAccess(configBuilder);

        mSampleIds = Lists.newArrayList();

        if(configBuilder.hasValue(SAMPLE))
        {
            mSampleIds.add(configBuilder.getValue(SAMPLE));
        }
        else
        {
            mSampleIds.addAll(ConfigUtils.loadSampleIdsFile(configBuilder));
        }

        mPurpleDir = checkAddDirSeparator(configBuilder.getValue(PURPLE_DIR_CFG));
        mLinxDir = checkAddDirSeparator(configBuilder.getValue(LINX_DIR_CFG));
        mThreads = parseThreads(configBuilder);
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

        final List<Callable<Void>> callableList = sampleTasks.stream().collect(Collectors.toList());
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

    public class SampleCheckerTask implements Callable<Void>
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
        public Void call()
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
            return null;
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

            SmallVariantFactory somaticVariantFactory = new SmallVariantFactory();

            try
            {
                List<SmallVariant> somaticVariants = somaticVariantFactory.fromVCFFile(sampleId, somaticVcf);
                checkFileVsDatabase(sampleId, "SomaticVariant", somaticVariants.size(), "somaticVariant");

                final String geneCopyNumberFile = GeneCopyNumberFile.generateFilename(purpleDir, sampleId);
                final String copyNumberFile = PurpleCopyNumberFile.generateFilenameForReading(purpleDir, sampleId);
                final String svVcf = PurpleCommon.purpleSomaticSvFile(purpleDir, sampleId);

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

                int fileDriverCount = (int) driverCatalog.stream().filter(x -> DRIVERS_LINX_SOMATIC.contains(x.driver())).count();
                int dbDriverCount = (int) dbDrivers.stream().filter(x -> DRIVERS_LINX_SOMATIC.contains(x.driver())).count();

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

    public static void main(@NotNull final String[] args) throws Exception
    {
        ConfigBuilder configBuilder = new ConfigBuilder();

        configBuilder.addConfigItem(SAMPLE, false, "Tumor sample ID");
        addSampleIdFile(configBuilder, false);
        configBuilder.addPath(PURPLE_DIR_CFG, true, PURPLE_DIR_DESC);
        configBuilder.addPath(LINX_DIR_CFG, true, LINX_DIR_DESC);
        addThreadOptions(configBuilder);
        addDatabaseCmdLineArgs(configBuilder, true);
        ConfigUtils.addLoggingOptions(configBuilder);
        FileWriterUtils.addOutputOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        setLogLevel(configBuilder);
        FileDbLoadChecker fileDbLoadChecker = new FileDbLoadChecker(configBuilder);
        fileDbLoadChecker.run();
    }
}
