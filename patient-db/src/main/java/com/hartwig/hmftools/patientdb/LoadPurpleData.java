package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.common.purple.PurpleCommon.purpleSomaticSvFile;
import static com.hartwig.hmftools.common.purple.PurpleCommon.purpleSomaticVcfFile;
import static com.hartwig.hmftools.common.sv.StructuralVariantData.convertSvData;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.patientdb.CommonUtils.LOGGER;
import static com.hartwig.hmftools.patientdb.CommonUtils.logVersion;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.FittedPurityRangeFile;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariantFactory;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.sv.StructuralVariantFileLoader;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;
import com.hartwig.hmftools.patientdb.dao.BufferedWriter;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class LoadPurpleData
{
    private static final String RNA_SAMPLE = "rna";
    private static final String DB_SAMPLE = "db_sample";

    private static final String SOMATIC_ONLY = "somatic_only";
    private static final String GERMLINE_ONLY = "germline_only";

    public static void main(@NotNull String[] args)
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

        try(DatabaseAccess dbAccess = databaseAccess(configBuilder))
        {
            String sampleId = configBuilder.getValue(SAMPLE);
            String dbSampleId = configBuilder.hasValue(DB_SAMPLE) ? configBuilder.getValue(DB_SAMPLE) : sampleId;
            String referenceId = configBuilder.getValue(REFERENCE);
            String rnaId = configBuilder.getValue(RNA_SAMPLE);

            if(sampleId == null && referenceId == null)
            {
                LOGGER.error("missing sample or reference ID config");
                System.exit(1);
            }

            if(sampleId == null)
                sampleId = referenceId;

            String purpleDir = checkAddDirSeparator(configBuilder.getValue(PURPLE_DIR_CFG));

            if(!Files.exists(Paths.get(purpleDir)))
            {
                LOGGER.error("invalid Purple data directory({})", purpleDir);
                System.exit(1);
            }

            boolean loadGermline = !configBuilder.hasFlag(SOMATIC_ONLY);
            boolean loadSomatic = !configBuilder.hasFlag(GERMLINE_ONLY);

            LOGGER.info("loading sample({}) {} Purple data from {}",
                    sampleId,
                    loadGermline & loadSomatic ? "somatic and germline" : (loadSomatic ? "somatic" : "germline"),
                    purpleDir);

            final String sample = sampleId;

            dbAccess.context().transaction(tr ->
            {
                loadCommonData(dbSampleId, sample, dbAccess, purpleDir);

                if(loadSomatic)
                    loadSomaticData(dbSampleId, sample, referenceId, rnaId, dbAccess, purpleDir);

                if(loadGermline)
                    loadGermlineData(dbSampleId, sample, referenceId, rnaId, dbAccess, purpleDir);
            });

            LOGGER.info("Purple data loading complete");
        }
        catch(Exception e)
        {
            LOGGER.error("Failed to load Purple data", e);
            System.exit(1);
        }
    }

    private static void loadCommonData(
            final String dbSample, final String sampleId, final DatabaseAccess dbAccess, final String purpleDir) throws Exception
    {
        LOGGER.info("loading common data");
        PurityContext purityContext = PurityContextFile.read(purpleDir, sampleId);

        dbAccess.writePurity(dbSample, purityContext, purityContext.qc());
    }

    private static void loadSomaticData(
            final String dbSampleId, final String sampleId, final String referenceId, final String rnaId,
            final DatabaseAccess dbAccess, final String purpleDir) throws Exception
    {
        // check all somatic files exist before attempting to load
        final String geneCopyNumberFile = GeneCopyNumberFile.generateFilename(purpleDir, sampleId);
        final String copyNumberFile = PurpleCopyNumberFile.generateFilenameForReading(purpleDir, sampleId);
        final String somaticDriversFile = DriverCatalogFile.generateSomaticFilename(purpleDir, sampleId);
        final String somaticVcf = purpleSomaticVcfFile(purpleDir, sampleId);
        final String svVcf = purpleSomaticSvFile(purpleDir, sampleId);

        // skip loading if any files are missing
        List<String> requiredFiles = Lists.newArrayList(
                geneCopyNumberFile, copyNumberFile, somaticDriversFile, somaticVcf, svVcf);

        if(requiredFiles.stream().noneMatch(x -> Files.exists(Paths.get(x))))
        {
            LOGGER.info("skipping somatic data - no files present");
            return;
        }

        if(hasMissingFiles(requiredFiles, "somatic"))
            System.exit(1);

        List<GeneCopyNumber> geneCopyNumbers = GeneCopyNumberFile.read(geneCopyNumberFile);
        List<PurpleCopyNumber> copyNumbers = PurpleCopyNumberFile.read(copyNumberFile);
        List<FittedPurity> bestFitPerPurity = FittedPurityRangeFile.readBestFitPerPurity(purpleDir, sampleId);
        List<DriverCatalog> somaticDriverCatalog = DriverCatalogFile.read(somaticDriversFile);

        dbAccess.writeBestFitPerPurity(dbSampleId, bestFitPerPurity);
        dbAccess.writeCopynumbers(dbSampleId, copyNumbers);
        dbAccess.writeGeneCopyNumbers(dbSampleId, geneCopyNumbers);
        dbAccess.writePurpleDriverCatalog(dbSampleId, somaticDriverCatalog, null);

        LOGGER.info("loading geneCopyNumber({}) copyNumber({}) purityFits({}) somaticDrivers({})",
                geneCopyNumbers.size(), copyNumbers.size(), bestFitPerPurity.size(), somaticDriverCatalog.size());

        List<StructuralVariant> variants = StructuralVariantFileLoader.fromFile(svVcf, new AlwaysPassFilter());
        List<EnrichedStructuralVariant> enrichedVariants = new EnrichedStructuralVariantFactory().enrich(variants);

        // Generate a unique ID for each SV record
        int svId = 0;

        List<StructuralVariantData> structuralVariants = Lists.newArrayList();
        for (EnrichedStructuralVariant variant : enrichedVariants)
        {
            structuralVariants.add(convertSvData(variant, svId++));
        }

        LOGGER.info("loading {} SVs", structuralVariants.size());
        dbAccess.writeStructuralVariants(dbSampleId, structuralVariants);

        BufferedWriter<SomaticVariant> somaticWriter = dbAccess.somaticVariantWriter(dbSampleId);

        SomaticVariantFactory somaticVariantFactory = new SomaticVariantFactory();

        somaticVariantFactory.fromVCFFile(sampleId, referenceId, rnaId, somaticVcf, referenceId != null, somaticWriter);
        somaticWriter.close();

        LOGGER.info("loaded {} somatic variants, filtered({})",
                somaticVariantFactory.getCreatedCount(), somaticVariantFactory.getFilteredCount());
    }

    private static void loadGermlineData(
            final String dbSampleId, final String sampleId, final String referenceId, final String rnaId,
            final DatabaseAccess dbAccess, final String purpleDir) throws Exception
    {
        final String germlineVcf = PurpleCommon.purpleGermlineVcfFile(purpleDir, sampleId);
        final String germlineDeletionsFile = GermlineDeletion.generateFilename(purpleDir, sampleId);
        final String germlineDriverFile = DriverCatalogFile.generateGermlineFilename(purpleDir, sampleId);

        // skip loading if any files are missing
        List<String> requiredFiles = Lists.newArrayList(
                germlineVcf, germlineDeletionsFile, germlineDriverFile);

        if(requiredFiles.stream().noneMatch(x -> Files.exists(Paths.get(x))))
        {
            LOGGER.info("skipping germline data - no files present");
            return;
        }

        if(hasMissingFiles(requiredFiles, "germline"))
            System.exit(1);

        List<GermlineDeletion> germlineDeletions = GermlineDeletion.read(germlineDeletionsFile);

        List<DriverCatalog> germlineDriverCatalog = DriverCatalogFile.read(germlineDriverFile);

        LOGGER.info("loading germline drivers({}) deletions({})", germlineDriverCatalog.size(), germlineDeletions.size());

        dbAccess.writeGermlineDeletions(dbSampleId, germlineDeletions);
        dbAccess.writePurpleDriverCatalog(dbSampleId, null, germlineDriverCatalog);

        int variantCount = 0;

        try(AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(germlineVcf, new VCFCodec(), false);

        BufferedWriter<VariantContext> dbWriter = dbAccess.germlineVariantWriter(dbSampleId, referenceId, rnaId))
        {
            dbWriter.initialise();

            for(VariantContext context : reader.iterator())
            {
                dbWriter.accept(context);
                ++variantCount;
            }
        }

        LOGGER.info("loaded {} germline variants", variantCount);
    }

    public static boolean hasMissingFiles(final List<String> requiredFiles, final String sourceType)
    {
        if(requiredFiles.stream().allMatch(x -> Files.exists(Paths.get(x))))
            return false;

        LOGGER.warn("skipping {} data - files are missing:", sourceType);

        for(String filename : requiredFiles)
        {
            if(!Files.exists(Paths.get(filename)))
            {
                LOGGER.warn("missing file {}", filename);
            }
        }

        return true;
    }

    private static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, SAMPLE_DESC);
        configBuilder.addConfigItem(DB_SAMPLE, "ID of the sample in the database (optional). Defaults to the sample ID.");
        configBuilder.addConfigItem(REFERENCE, REFERENCE_DESC);
        configBuilder.addConfigItem(RNA_SAMPLE, "RNA sample ID");
        configBuilder.addConfigItem(PURPLE_DIR_CFG, true, PURPLE_DIR_DESC);
        configBuilder.addFlag(SOMATIC_ONLY, "Only load somatic data");
        configBuilder.addFlag(GERMLINE_ONLY, "Only load germline data");
        addDatabaseCmdLineArgs(configBuilder, true);
        ConfigUtils.addLoggingOptions(configBuilder);
    }
}
