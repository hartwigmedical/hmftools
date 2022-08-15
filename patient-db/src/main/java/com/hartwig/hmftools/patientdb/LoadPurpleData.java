package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.common.purple.PurpleCommon.purpleSomaticVcfFile;
import static com.hartwig.hmftools.common.purple.PurpleCommon.purpleSvFile;
import static com.hartwig.hmftools.common.sv.StructuralVariantData.convertSvData;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.FittedPurityRangeFile;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariantFactory;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.sv.StructuralVariantFileLoader;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;
import com.hartwig.hmftools.patientdb.dao.BufferedWriter;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class LoadPurpleData 
{
    private static final Logger LOGGER = LogManager.getLogger(LoadPurpleData.class);

    private static final String SAMPLE = "sample";
    private static final String REFERENCE = "reference";
    private static final String RNA_SAMPLE = "rna";

    private static final String PURPLE_DIR = "purple_dir";

    private static final String SOMATIC_ONLY = "somatic_only";
    private static final String GERMLINE_ONLY = "germline_only";

    public static void main(@NotNull String[] args)
    {
        Options options = createOptions();

        try
        {
            CommandLine cmd = new DefaultParser().parse(options, args);
            DatabaseAccess dbAccess = databaseAccess(cmd);

            String sampleId = cmd.getOptionValue(SAMPLE);
            String referenceId = cmd.getOptionValue(REFERENCE);
            String rnaId = cmd.getOptionValue(RNA_SAMPLE);

            if(sampleId == null && referenceId == null)
            {
                LOGGER.error("missing sample or reference ID config");
                System.exit(1);
            }

            if(sampleId == null)
                sampleId = referenceId;

            String purpleDir = checkAddDirSeparator(cmd.getOptionValue(PURPLE_DIR));

            if(!Files.exists(Paths.get(purpleDir)))
            {
                LOGGER.error("invalid Purple data directory({})", purpleDir);
                System.exit(1);
            }

            boolean loadGermline = !cmd.hasOption(SOMATIC_ONLY);
            boolean loadSomatic = !cmd.hasOption(GERMLINE_ONLY);

            LOGGER.info("loading sample({}) {} Purple data from {}",
                    sampleId,
                    loadGermline & loadSomatic ? "somatic and germline" : (loadSomatic ? "somatic" : "germline"),
                    purpleDir);

            final String sample = sampleId;

            dbAccess.context().transaction(tr ->
            {
                loadCommonData(sample, dbAccess, purpleDir);

                if(loadSomatic)
                    loadSomaticData(sample, referenceId, rnaId, dbAccess, purpleDir);

                if(loadGermline)
                    loadGermlineData(sample, referenceId, rnaId, dbAccess, purpleDir);
            });

            LOGGER.info("Purple data loading complete");
        }
        catch(Exception e)
        {
            LOGGER.error("failed to load Purple data: {}", e.toString());
            e.printStackTrace();
            System.exit(1);
        }
    }

    private static void loadCommonData(
            final String sampleId, final DatabaseAccess dbAccess, final String purpleDir) throws Exception
    {
        LOGGER.info("loading common data");
        PurityContext purityContext = PurityContextFile.read(purpleDir, sampleId);

        dbAccess.writePurity(sampleId, purityContext, purityContext.qc());
    }

    private static void loadSomaticData(
            final String sampleId, final String referenceId, final String rnaId,
            final DatabaseAccess dbAccess, final String purpleDir) throws Exception
    {
        // check all somatic files exist before attempting to load
        final String geneCopyNumberFile = GeneCopyNumberFile.generateFilenameForReading(purpleDir, sampleId);
        final String copyNumberFile = PurpleCopyNumberFile.generateFilenameForReading(purpleDir, sampleId);
        final String purityRangeFile = FittedPurityRangeFile.generateFilenameForReading(purpleDir, sampleId);
        final String somaticDriversFile = DriverCatalogFile.generateSomaticFilename(purpleDir, sampleId);
        final String somaticVcf = purpleSomaticVcfFile(purpleDir, sampleId);
        final String svVcf = purpleSvFile(purpleDir, sampleId);

        // skip loading if any files are missing
        List<String> requiredFiles = Lists.newArrayList(
                geneCopyNumberFile, copyNumberFile, somaticDriversFile, purityRangeFile, somaticVcf, svVcf);

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

        dbAccess.writeBestFitPerPurity(sampleId, bestFitPerPurity);
        dbAccess.writeCopynumbers(sampleId, copyNumbers);
        dbAccess.writeGeneCopyNumbers(sampleId, geneCopyNumbers);
        dbAccess.writePurpleDriverCatalog(sampleId, somaticDriverCatalog, null);

        LOGGER.info("loading geneCopyNumber({}) copyNumber({}) purityFits({}) somaticDrivers({})",
                geneCopyNumbers.size(), copyNumbers.size(), bestFitPerPurity.size(), somaticDriverCatalog.size());

        List<StructuralVariant> variants = StructuralVariantFileLoader.fromFile(svVcf, new AlwaysPassFilter());
        List<EnrichedStructuralVariant> enrichedVariants = new EnrichedStructuralVariantFactory().enrich(variants);

        // Generate a unique ID for each SV record
        int svId = 0;

        List<StructuralVariantData> structuralVariants = com.google.common.collect.Lists.newArrayList();
        for (EnrichedStructuralVariant variant : enrichedVariants)
        {
            structuralVariants.add(convertSvData(variant, svId++));
        }

        LOGGER.info("loading {} SVs", structuralVariants.size());
        dbAccess.writeStructuralVariants(sampleId, structuralVariants);

        BufferedWriter<SomaticVariant> somaticWriter = dbAccess.somaticVariantWriter(sampleId);

        SomaticVariantFactory somaticVariantFactory = new SomaticVariantFactory();

        somaticVariantFactory.fromVCFFile(sampleId, referenceId, rnaId, somaticVcf, referenceId != null, somaticWriter);
        somaticWriter.close();

        LOGGER.info("loaded {} somatic variants, filtered({})",
                somaticVariantFactory.getCreatedCount(), somaticVariantFactory.getFilteredCount());
    }

    private static void loadGermlineData(
            final String sampleId, final String referenceId, final String rnaId,
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

        LOGGER.info("loading germlime drivers({}) deletions({})", germlineDriverCatalog.size(), germlineDeletions.size());

        dbAccess.writeGermlineDeletions(sampleId, germlineDeletions);
        dbAccess.writePurpleDriverCatalog(sampleId, null, germlineDriverCatalog);

        int variantCount = 0;

        try(AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(germlineVcf, new VCFCodec(), false);

        BufferedWriter<VariantContext> dbWriter = dbAccess.germlineVariantWriter(sampleId, referenceId, rnaId))
        {
            dbWriter.initialise();

            for(VariantContext context : reader.iterator())
            {
                dbWriter.accept(context);
                ++variantCount;
            }
        }

        LOGGER.info("loaded {} germlime variants", variantCount);
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

    @NotNull
    private static Options createOptions()
    {
        Options options = new Options();
        options.addOption(SAMPLE, true, "Tumor sample ID");
        options.addOption(REFERENCE, true, "Reference ID");
        options.addOption(RNA_SAMPLE, true, "RNA sample ID");
        options.addOption(PURPLE_DIR, true, "Path to the Purple directory");
        options.addOption(SOMATIC_ONLY, false, "Only load somatic data");
        options.addOption(GERMLINE_ONLY, false, "Only load germline data");
        addDatabaseCmdLineArgs(options);
        return options;
    }
}
