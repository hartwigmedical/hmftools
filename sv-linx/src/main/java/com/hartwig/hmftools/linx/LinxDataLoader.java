package com.hartwig.hmftools.linx;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantData.convertSvData;
import static com.hartwig.hmftools.linx.LinxConfig.SAMPLE;
import static com.hartwig.hmftools.linx.LinxConfig.SV_DATA_DIR;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;
import static com.hartwig.hmftools.patientdb.dao.DatabaseUtil.valueNotNull;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariantFactory;
import com.hartwig.hmftools.common.variant.structural.ImmutableStructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFile;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFileLoader;
import com.hartwig.hmftools.common.variant.structural.linx.LinxBreakend;
import com.hartwig.hmftools.common.variant.structural.linx.LinxCluster;
import com.hartwig.hmftools.common.variant.structural.linx.LinxDriver;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;
import com.hartwig.hmftools.common.variant.structural.linx.LinxLink;
import com.hartwig.hmftools.common.variant.structural.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.variant.structural.linx.LinxViralInsertion;
import com.hartwig.hmftools.linx.types.GermlineFilter;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.patientdb.dao.DatabaseUtil;
import com.hartwig.hmftools.patientdb.dao.StructuralVariantFusionDAO;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LinxDataLoader
{
    private static final Logger LOGGER = LogManager.getLogger(LinxDataLoader.class);

    public static final String VCF_FILE = "sv_vcf";

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final DatabaseAccess dbAccess = createDatabaseAccess(cmd);

        if(dbAccess == null)
        {
            LOGGER.error("failed to create DB connection");
            return;
        }

        final String sampleId = cmd.getOptionValue(SAMPLE);
        final String svDataPath = cmd.getOptionValue(SV_DATA_DIR);

        loadLinxData(dbAccess, sampleId, svDataPath);

        LOGGER.info("sample({}) data loading complete", sampleId);
    }

    private static void loadLinxData(final DatabaseAccess dbAccess, final String sampleId, final String svDataOutputDir)
    {
        try
        {
            List<LinxSvAnnotation> linxSvData = LinxSvAnnotation.read(LinxSvAnnotation.generateFilename(svDataOutputDir, sampleId));
            LOGGER.info("sample({}) loading {} SV annotation records", sampleId, linxSvData.size());
            dbAccess.writeSvLinxData(sampleId, linxSvData);

            List<LinxCluster> clusterData = LinxCluster.read(LinxCluster.generateFilename(svDataOutputDir, sampleId));
            LOGGER.info("sample({}) loading {} SV cluster records", sampleId, clusterData.size());
            dbAccess.writeSvClusters(sampleId, clusterData);

            List<LinxLink> linksData = LinxLink.read(LinxLink.generateFilename(svDataOutputDir, sampleId));
            LOGGER.info("sample({}) loading {} SV links records", sampleId, linksData.size());
            dbAccess.writeSvLinks(sampleId, linksData);

            String viralInsertFilename = LinxViralInsertion.generateFilename(svDataOutputDir, sampleId);
            if(Files.exists(Paths.get(viralInsertFilename)))
            {
                List<LinxViralInsertion> viralInserts = LinxViralInsertion.read(viralInsertFilename);

                if (!viralInserts.isEmpty())
                {
                    LOGGER.info("sample({}) loading {} SV viral inserts records", sampleId, viralInserts.size());
                }

                dbAccess.writeSvViralInserts(sampleId, viralInserts);
            }

            final String fusionsFilename = LinxFusion.generateFilename(svDataOutputDir, sampleId);
            final String breakendsFilename = LinxBreakend.generateFilename(svDataOutputDir, sampleId);

            if(Files.exists(Paths.get(breakendsFilename)))
            {
                List<LinxBreakend> breakends = LinxBreakend.read(breakendsFilename);

                List<LinxFusion> fusions = Files.exists(Paths.get(fusionsFilename)) ?
                        LinxFusion.read(fusionsFilename) : Lists.newArrayList();

                LOGGER.info("sample({}) loading {} breakends and {} fusion records", sampleId, breakends.size(), fusions.size());

                final StructuralVariantFusionDAO annotationDAO = new StructuralVariantFusionDAO(dbAccess.context());
                annotationDAO.writeBreakendsAndFusions(sampleId, breakends, fusions);
            }

            final String driverCatalogFilename = DriverCatalogFile.generateSomaticFilenameForReading(svDataOutputDir, sampleId);

            if(Files.exists(Paths.get(driverCatalogFilename)))
            {
                List<DriverCatalog> drivers = DriverCatalogFile.read(driverCatalogFilename);
                LOGGER.info("sample({}) loading {} driver catalog records", sampleId, drivers.size());
                dbAccess.writeLinxDriverCatalog(sampleId, drivers);
            }

            final String driversFilename = LinxDriver.generateFilename(svDataOutputDir, sampleId);

            if(Files.exists(Paths.get(driversFilename)))
            {
                List<LinxDriver> drivers = LinxDriver.read(driversFilename);
                LOGGER.info("sample({}) loading {} SV driver records", sampleId, drivers.size());
                dbAccess.writeSvDrivers(sampleId, drivers);
            }

        }
        catch(IOException e)
        {
            LOGGER.error("failed to load SV data files: {}", e.toString());
        }
    }


    public static List<StructuralVariantData> loadSvDataFromVcf(final String vcfFile)
    {
        final List<StructuralVariantData> svDataList = Lists.newArrayList();

        try
        {
            final List<StructuralVariant> variants = StructuralVariantFileLoader.fromFile(vcfFile, new AlwaysPassFilter());
            final List<EnrichedStructuralVariant> enrichedVariants = new EnrichedStructuralVariantFactory().enrich(variants);

            // generate a unique ID for each SV record
            int svId = 0;

            for (EnrichedStructuralVariant var : enrichedVariants)
            {
                svDataList.add(convertSvData(var, svId++));
            }
        }
        catch(IOException e)
        {
            LOGGER.error("failed to load SVs from VCF: {}", e.toString());
        }

        LOGGER.info("loaded {} SV data records from VCF file: {}", svDataList.size(), vcfFile);

        return svDataList;
    }

    public static List<StructuralVariantData> loadSvDataFromGermlineVcf(final String vcfFile)
    {
        final List<StructuralVariantData> svDataList = Lists.newArrayList();

        try
        {
            final List<StructuralVariant> variants = StructuralVariantFileLoader.fromFile(vcfFile, new GermlineFilter(true));

            int svId = 0;

            for (StructuralVariant var : variants)
            {
                svDataList.add(convertGermlineSvData(var, svId++));
            }
        }
        catch(IOException e)
        {
            LOGGER.error("failed to load SVs from VCF: {}", e.toString());
        }

        LOGGER.info("loaded {} germline SV data records from VCF file: {}", svDataList.size(), vcfFile);

        return svDataList;
    }

    public static List<StructuralVariantData> loadSvDataFromSvFile(final String sampleId, final String svDataPath)
    {
        try
        {
            final String svDataFile = StructuralVariantFile.generateFilename(svDataPath, sampleId);
            return StructuralVariantFile.read(svDataFile);
        }
        catch(IOException e)
        {
            LOGGER.error("failed to load SV data: {}", e.toString());
            return Lists.newArrayList();
        }
    }

    public static StructuralVariantData convertGermlineSvData(final StructuralVariant var, int svId)
    {
        return ImmutableStructuralVariantData.builder()
                .id(svId)
                .startChromosome(var.chromosome(true))
                .endChromosome(var.end() == null ? "0" : var.chromosome(false))
                .startPosition(var.position(true).intValue())
                .endPosition(var.end() == null ? -1 : var.position(false).intValue())
                .startOrientation(var.orientation(true))
                .endOrientation(var.end() == null ? (byte) 0 : var.orientation(false))
                .startHomologySequence(var.start().homology())
                .endHomologySequence(var.end() == null ? "" : var.end().homology())
                .junctionCopyNumber(1)
                .startAF(DatabaseUtil.valueNotNull(var.start().alleleFrequency()))
                .endAF(var.end() == null ? 0 : DatabaseUtil.valueNotNull(var.end().alleleFrequency()))
                .adjustedStartAF(DatabaseUtil.valueNotNull(var.start().alleleFrequency()))
                .adjustedEndAF(var.end() == null ? 0 : DatabaseUtil.valueNotNull(var.end().alleleFrequency()))
                .adjustedStartCopyNumber(DatabaseUtil.valueNotNull(1))
                .adjustedEndCopyNumber(var.end() == null ? 0 : 1)
                .adjustedStartCopyNumberChange(1)
                .adjustedEndCopyNumberChange(var.end() == null ? 0 : 1)
                .insertSequence(var.insertSequence())
                .type(var.type())
                .filter(var.filter())
                .imprecise(var.imprecise())
                .qualityScore(DatabaseUtil.valueNotNull(var.qualityScore()))
                .event(valueNotNull(var.event()))
                .startTumorVariantFragmentCount(DatabaseUtil.valueNotNull(var.start().tumorVariantFragmentCount()))
                .startTumorReferenceFragmentCount(DatabaseUtil.valueNotNull(var.start().tumorReferenceFragmentCount()))
                .startNormalVariantFragmentCount(DatabaseUtil.valueNotNull(var.start().normalVariantFragmentCount()))
                .startNormalReferenceFragmentCount(DatabaseUtil.valueNotNull(var.start().normalReferenceFragmentCount()))
                .endTumorVariantFragmentCount(0)
                .endTumorReferenceFragmentCount(0)
                .endNormalVariantFragmentCount(var.end() == null ? 0 : DatabaseUtil.valueNotNull(var.end().normalVariantFragmentCount()))
                .endNormalReferenceFragmentCount(var.end() == null ? 0 : DatabaseUtil.valueNotNull(var.end().normalReferenceFragmentCount()))
                .startIntervalOffsetStart(DatabaseUtil.valueNotNull(var.start().startOffset()))
                .startIntervalOffsetEnd(DatabaseUtil.valueNotNull(var.start().endOffset()))
                .endIntervalOffsetStart(var.end() == null ? 0 : DatabaseUtil.valueNotNull(var.end().startOffset()))
                .endIntervalOffsetEnd(var.end() == null ? 0 : DatabaseUtil.valueNotNull(var.end().endOffset()))
                .inexactHomologyOffsetStart(DatabaseUtil.valueNotNull(var.start().inexactHomologyOffsetStart()))
                .inexactHomologyOffsetEnd(DatabaseUtil.valueNotNull(var.start().inexactHomologyOffsetEnd()))
                .startLinkedBy(valueNotNull(var.startLinkedBy()))
                .endLinkedBy(valueNotNull(var.endLinkedBy()))
                .vcfId(valueNotNull(var.id()))
                .startRefContext("") // getValueNotNull(var.start().refGenomeContext())
                .endRefContext(var.end() == null ? "" : "") // getValueNotNull(var.end().refGenomeContext())
                .recovered(var.recovered())
                .recoveryMethod((valueNotNull(var.recoveryMethod())))
                .recoveryFilter(valueNotNull(var.recoveryFilter()))
                .insertSequenceAlignments(valueNotNull(var.insertSequenceAlignments()))
                .insertSequenceRepeatClass(valueNotNull(var.insertSequenceRepeatClass()))
                .insertSequenceRepeatType(valueNotNull(var.insertSequenceRepeatType()))
                .insertSequenceRepeatOrientation(DatabaseUtil.valueNotNull(var.insertSequenceRepeatOrientation()))
                .insertSequenceRepeatCoverage(DatabaseUtil.valueNotNull(var.insertSequenceRepeatCoverage()))
                .startAnchoringSupportDistance(var.start().anchoringSupportDistance())
                .endAnchoringSupportDistance(var.end() == null ? 0 : var.end().anchoringSupportDistance())
                .build();
    }

    @NotNull
    private static Options createBasicOptions()
    {
        final Options options = new Options();
        addDatabaseCmdLineArgs(options);
        options.addOption(SAMPLE, true, "Name of the tumor sample. This should correspond to the value used in PURPLE");
        options.addOption(VCF_FILE, true, "Path to the PURPLE structural variant VCF file");
        options.addOption(SV_DATA_DIR, true, "Directory to read or write SV data");

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
