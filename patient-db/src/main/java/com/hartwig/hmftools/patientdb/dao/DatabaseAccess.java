package com.hartwig.hmftools.patientdb.dao;

import java.io.IOException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.EnumSet;
import java.util.List;

import com.hartwig.hmftools.common.amber.AmberAnonymous;
import com.hartwig.hmftools.common.cider.Cdr3LocusSummary;
import com.hartwig.hmftools.common.cider.Cdr3Sequence;
import com.hartwig.hmftools.common.cuppa.CuppaPredictions;
import com.hartwig.hmftools.common.teal.TelomereLength;
import com.hartwig.hmftools.patientdb.amber.AmberMapping;
import com.hartwig.hmftools.patientdb.amber.AmberPatient;
import com.hartwig.hmftools.patientdb.amber.AmberSample;
import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.flagstat.Flagstat;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxCluster;
import com.hartwig.hmftools.common.linx.LinxDriver;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.linx.LinxGermlineSv;
import com.hartwig.hmftools.common.linx.LinxLink;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.metrics.WGSMetricWithQC;
import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleQC;
import com.hartwig.hmftools.common.sigs.SignatureAllocation;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.VirusBreakend;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.jooq.DSLContext;
import org.jooq.SQLDialect;
import org.jooq.conf.MappedSchema;
import org.jooq.conf.RenderMapping;
import org.jooq.conf.Settings;
import org.jooq.impl.DSL;

import htsjdk.variant.variantcontext.VariantContext;

public class DatabaseAccess implements AutoCloseable
{
    private static final Logger LOGGER = LogManager.getLogger(DatabaseAccess.class);
    private static final String DEV_CATALOG = "hmfpatients_test";

    public static final double MIN_SAMPLE_PURITY = 0.195;

    public static final String DB_USER = "db_user";
    public static final String DB_PASS = "db_pass";
    public static final String DB_URL = "db_url";

    private static final String DB_USER_DESC = "Database username";
    private static final String DB_PASS_DESC = "Database password";
    private static final String DB_URL_DESC = "Database url";

    public static final String DB_DEFAULT_ARGS = "?serverTimezone=UTC&useSSL=false";

    @NotNull
    private final Connection connection;
    @NotNull
    private final DSLContext context;
    @NotNull
    private final DriverGenePanelDAO driverGenePanelDAO;
    @NotNull
    private final CanonicalTranscriptDAO canonicalTranscriptDAO;
    @NotNull
    private final MetricDAO metricDAO;
    @NotNull
    private final FlagstatDAO flagstatDAO;
    @NotNull
    private final SnpCheckDAO snpCheckDAO;
    @NotNull
    private final SomaticVariantDAO somaticVariantDAO;
    @NotNull
    private final GermlineVariantDAO germlineVariantDAO;
    @NotNull
    private final AmberDAO amberDAO;
    @NotNull
    private final PurityDAO purityDAO;
    @NotNull
    private final CopyNumberDAO copyNumberDAO;
    @NotNull
    private final GeneCopyNumberDAO geneCopyNumberDAO;
    @NotNull
    private final StructuralVariantDAO structuralVariantDAO;
    @NotNull
    private final StructuralVariantClusterDAO structuralVariantClusterDAO;
    @NotNull
    private final StructuralVariantFusionDAO structuralVariantFusionDAO;
    @NotNull
    private final DriverCatalogDAO driverCatalogDAO;
    @NotNull
    private final PeachDAO peachDAO;
    @NotNull
    private final SignatureDAO signatureDAO;
    @NotNull
    private final CuppaDAO cuppaDAO;
    @NotNull
    private final ChordDAO chordDAO;
    @NotNull
    private final VirusBreakendDAO virusBreakendDAO;
    @NotNull
    private final VirusInterpreterDAO virusInterpreterDAO;
    @NotNull
    private final CiderDAO ciderDAO;
    @NotNull
    private final TealDAO tealDAO;

    public DatabaseAccess(@NotNull final String userName, @NotNull final String password, @NotNull final String url) throws SQLException
    {
        System.setProperty("org.jooq.no-logo", "true");
        System.setProperty("org.jooq.no-tips", "true");

        this.connection = DriverManager.getConnection(url, userName, password);
        String catalog = connection.getCatalog();
        LOGGER.debug("Connecting to database '{}'", catalog);
        this.context = DSL.using(connection, SQLDialect.MYSQL, settings(catalog));

        this.driverGenePanelDAO = new DriverGenePanelDAO(context);
        this.canonicalTranscriptDAO = new CanonicalTranscriptDAO(context);
        this.metricDAO = new MetricDAO(context);
        this.flagstatDAO = new FlagstatDAO(context);
        this.snpCheckDAO = new SnpCheckDAO(context);
        this.somaticVariantDAO = new SomaticVariantDAO(context);
        this.germlineVariantDAO = new GermlineVariantDAO(context);
        this.amberDAO = new AmberDAO(context);
        this.purityDAO = new PurityDAO(context);
        this.copyNumberDAO = new CopyNumberDAO(context);
        this.geneCopyNumberDAO = new GeneCopyNumberDAO(context);
        this.structuralVariantDAO = new StructuralVariantDAO(context);
        this.structuralVariantClusterDAO = new StructuralVariantClusterDAO(context);
        this.structuralVariantFusionDAO = new StructuralVariantFusionDAO(context);
        this.driverCatalogDAO = new DriverCatalogDAO(context);
        this.peachDAO = new PeachDAO(context);
        this.signatureDAO = new SignatureDAO(context);
        this.cuppaDAO = new CuppaDAO(context);
        this.chordDAO = new ChordDAO(context);
        this.virusBreakendDAO = new VirusBreakendDAO(context);
        this.virusInterpreterDAO = new VirusInterpreterDAO(context);
        this.ciderDAO = new CiderDAO(context);
        this.tealDAO = new TealDAO(context);
    }

    public static void addDatabaseCmdLineArgs(final ConfigBuilder configBuilder, boolean isRequired)
    {
        configBuilder.addConfigItem(DB_USER, isRequired, DB_USER_DESC);
        configBuilder.addConfigItem(DB_PASS, isRequired, DB_PASS_DESC);
        configBuilder.addConfigItem(DB_URL, isRequired, DB_URL_DESC);
    }

    public static void addDatabaseCmdLineArgs(@NotNull Options options)
    {
        addDatabaseCmdLineArgs(options, false);
    }

    public static void addDatabaseCmdLineArgs(@NotNull Options options, boolean isRequired)
    {
        options.addOption(Option.builder(DB_USER).desc("Database username").hasArg(true).required(isRequired).build());
        options.addOption(Option.builder(DB_PASS).desc("Database password").hasArg(true).required(isRequired).build());
        options.addOption(Option.builder(DB_URL).desc("Database url").hasArg(true).required(isRequired).build());
    }

    public static boolean hasDatabaseConfig(final ConfigBuilder configBuilder)
    {
        return configBuilder.hasValue(DB_URL) && configBuilder.hasValue(DB_USER) && configBuilder.hasValue(DB_PASS);
    }

    public static DatabaseAccess databaseAccess(final ConfigBuilder configBuilder) throws SQLException
    {
        return databaseAccess(configBuilder, false);
    }

    public static DatabaseAccess databaseAccess(final ConfigBuilder configBuilder, boolean applyDefaultArgs) throws SQLException
    {
        return databaseAccess(
                configBuilder.getValue(DB_USER), configBuilder.getValue(DB_PASS), configBuilder.getValue(DB_URL), applyDefaultArgs);
    }

    public static boolean hasDatabaseConfig(@NotNull CommandLine cmd)
    {
        return cmd.hasOption(DB_URL) && cmd.hasOption(DB_USER) && cmd.hasOption(DB_PASS);
    }

    @NotNull
    public static DatabaseAccess databaseAccess(@NotNull CommandLine cmd) throws SQLException
    {
        return databaseAccess(cmd, false);
    }

    public static DatabaseAccess databaseAccess(@NotNull CommandLine cmd, boolean applyDefaultArgs) throws SQLException
    {
        return databaseAccess(cmd.getOptionValue(DB_USER), cmd.getOptionValue(DB_PASS), cmd.getOptionValue(DB_URL), applyDefaultArgs);
    }

    private static DatabaseAccess databaseAccess(
            final String userName, final String password, final String databaseUrl, boolean applyDefaultArgs) throws SQLException
    {
        String jdbcUrl = "jdbc:" + databaseUrl;

        if(applyDefaultArgs && !jdbcUrl.contains("serverTimezone") && !jdbcUrl.contains("useSSL"))
        {
            jdbcUrl += DB_DEFAULT_ARGS;
        }

        return new DatabaseAccess(userName, password, jdbcUrl);
    }

    @Nullable
    public static DatabaseAccess createDatabaseAccess(final ConfigBuilder configBuilder)
    {
        if(!hasDatabaseConfig(configBuilder))
            return null;

        try
        {
            return databaseAccess(configBuilder, true);
        }
        catch(SQLException e)
        {
            LOGGER.error("DB connection failed: {}", e.toString());
            return null;
        }
    }

    @Nullable
    public static DatabaseAccess createDatabaseAccess(@NotNull CommandLine cmd)
    {
        if(!hasDatabaseConfig(cmd))
            return null;

        try
        {
            return databaseAccess(cmd, true);
        }
        catch(SQLException e)
        {
            LOGGER.error("DB connection failed: {}", e.toString());
            return null;
        }
    }

    @NotNull
    public DSLContext context()
    {
        return context;
    }

    @Override
    public void close()
    {
        try
        {
            connection.close();
        }
        catch(SQLException e)
        {
            LOGGER.error("DB connection close failed: {}", e.toString());
        }
    }

    @Nullable
    private static Settings settings(@NotNull String catalog)
    {
        if(catalog.equals(DEV_CATALOG))
        {
            return null;
        }

        return new Settings().withRenderMapping(new RenderMapping().withSchemata(new MappedSchema().withInput(DEV_CATALOG)
                .withOutput(catalog)));
    }

    @NotNull
    public BufferedWriter<VariantContext> germlineVariantWriter(String tumorSample, String referenceSample, String rnaSample)
    {
        return germlineVariantDAO.writer(tumorSample, referenceSample, rnaSample);
    }

    @NotNull
    public List<String> readPurpleSampleList()
    {
        return purityDAO.getSampleIds();
    }

    @NotNull
    public List<String> readPurpleSampleListPassingQC(double minPurity)
    {
        return purityDAO.getSamplesPassingQC(minPurity);
    }

    @Nullable
    public PurityContext readPurityContext(@NotNull String sampleId)
    {
        return purityDAO.readPurityContext(sampleId);
    }

    @NotNull
    public List<PurpleCopyNumber> readCopynumbers(@NotNull String sample)
    {
        return copyNumberDAO.read(sample);
    }

    @NotNull
    public List<GeneCopyNumber> readGeneCopynumbers(@NotNull String sample, @NotNull List<String> genes)
    {
        return geneCopyNumberDAO.readCopyNumbers(sample, genes);
    }

    @NotNull
    public List<GermlineDeletion> readGermlineDeletions(@NotNull String sample)
    {
        return geneCopyNumberDAO.readGermlineDeletions(sample);
    }

    @NotNull
    public List<String> readSomaticVariantSampleList()
    {
        return somaticVariantDAO.getSamplesList();
    }

    @NotNull
    public List<SomaticVariant> readSomaticVariants(@NotNull String sample, VariantType type)
    {
        return somaticVariantDAO.read(sample, type);
    }

    @NotNull
    public List<StructuralVariantData> readStructuralVariantData(@NotNull String sample)
    {
        return structuralVariantDAO.read(sample);
    }

    @NotNull
    public List<LinxCluster> readClusters(@NotNull String sample)
    {
        return structuralVariantClusterDAO.readClusters(sample);
    }

    @NotNull
    public List<LinxSvAnnotation> readSvAnnotations(@NotNull String sample)
    {
        return structuralVariantClusterDAO.readAnnotations(sample);
    }

    @NotNull
    public List<DriverCatalog> readDriverCatalog(@NotNull String sample)
    {
        return driverCatalogDAO.readDriverData(sample);
    }

    @NotNull
    public List<LinxDriver> readSvDriver(@NotNull String sample)
    {
        return structuralVariantClusterDAO.readSvDrivers(sample);
    }

    @NotNull
    public List<SignatureAllocation> readSignatureAllocations(@NotNull String sample)
    {
        return signatureDAO.readAllocations(sample);
    }

    @NotNull
    public List<LinxFusion> readFusions(@NotNull String sample)
    {
        return structuralVariantFusionDAO.readFusions(sample);
    }

    @NotNull
    public List<LinxBreakend> readBreakends(@NotNull String sample)
    {
        return structuralVariantFusionDAO.readBreakends(sample);
    }

    public void writeCanonicalTranscripts(final String refGenomeVersion, final List<GeneData> geneDataList,
            final List<TranscriptData> transcripts)
    {
        canonicalTranscriptDAO.write(refGenomeVersion, geneDataList, transcripts);
    }

    public void writePurity(@NotNull String sampleId, @NotNull PurityContext context, @NotNull PurpleQC checks)
    {
        purityDAO.write(sampleId, context, checks);
    }

    public void writeBestFitPerPurity(@NotNull String sampleId, @NotNull List<FittedPurity> bestFitPerPurity)
    {
        purityDAO.write(sampleId, bestFitPerPurity);
    }

    public void writeCopynumbers(@NotNull String sample, @NotNull List<PurpleCopyNumber> copyNumbers)
    {
        copyNumberDAO.writeCopyNumber(sample, copyNumbers);
    }

    public void writeAmberMapping(@NotNull String sample, List<AmberMapping> mapping)
    {
        amberDAO.writeMapping(sample, mapping);
    }

    public void writeAmberPatients(List<AmberPatient> mapping)
    {
        amberDAO.writePatients(mapping);
    }

    public void writeAmberAnonymous(List<AmberAnonymous> mapping)
    {
        amberDAO.writeAnonymous(mapping);
    }

    public void writeAmberSample(@NotNull AmberSample identity)
    {
        amberDAO.writeIdentity(identity);
    }

    @NotNull
    public List<AmberSample> readAmberSamples()
    {
        return amberDAO.readSamples();
    }

    @NotNull
    public List<AmberAnonymous> readAmberAnonymous()
    {
        return amberDAO.readAnonymous();
    }

    @NotNull
    public List<AmberPatient> readAmberPatients()
    {
        return amberDAO.readPatients();
    }

    public void truncateAmberPatients()
    {
        amberDAO.truncatePatients();
    }

    public void truncateAmberMappings()
    {
        amberDAO.truncateMappings();
    }

    @NotNull
    public BufferedWriter<SomaticVariant> somaticVariantWriter(@NotNull final String sampleId)
    {
        return somaticVariantDAO.writer(sampleId);
    }

    public void writeStructuralVariants(@NotNull String sampleId, @NotNull List<StructuralVariantData> variants)
    {
        structuralVariantDAO.write(sampleId, variants);
    }

    public void writeSvClusters(@NotNull String sample, @NotNull List<LinxCluster> clusters)
    {
        structuralVariantClusterDAO.writeClusters(sample, clusters);
    }

    public void writeSvLinxData(@NotNull String sample, @NotNull List<LinxSvAnnotation> svData)
    {
        structuralVariantClusterDAO.writeSvData(sample, svData);
    }

    public void writeSvLinks(@NotNull String sample, @NotNull List<LinxLink> links)
    {
        structuralVariantClusterDAO.writeLinks(sample, links);
    }

    public void writeSvDrivers(@NotNull String sample, @NotNull List<LinxDriver> drivers)
    {
        structuralVariantClusterDAO.writeDrivers(sample, drivers);
    }

    public void writeGermlineSVs(@NotNull String sample, @NotNull List<LinxGermlineSv> germlineSvs)
    {
        germlineVariantDAO.writeGermlineSVs(sample, germlineSvs);
    }

    public void writeGermlineBreakends(@NotNull String sample, @NotNull List<LinxBreakend> germlineBreakends)
    {
        germlineVariantDAO.writeGermlineBreakends(sample, germlineBreakends);
    }

    public void writeSignatures(@NotNull String sample, @NotNull List<SignatureAllocation> sigAllocations)
    {
        signatureDAO.write(sample, sigAllocations);
    }

    public void writeGeneCopyNumbers(@NotNull String sample, @NotNull List<GeneCopyNumber> geneCopyNumbers)
    {
        geneCopyNumberDAO.writeCopyNumber(sample, geneCopyNumbers);
    }

    public void writeGermlineDeletions(@NotNull String sample, @NotNull List<GermlineDeletion> deletions)
    {
        geneCopyNumberDAO.writeGermlineDeletions(sample, deletions);
    }

    public void writeLinxDriverCatalog(@NotNull String sample, @NotNull List<DriverCatalog> driverCatalog,
            final EnumSet<DriverType> driverTypes)
    {
        driverCatalogDAO.writeLinxDrivers(sample, driverCatalog, driverTypes);
    }

    public void writePurpleDriverCatalog(@NotNull String sample, @Nullable List<DriverCatalog> somaticCatalog,
            @Nullable List<DriverCatalog> germlineCatalog)
    {
        driverCatalogDAO.writePurpleDrivers(sample, somaticCatalog, germlineCatalog);
    }

    public void writeMetrics(@NotNull String sample, @NotNull WGSMetricWithQC metrics)
    {
        metricDAO.writeMetrics(sample, metrics);
    }

    public void writeFlagstats(@NotNull String sample, @NotNull Flagstat refFlagstat, @NotNull Flagstat tumorFlagstat)
    {
        flagstatDAO.writeFlagstats(sample, refFlagstat, tumorFlagstat);
    }

    public void writePeach(@NotNull String sample, @NotNull List<PeachGenotype> peachGenotypes)
    {
        peachDAO.writePeach(sample, peachGenotypes);
    }

    public void writeCuppa(@NotNull String sample, @NotNull CuppaPredictions cuppaPredictions, int topNProbs) throws IOException {
        cuppaDAO.writeCuppa2(sample, cuppaPredictions, topNProbs);
    }

    public void writeVirusBreakend(@NotNull String sample, @NotNull List<VirusBreakend> virusBreakends)
    {
        virusBreakendDAO.writeVirusBreakend(sample, virusBreakends);
    }

    public void writeVirusInterpreter(@NotNull String sample, @NotNull List<AnnotatedVirus> virusAnnotations)
    {
        virusInterpreterDAO.writeVirusInterpreter(sample, virusAnnotations);
    }

    public void writeChord(@NotNull String sample, @NotNull ChordData chordData)
    {
        chordDAO.writeChord(sample, chordData);
    }

    public ChordData readChord(final String sampleId)
    {
        return chordDAO.readChord(sampleId);
    }

    public void writeSnpCheck(@NotNull String sample, boolean isPass)
    {
        snpCheckDAO.write(sample, isPass);
    }

    public void writeGenePanel(@NotNull final List<DriverGene> driverGenes)
    {
        driverGenePanelDAO.writeDriverGenes(driverGenes);
    }

    public void writeCdr3Sequences(@NotNull String sample, @NotNull List<Cdr3Sequence> cdr3Sequences)
    {
        ciderDAO.writeCdr3Sequence(sample, cdr3Sequences);
    }

    public void writeCdr3LocusSummaries(@NotNull String sample, @NotNull List<Cdr3LocusSummary> locusSummaries)
    {
        ciderDAO.writeLocusSummaries(sample, locusSummaries);
    }

    public void writeTelomereLength(@NotNull String sample, @Nullable TelomereLength germlineTelomereLength, @Nullable TelomereLength somaticTelomereLength)
    {
        tealDAO.writeTelomereLength(sample, germlineTelomereLength, somaticTelomereLength);
    }

    public void deletePipelineDataForSample(@NotNull String sample)
    {
        LOGGER.info("Deleting metric data for sample: {}", sample);
        metricDAO.deleteMetricForSample(sample);

        LOGGER.info("Deleting flagstat data for sample: {}", sample);
        flagstatDAO.deleteFlagstatsForSample(sample);

        LOGGER.info("Deleting AMBER data for sample: {}", sample);
        amberDAO.deleteAmberRecordsForSample(sample);

        LOGGER.info("Deleting purity data for sample: {}", sample);
        purityDAO.deletePurityForSample(sample);

        LOGGER.info("Deleting copy number data for sample: {}", sample);
        copyNumberDAO.deleteCopyNumberForSample(sample);

        LOGGER.info("Deleting gene copy numbers for sample: {}", sample);
        geneCopyNumberDAO.deleteGeneCopyNumberForSample(sample);
        geneCopyNumberDAO.deleteGermlineDeletionsForSample(sample);

        LOGGER.info("Deleting somatic variants for sample: {}", sample);
        somaticVariantDAO.deleteSomaticVariantForSample(sample);

        LOGGER.info("Deleting germline variant data for sample: {}", sample);
        germlineVariantDAO.deleteGermlineVariantsForSample(sample);
        germlineVariantDAO.deleteGermlineStructuralVariantsForSample(sample);

        LOGGER.info("Deleting structural variant annotation data for sample: {}", sample);
        structuralVariantFusionDAO.deleteAnnotationsForSample(sample);

        LOGGER.info("Deleting structural variant cluster data for sample: {}", sample);
        structuralVariantClusterDAO.deleteClusterDataForSample(sample);

        LOGGER.info("Deleting structural variants for sample: {}", sample);
        structuralVariantDAO.deleteStructuralVariantsForSample(sample);

        LOGGER.info("Deleting signatures for sample: {}", sample);
        signatureDAO.deleteSignatureDataForSample(sample);

        LOGGER.info("Deleting driver catalog for sample: {}", sample);
        driverCatalogDAO.deleteForSample(sample);

        LOGGER.info("Deleting CHORD data for sample: {}", sample);
        chordDAO.deleteChordForSample(sample);

        LOGGER.info("Deleting PEACH data for sample: {}", sample);
        peachDAO.deletePeachForSample(sample);

        LOGGER.info("Deleting virus breakend data for sample: {}", sample);
        virusBreakendDAO.deleteVirusBreakendForSample(sample);

        LOGGER.info("Deleting virus annotation data for sample: {}", sample);
        virusInterpreterDAO.deleteVirusAnnotationForSample(sample);

        LOGGER.info("Deleting CUPPA result for sample: {}", sample);
        cuppaDAO.deleteCuppaForSample(sample);

        LOGGER.info("Deleting CIDER data for sample: {}", sample);
        ciderDAO.deleteCiderDataForSample(sample);

        LOGGER.info("Deleting TEAL data for sample: {}", sample);
        tealDAO.deleteTealDataForSample(sample);

        LOGGER.info("All data for sample '{}' has been deleted", sample);
    }
}

