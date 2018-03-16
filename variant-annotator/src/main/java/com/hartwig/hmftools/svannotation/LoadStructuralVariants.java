package com.hartwig.hmftools.svannotation;

import java.io.IOException;
import java.sql.SQLException;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;
import com.hartwig.hmftools.common.cosmic.fusions.CosmicFusionModel;
import com.hartwig.hmftools.common.cosmic.fusions.CosmicFusions;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariantFactory;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory;
import com.hartwig.hmftools.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalysis;
import com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalyzer;
import com.hartwig.hmftools.svannotation.analysis.StructuralVariantClustering;
import com.hartwig.hmftools.svannotation.analysis.SvClusterData;
import com.hartwig.hmftools.svannotation.analysis.SvClusteringConfig;
import com.hartwig.hmftools.svannotation.dao.StructuralVariantAnnotationDAO;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class LoadStructuralVariants {

    private static final Logger LOGGER = LogManager.getLogger(LoadStructuralVariants.class);

    private static final String SAMPLE = "sample";
    private static final String VCF_FILE = "vcf_file";
    private static final String ENSEMBL_DB = "ensembl_db";
    private static final String FUSION_CSV = "fusion_csv";
    private static final String CLUSTER_SVS = "cluster_svs";
    private static final String SKIP_ANNOTATIONS = "skip_annotations";
    private static final String LOAD_FROM_DB = "load_from_db";
    private static final String CLUSTER_BASE_DISTANCE = "cluster_bases";
    private static final String DATA_OUTPUT_PATH = "data_output_path";
    private static final String LOG_DEBUG = "log_debug";
    private static final String WRITE_FILTERED_SVS = "log_filters";
    private static final String SV_PON_FILE = "sv_pon_file";
    private static final String FRAGILE_SITE_FILE = "fragile_site_file";
    private static final String LINE_ELEMENT_FILE = "line_element_file";

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    public static void main(@NotNull final String[] args) throws ParseException, IOException, SQLException {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        boolean loadFromDB = cmd.hasOption(LOAD_FROM_DB);
        final String tumorSample = cmd.getOptionValue(SAMPLE);
        boolean runClustering = cmd.hasOption(CLUSTER_SVS);
        boolean createFilteredPON = cmd.hasOption(WRITE_FILTERED_SVS);

        final DatabaseAccess dbAccess = !createFilteredPON ? databaseAccess(cmd) : null;

        if (cmd.hasOption(LOG_DEBUG)) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        StructuralVariantClustering svClusterer = null;

        if (runClustering) {
            LOGGER.info("will run clustering logic");

            SvClusteringConfig clusteringConfig = new SvClusteringConfig();
            clusteringConfig.setOutputCsvPath(cmd.getOptionValue(DATA_OUTPUT_PATH));
            clusteringConfig.setBaseDistance(Integer.parseInt(cmd.getOptionValue(CLUSTER_BASE_DISTANCE, "0")));
            clusteringConfig.setUseCombinedOutputFile(tumorSample.equals("*"));
            clusteringConfig.setSvPONFile(cmd.getOptionValue(SV_PON_FILE));
            clusteringConfig.setFragileSiteFile(cmd.getOptionValue(FRAGILE_SITE_FILE));
            clusteringConfig.setLineElementFile(cmd.getOptionValue(LINE_ELEMENT_FILE));
            svClusterer = new StructuralVariantClustering(clusteringConfig);
        }

        if (createFilteredPON) {
            LOGGER.info("reading VCF file including filtered SVs");

            FilteredSVWriter filteredSvWriter = new FilteredSVWriter(cmd.getOptionValue(VCF_FILE), cmd.getOptionValue(DATA_OUTPUT_PATH));
            filteredSvWriter.processVcfFiles();

            LOGGER.info("reads complete");
            return;
        }

        if (!loadFromDB) {
            boolean skipAnnotations = cmd.hasOption(SKIP_ANNOTATIONS);
            LOGGER.info("reading VCF File");
            final List<StructuralVariant> variants = readFromVcf(cmd.getOptionValue(VCF_FILE), true);

            LOGGER.info("enriching structural variants based on purple data");
            final List<EnrichedStructuralVariant> enrichedVariantWithoutPrimaryId =
                    enrichStructuralVariants(variants, dbAccess, tumorSample);

            LOGGER.info("persisting variants to database");
            dbAccess.writeStructuralVariants(tumorSample, enrichedVariantWithoutPrimaryId);

            // NEVA: We read after we write to populate the primaryId field
            final List<EnrichedStructuralVariant> enrichedVariants = dbAccess.readStructuralVariants(tumorSample);

            LOGGER.info("initialising MqSql annotator");
            final VariantAnnotator annotator = MySQLAnnotator.make("jdbc:" + cmd.getOptionValue(ENSEMBL_DB));

            LOGGER.info("loading Cosmic Fusion data");
            final CosmicFusionModel cosmicGeneFusions = CosmicFusions.readFromCSV(cmd.getOptionValue(FUSION_CSV));

            final StructuralVariantAnalyzer analyzer =
                    new StructuralVariantAnalyzer(annotator, HmfGenePanelSupplier.hmfPanelGeneList(), cosmicGeneFusions);
            LOGGER.info("analyzing structural variants for impact via disruptions and fusions");
            final StructuralVariantAnalysis analysis = analyzer.run(enrichedVariants, skipAnnotations);

            if (runClustering) {
                svClusterer.loadFromEnrichedSVs(tumorSample, enrichedVariants);
                svClusterer.runClustering();
            }

            LOGGER.info("persisting annotations to database");
            final StructuralVariantAnnotationDAO annotationDAO = new StructuralVariantAnnotationDAO(dbAccess.getContext());
            annotationDAO.write(analysis);
        } else {
            // CHSH: only run the SV clustering routines, taking source data from the SV table(s)
            // KODU: Below assert feels somewhat risky!?
            assert runClustering;

            List<String> samplesList = Lists.newArrayList();

            if (tumorSample.isEmpty() || tumorSample.equals("*")) {
                samplesList = getStructuralVariantSamplesList(dbAccess);
            } else if (tumorSample.contains(",")) {
                String[] tumorList = tumorSample.split(",");
                samplesList = Arrays.stream(tumorList).collect(Collectors.toList());
            } else {
                samplesList.add(tumorSample);
            }

            int count = 0;
            for (final String sample : samplesList) {
                ++count;
                LOGGER.info("clustering for sample({}), total({})", sample, count);

                List<SvClusterData> svClusterData = queryStructuralVariantData(dbAccess, sample);

                svClusterer.loadFromDatabase(sample, svClusterData);
                svClusterer.runClustering();
            }
        }

        LOGGER.info("run complete");
    }

    @NotNull
    private static List<StructuralVariant> readFromVcf(@NotNull String vcfFileLocation, final boolean filterOnPasses) throws IOException {
        final StructuralVariantFactory factory = new StructuralVariantFactory(filterOnPasses);
        try (final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(vcfFileLocation,
                new VCFCodec(),
                false)) {
            reader.iterator().forEach(factory::addVariantContext);
        }
        return factory.results();
    }

    @NotNull
    private static List<EnrichedStructuralVariant> enrichStructuralVariants(@NotNull List<StructuralVariant> variants,
            @NotNull DatabaseAccess dbAccess, @NotNull String tumorSample) {
        final PurityContext purityContext = dbAccess.readPurityContext(tumorSample);

        if (purityContext == null) {
            LOGGER.warn("Unable to retrieve purple data. Enrichment may be incomplete.");
        }

        final PurityAdjuster purityAdjuster = purityContext == null
                ? new PurityAdjuster(Gender.FEMALE, 1, 1)
                : new PurityAdjuster(purityContext.gender(), purityContext.bestFit().purity(), purityContext.bestFit().normFactor());

        final List<PurpleCopyNumber> copyNumberList = dbAccess.readCopynumbers(tumorSample);
        final Multimap<String, PurpleCopyNumber> copyNumbers = Multimaps.index(copyNumberList, PurpleCopyNumber::chromosome);
        return EnrichedStructuralVariantFactory.enrich(variants, purityAdjuster, copyNumbers);
    }

    @NotNull
    private static List<SvClusterData> queryStructuralVariantData(@NotNull DatabaseAccess dbAccess, @NotNull String sampleId) {
        List<SvClusterData> svClusterDataItems = Lists.newArrayList();

        List<StructuralVariantData> svRecords = dbAccess.readStructuralVariantData(sampleId);

        for (final StructuralVariantData svRecord : svRecords) {
            svClusterDataItems.add(new SvClusterData(svRecord));
        }

        return svClusterDataItems;
    }

    @NotNull
    private static List<String> getStructuralVariantSamplesList(@NotNull DatabaseAccess dbAccess) {
        return dbAccess.getStructuralVariantSampleList("");
    }

    @NotNull
    private static Options createBasicOptions() {
        final Options options = new Options();
        options.addOption(VCF_FILE, true, "Path to the vcf file.");
        options.addOption(SAMPLE, true, "Tumor sample.");
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");
        options.addOption(FUSION_CSV, true, "Path towards a CSV containing white-listed gene fusions.");
        options.addOption(ENSEMBL_DB, true, "Annotate structural variants using this Ensembl DB URI");
        options.addOption(LOAD_FROM_DB, false, "Use existing SV data, skip loading from VCF");
        options.addOption(CLUSTER_SVS, false, "Whether to run clustering logic");
        options.addOption(DATA_OUTPUT_PATH, true, "CSV output directory");
        options.addOption(SKIP_ANNOTATIONS, false, "Skip annotations, including Ensemble DB data sync, for testing only)");
        options.addOption(CLUSTER_BASE_DISTANCE, true, "Clustering base distance, defaults to 1000");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");
        options.addOption(WRITE_FILTERED_SVS, false, "Includes filtered SVs and writes all to file for PON creation");
        options.addOption(SV_PON_FILE, true, "PON file for SVs");
        options.addOption(LINE_ELEMENT_FILE, true, "Line Elements file for SVs");
        options.addOption(FRAGILE_SITE_FILE, true, "Fragile Site file for SVs");

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static DatabaseAccess databaseAccess(@NotNull final CommandLine cmd) throws SQLException {
        final String userName = cmd.getOptionValue(DB_USER);
        final String password = cmd.getOptionValue(DB_PASS);
        final String databaseUrl = cmd.getOptionValue(DB_URL);
        final String jdbcUrl = "jdbc:" + databaseUrl;
        return new DatabaseAccess(userName, password, jdbcUrl);
    }
}
