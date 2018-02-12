package com.hartwig.hmftools.svannotation;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;
import com.hartwig.hmftools.common.cosmic.fusions.CosmicFusionModel;
import com.hartwig.hmftools.common.cosmic.fusions.CosmicFusions;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.variant.structural.*;
import com.hartwig.hmftools.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.svannotation.analysis.*;
import com.hartwig.hmftools.svannotation.dao.StructuralVariantAnnotationDAO;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
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

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    public static void main(@NotNull final String[] args) throws ParseException, IOException, HartwigException, SQLException {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final DatabaseAccess dbAccess = databaseAccess(cmd);

        boolean loadFromDB = cmd.hasOption(LOAD_FROM_DB);
        final String tumorSample = cmd.getOptionValue(SAMPLE);
        boolean runClustering = cmd.hasOption(CLUSTER_SVS);

        if(!loadFromDB) {

            boolean skipAnnotations = cmd.hasOption(SKIP_ANNOTATIONS);
            LOGGER.info("reading VCF File");
            final List<StructuralVariant> variants = readFromVcf(cmd.getOptionValue(VCF_FILE));

            LOGGER.info("enriching structural variants based on purple data");
            final List<EnrichedStructuralVariant> enrichedVariants = enrichStructuralVariants(variants, dbAccess, tumorSample);

            LOGGER.info("initialising MqSql annotator");
            final VariantAnnotator annotator = MySQLAnnotator.make("jdbc:" + cmd.getOptionValue(ENSEMBL_DB));

            LOGGER.info("loading Cosmic Fusion data");
            final CosmicFusionModel cosmicGeneFusions = CosmicFusions.readFromCSV(cmd.getOptionValue(FUSION_CSV));

            final StructuralVariantAnalyzer analyzer = new StructuralVariantAnalyzer(annotator, HmfGenePanelSupplier.hmfPanelGeneList(), cosmicGeneFusions);
            LOGGER.info("analyzing structural variants for impact via disruptions and fusions");
            final StructuralVariantAnalysis analysis = analyzer.run(enrichedVariants, skipAnnotations);

            if(runClustering)
            {
                LOGGER.info("running clustering logic");

                SvClusteringConfig clusteringConfig = new SvClusteringConfig();
                clusteringConfig.setOutputCsvFile(tumorSample);
                clusteringConfig.setBaseDistance(Integer.parseInt(cmd.getOptionValue(CLUSTER_BASE_DISTANCE, "0")));
                StructuralVariantClustering svClusterer = new StructuralVariantClustering(clusteringConfig);

                svClusterer.loadFromEnrichedSVs(tumorSample, enrichedVariants);
                svClusterer.runClustering();
            }

            LOGGER.info("persisting variants and annotations to database");
            dbAccess.writeStructuralVariants(tumorSample, enrichedVariants);

            final StructuralVariantAnnotationDAO annotationDAO = new StructuralVariantAnnotationDAO(dbAccess.getContext());
            annotationDAO.write(analysis);
        }
        else
        {
            LOGGER.info("running clustering logic");

            SvClusteringConfig clusteringConfig = new SvClusteringConfig();
            clusteringConfig.setOutputCsvFile(tumorSample);
            clusteringConfig.setBaseDistance(Integer.parseInt(cmd.getOptionValue(CLUSTER_BASE_DISTANCE, "0")));
            StructuralVariantClustering svClusterer = new StructuralVariantClustering(clusteringConfig);

            List<SvClusterData> svClusterData = queryStructuralVariantData(dbAccess, tumorSample);
            svClusterer.loadFromDatabase(tumorSample, svClusterData);
            svClusterer.runClustering();
        }

        LOGGER.info("run complete");
    }

    @NotNull
    private static List<StructuralVariant> readFromVcf(@NotNull String vcfFileLocation) throws IOException {
        final StructuralVariantFactory factory = new StructuralVariantFactory();
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

    private static List<SvClusterData> queryStructuralVariantData(DatabaseAccess dbAccess, final String sampleId)
    {
        List<SvClusterData> svClusterDataItems = Lists.newArrayList();

        List<StructuralVariantData> svRecords = dbAccess.readStructuralVariants(sampleId);

        for(final StructuralVariantData svRecord : svRecords)
        {
            svClusterDataItems.add(SvClusterData.from(svRecord));
        }

        return svClusterDataItems;
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
        options.addOption(SKIP_ANNOTATIONS, false, "Skip annotations, including Ensemble DB data sync, for testing only)");
        options.addOption(CLUSTER_BASE_DISTANCE, true, "Clustering base distance, defaults to 1000");
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
        final String databaseUrl = cmd.getOptionValue(DB_URL);  //e.g. mysql://localhost:port/database";
        final String jdbcUrl = "jdbc:" + databaseUrl;
        return new DatabaseAccess(userName, password, jdbcUrl);
    }
}
