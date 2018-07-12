package com.hartwig.hmftools.svannotation;

import static com.hartwig.hmftools.svanalysis.annotators.SvPONAnnotator.PON_FILTER_PASS;
import static com.hartwig.hmftools.svanalysis.annotators.SvPONAnnotator.PON_FILTER_PON;

import java.io.FileInputStream;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;
import com.hartwig.hmftools.common.fusions.KnownFusionsModel;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariantFactory;
import com.hartwig.hmftools.common.variant.structural.ImmutableEnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory;
import com.hartwig.hmftools.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.svanalysis.annotators.SvPONAnnotator;
import com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalysis;
import com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalyzer;
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
    private static final String ENSEMBL_DB_LOCAL = "local_ensembl";
    private static final String FUSION_PAIRS_CSV = "fusion_pairs_csv";
    private static final String PROMISCUOUS_FIVE_CSV = "promiscuous_five_csv";
    private static final String PROMISCUOUS_THREE_CSV = "promiscuous_three_csv";

    private static final String SOURCE_SVS_FROM_DB = "source_svs_from_db";
    private static final String LOG_DEBUG = "log_debug";
    private static final String SV_PON_FILE = "sv_pon_file";

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    public static void main(@NotNull final String[] args) throws ParseException, IOException, SQLException {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        final String tumorSample = cmd.getOptionValue(SAMPLE);

        LOGGER.info("annotating variants for sample({})", tumorSample);

        final DatabaseAccess dbAccess = databaseAccess(cmd);

        if (cmd.hasOption(LOG_DEBUG)) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        SvPONAnnotator svPONAnnotator = null;

        if (cmd.hasOption(SV_PON_FILE)) {
            svPONAnnotator = new SvPONAnnotator();
            svPONAnnotator.loadPonFile(cmd.getOptionValue(SV_PON_FILE));
        }

        boolean sourceSVsFromDB = cmd.hasOption(SOURCE_SVS_FROM_DB);

        List<EnrichedStructuralVariant> svList;

        if (sourceSVsFromDB) {
            svList = dbAccess.readStructuralVariants(tumorSample);

            if (svList.isEmpty()) {
                LOGGER.info("no SVs loaded from DB");
                return;
            }
        } else {
            LOGGER.info("reading VCF File");
            final List<StructuralVariant> variants = readFromVcf(cmd.getOptionValue(VCF_FILE), true);

            LOGGER.info("enriching structural variants based on purple data");
            svList = enrichStructuralVariants(variants, dbAccess, tumorSample);
        }

        List<EnrichedStructuralVariant> updatedSVs = Lists.newArrayList();

        if (svPONAnnotator != null && svPONAnnotator.hasEntries()) {

            for (EnrichedStructuralVariant variant : svList) {
                int ponCount = svPONAnnotator.getPonOccurenceCount(variant.chromosome(true),
                        variant.chromosome(false),
                        variant.position(true),
                        variant.position(false),
                        variant.orientation(true),
                        variant.orientation(false),
                        variant.type().toString());

                String filter = ponCount > 1 ? PON_FILTER_PON : PON_FILTER_PASS;

                final ImmutableEnrichedStructuralVariant updatedSV =
                        ImmutableEnrichedStructuralVariant.builder().from(variant).filter(filter).build();
                updatedSVs.add(updatedSV);
            }
        } else {
            updatedSVs = svList;
        }

        LOGGER.info("persisting {} SVs to database", updatedSVs.size());

        dbAccess.writeStructuralVariants(tumorSample, updatedSVs);

        if (!sourceSVsFromDB) {
            // NEVA: We read after we write to populate the primaryId field
            final List<EnrichedStructuralVariant> enrichedVariants = dbAccess.readStructuralVariants(tumorSample);

            DatabaseAccess ensembleDBConn = null;

            if (cmd.hasOption(ENSEMBL_DB_LOCAL)) {

                // CHSH: The same credential work for both hmf and ensembl DBs
                final String ensembleJdbcUrl = "jdbc:" + cmd.getOptionValue(ENSEMBL_DB);
                final String ensembleUser = cmd.getOptionValue(DB_USER);
                final String ensemblePassword = cmd.getOptionValue(DB_PASS);

                LOGGER.debug("connecting to local ensembl DB: {}", cmd.getOptionValue(ENSEMBL_DB));

                try {
                    ensembleDBConn = new DatabaseAccess(ensembleUser, ensemblePassword, ensembleJdbcUrl);
                } catch (Exception e) {
                    LOGGER.warn("Ensembl DB connection failed: {}", e.toString(), e.getMessage());
                    return;
                }
            }

            final VariantAnnotator annotator = ensembleDBConn != null
                    ? new MySQLAnnotator(ensembleDBConn.context())
                    : MySQLAnnotator.make("jdbc:" + cmd.getOptionValue(ENSEMBL_DB));

            LOGGER.info("loading Fusion data");
            final KnownFusionsModel knownFusionsModel =
                    KnownFusionsModel.fromInputStreams(new FileInputStream(cmd.getOptionValue(FUSION_PAIRS_CSV)),
                            new FileInputStream(cmd.getOptionValue(PROMISCUOUS_FIVE_CSV)),
                            new FileInputStream(cmd.getOptionValue(PROMISCUOUS_THREE_CSV)));

            final StructuralVariantAnalyzer analyzer =
                    new StructuralVariantAnalyzer(annotator, HmfGenePanelSupplier.hmfPanelGeneList(), knownFusionsModel);
            LOGGER.info("analyzing structural variants for impact via disruptions and fusions");
            final StructuralVariantAnalysis analysis = analyzer.run(enrichedVariants);

            LOGGER.info("persisting annotations to database");
            final StructuralVariantAnnotationDAO annotationDAO = new StructuralVariantAnnotationDAO(dbAccess.context());
            annotationDAO.write(analysis);
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
    private static Options createBasicOptions() {
        final Options options = new Options();
        options.addOption(VCF_FILE, true, "Path to the vcf file.");
        options.addOption(SAMPLE, true, "Tumor sample.");
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");
        options.addOption(FUSION_PAIRS_CSV, true, "Path towards a CSV containing white-listed gene fusion pairs.");
        options.addOption(PROMISCUOUS_FIVE_CSV, true, "Path towards a CSV containing white-listed promiscuous 5' genes.");
        options.addOption(PROMISCUOUS_THREE_CSV, true, "Path towards a CSV containing white-listed promiscuous 3' genes.");
        options.addOption(ENSEMBL_DB, true, "Annotate structural variants using this Ensembl DB URI");
        options.addOption(ENSEMBL_DB_LOCAL, false, "Connect to local Ensembl DB");
        options.addOption(SOURCE_SVS_FROM_DB, false, "Skip annotations, including Ensemble DB data sync, for testing only)");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");
        options.addOption(SV_PON_FILE, true, "PON file for SVs");
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
