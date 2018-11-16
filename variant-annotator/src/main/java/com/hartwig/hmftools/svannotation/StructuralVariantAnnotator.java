package com.hartwig.hmftools.svannotation;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.variant.structural.annotation.SvPONAnnotator.PON_FILTER_PASS;
import static com.hartwig.hmftools.common.variant.structural.annotation.SvPONAnnotator.PON_FILTER_PON;
import static com.hartwig.hmftools.svannotation.analysis.SvFusionAnalyser.FUSION_PAIRS_CSV;
import static com.hartwig.hmftools.svannotation.analysis.SvFusionAnalyser.PROMISCUOUS_FIVE_CSV;
import static com.hartwig.hmftools.svannotation.analysis.SvFusionAnalyser.PROMISCUOUS_THREE_CSV;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.dnds.DndsDriverGeneLikelihoodSupplier;
import com.hartwig.hmftools.common.fusions.KnownFusionsModel;
import com.hartwig.hmftools.common.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariantFactory;
import com.hartwig.hmftools.common.variant.structural.ImmutableEnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneDisruption;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.StructuralVariantAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.SvGeneTranscriptCollection;
import com.hartwig.hmftools.common.variant.structural.annotation.SvPONAnnotator;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.svannotation.analysis.ImmutableStructuralVariantAnalysis;
import com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalysis;
import com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalyzer;
import com.hartwig.hmftools.svannotation.analysis.SvFusionAnalyser;
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

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;
import htsjdk.variant.vcf.VCFCodec;

public class StructuralVariantAnnotator
{
    private static final String SAMPLE = "sample";
    private static final String VCF_FILE = "vcf_file";
    private static final String ENSEMBL_DB = "ensembl_db";
    private static final String ENSEMBL_DB_USER = "ensembl_user";
    private static final String ENSEMBL_DB_PASS = "ensembl_pass";
    private static final String REF_GENOME = "ref_genome";
    private static final String DATA_OUTPUT_DIR = "data_output_dir";

    private static final String SOURCE_SVS_FROM_DB = "source_svs_from_db";
    private static final String LOAD_ANNOTATIONS_FROM_FILE = "load_annotations";
    private static final String SKIP_DB_UPLOAD = "skip_db_upload";
    private static final String LOG_DEBUG = "log_debug";
    private static final String SV_PON_FILE = "sv_pon_file";

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    private static final String TEST_SV_LIMIT = "test_limit_sv_count";

    private String mSampleId;
    private String mDataPath;

    private final CommandLine mCmdLineArgs;
    private DatabaseAccess mDbAccess;
    private boolean mSourceSvFromDB;
    private SvGeneTranscriptCollection mSvGeneTranscriptCollection;

    // Let PON filtered SVs through since GRIDSS PON filtering is performed upstream
    private static Set<String> ALLOWED_FILTERS = Sets.newHashSet( "INFERRED", PON_FILTER_PON, PON_FILTER_PASS);

    private static final Logger LOGGER = LogManager.getLogger(StructuralVariantAnnotator.class);

    public StructuralVariantAnnotator(final CommandLine cmd)
    {
        mSampleId = "";
        mDataPath = "";
        mDbAccess = null;
        mSourceSvFromDB = false;
        mCmdLineArgs = cmd;
        mSvGeneTranscriptCollection = new SvGeneTranscriptCollection();
    }

    public boolean initialise()
    {
        mSampleId = mCmdLineArgs.getOptionValue(SAMPLE);

        try
        {
            mDbAccess = databaseAccess(mCmdLineArgs);
        }
        catch (SQLException e)
        {
            LOGGER.error("database connection failed: {}", e.toString());
            return false;
        }

        mSourceSvFromDB = mCmdLineArgs.hasOption(SOURCE_SVS_FROM_DB);

        mDataPath = mCmdLineArgs.hasOption(DATA_OUTPUT_DIR) ? mCmdLineArgs.getOptionValue(DATA_OUTPUT_DIR) : "";
        mSvGeneTranscriptCollection.setDataPath(mDataPath);

        return true;
    }

    public boolean run()
    {
        // create the annotator - typically a MySQL connection to the Ensembl database
        VariantAnnotator annotator = null;

        if(mCmdLineArgs.hasOption(ENSEMBL_DB))
        {
            try
            {
                final String ensembleUrl = "jdbc:" + mCmdLineArgs.getOptionValue(ENSEMBL_DB);
                final String ensembleUser = mCmdLineArgs.getOptionValue(ENSEMBL_DB_USER, "");
                final String ensemblePassword = mCmdLineArgs.getOptionValue(ENSEMBL_DB_PASS, "");

                LOGGER.debug("connecting to ensembl DB: {}", ensembleUrl);

                if (!ensembleUser.isEmpty() && !ensemblePassword.isEmpty())
                {
                    DatabaseAccess ensembleDBConn = new DatabaseAccess(ensembleUser, ensemblePassword, ensembleUrl);
                    annotator = new MySQLAnnotator(ensembleDBConn.context());
                }
                else
                {
                    MySQLAnnotator.make(ensembleUrl);
                }
            }
            catch (SQLException e)
            {
                LOGGER.warn("Ensembl DB connection failed: {}", e.toString());
                return false;
            }
        }
        else
        {
            annotator = new NullAnnotator();
        }

        LOGGER.debug("loading known fusion data");
        KnownFusionsModel knownFusionsModel = null;

        try
        {
            knownFusionsModel = KnownFusionsModel.fromInputStreams(
                    new FileInputStream(mCmdLineArgs.getOptionValue(FUSION_PAIRS_CSV)),
                    new FileInputStream(mCmdLineArgs.getOptionValue(PROMISCUOUS_FIVE_CSV)),
                    new FileInputStream(mCmdLineArgs.getOptionValue(PROMISCUOUS_THREE_CSV)));

        }
        catch(IOException e)
        {
            LOGGER.error("failed to load known fusion files");
            return false;
        }

        // TODO (KODU): Add potentially actionable genes, see also instantiation in patient reporter.
        final StructuralVariantAnalyzer svAnalyser = new StructuralVariantAnalyzer(annotator, tsgDriverGeneIDs(), knownFusionsModel);

        if(mSampleId.isEmpty() || mSampleId.equals("*"))
        {
            List<String> samplesList = mDbAccess.structuralVariantSampleList("");

            LOGGER.info("loaded {} samples from database", samplesList.size());

            for(final String sampleId : samplesList)
            {
                runSample(sampleId, svAnalyser);
            }
        }
        else
        {
            runSample(mSampleId, svAnalyser);
        }

        return true;
    }

    private void runSample(final String sampleId, final StructuralVariantAnalyzer svAnalyser)
    {
        LOGGER.info("annotating variants for sample({})", sampleId);

        List<EnrichedStructuralVariant> enrichedVariants;

        if(mSourceSvFromDB)
        {
            // optionally load existing SVs from the database rather than from VCF
            enrichedVariants = mDbAccess.readStructuralVariants(sampleId)
                    .stream().filter(x -> x.filter().equals("PASS")).collect(toList());

            if (enrichedVariants.isEmpty())
            {
                LOGGER.debug("sample({}) no SVs loaded from DB", sampleId);
                return;
            }

            LOGGER.debug("sample({}) loaded {} SVs from DB", sampleId, enrichedVariants.size());
        }
        else
        {
            List<EnrichedStructuralVariant> svList = loadSVsFromVCF(sampleId);

            LOGGER.info("sample({}) persisting {} SVs to database", sampleId, svList.size());

            mDbAccess.writeStructuralVariants(sampleId, svList);

            // re-read the data to get primaryId field as a foreign key for disruptions and fusions
            enrichedVariants = mDbAccess.readStructuralVariants(sampleId);;
        }

        List<StructuralVariantAnnotation> annotations;
        if(mCmdLineArgs.hasOption(LOAD_ANNOTATIONS_FROM_FILE) && !mDataPath.isEmpty())
        {
            mSvGeneTranscriptCollection.loadSampleGeneTranscripts(sampleId);

            annotations = createAnnotations(enrichedVariants);

            LOGGER.debug("loaded {} Ensembl annotations from file", annotations.size());
        }
        else if(svAnalyser.getGeneDataAnnotator() != null)
        {
            int testSvLimit = Integer.parseInt(mCmdLineArgs.getOptionValue(TEST_SV_LIMIT, "0"));

            if(testSvLimit > 0)
            {
                List<EnrichedStructuralVariant> subset = Lists.newArrayList(enrichedVariants.subList(0, testSvLimit));
                annotations = svAnalyser.findAnnotations(subset);
            }
            else
            {
                annotations = svAnalyser.findAnnotations(enrichedVariants);
            }

            LOGGER.debug("sample({}) matched {} annotations from Ensembl database", sampleId, annotations.size());

            // optionally persist to save having to look up this Ensembl data again
            if(!mDataPath.isEmpty())
                mSvGeneTranscriptCollection.writeAnnotations(sampleId, annotations);
        }
        else
        {
            LOGGER.error("Ensemble data not loaded from DB nor file");
            return;
        }

        if(!mCmdLineArgs.hasOption(SKIP_DB_UPLOAD))
        {
            LOGGER.debug("sample({}) finding disruptions and fusions", sampleId);

            final List<GeneFusion> fusions = svAnalyser.findFusions(annotations);
            final List<GeneDisruption> disruptions = svAnalyser.findDisruptions(annotations);

            LOGGER.debug("sample({}) found {} disruptions and {} fusions", sampleId, disruptions.size(), fusions.size());

            final StructuralVariantAnalysis analysis = ImmutableStructuralVariantAnalysis.of(annotations, fusions, disruptions);

            LOGGER.debug("persisting annotations to database");
            final StructuralVariantAnnotationDAO annotationDAO = new StructuralVariantAnnotationDAO(mDbAccess.context());

            annotationDAO.deleteAnnotationsForSample(sampleId);
            annotationDAO.write(analysis, sampleId);
        }
    }

    private List<EnrichedStructuralVariant> loadSVsFromVCF(final String sampleId)
    {
        List<EnrichedStructuralVariant> svList = Lists.newArrayList();

        try
        {
            LOGGER.debug("Loading indexed fasta reference file");
            final String fastaFileLocation = mCmdLineArgs.getOptionValue(REF_GENOME);
            final IndexedFastaSequenceFile indexedFastaSequenceFile = new IndexedFastaSequenceFile(new File(fastaFileLocation));

            LOGGER.debug("reading VCF File");
            final List<StructuralVariant> variants = readFromVcf(mCmdLineArgs.getOptionValue(VCF_FILE));

            LOGGER.debug("enriching structural variants based on purple data");
            svList = enrichStructuralVariants(sampleId,indexedFastaSequenceFile, variants);
        }
        catch(IOException e)
        {
            LOGGER.error("failed to read ref genome file");
            return svList;
        }

        if(mCmdLineArgs.hasOption(SV_PON_FILE))
        {
            svList = applyPONFilter(mCmdLineArgs.getOptionValue(SV_PON_FILE), svList);
        }

        return svList;
    }

    private List<EnrichedStructuralVariant> applyPONFilter(final String ponFilename, List<EnrichedStructuralVariant> svList)
    {
        SvPONAnnotator svPONAnnotator = new SvPONAnnotator();
        svPONAnnotator.loadPonFile(ponFilename);

        if (!svPONAnnotator.hasEntries())
            return svList;

        List<EnrichedStructuralVariant> updatedSVs = Lists.newArrayList();

        for (EnrichedStructuralVariant variant : svList)
        {
            boolean ponFiltered = false;
            if (variant.end() != null)
            {
                int ponCount = svPONAnnotator.getPonOccurenceCount(variant.chromosome(true),
                        variant.chromosome(false),
                        variant.position(true),
                        variant.position(false),
                        variant.orientation(true),
                        variant.orientation(false),
                        variant.type().toString());
                ponFiltered = ponCount > 1;
            }

            String filterString = variant.filter();
            assert filterString != null;

            Set<String> filterSet = Stream.of(filterString.split(";"))
                    .filter(s -> !s.equals("PASS"))
                    .filter(s -> !s.equals("."))
                    .filter(s -> !s.equals(""))
                    .collect(Collectors.toSet());

            if (ponFiltered)
            {
                filterSet.add(PON_FILTER_PON);
            }

            String filter = filterSet.size() == 0 ? "PASS" : filterSet.stream().sorted().collect(Collectors.joining(";"));

            final ImmutableEnrichedStructuralVariant updatedSV =
                    ImmutableEnrichedStructuralVariant.builder().from(variant).filter(filter).build();

            updatedSVs.add(updatedSV);
        }

        return updatedSVs;
    }

    @NotNull
    private static Set<String> tsgDriverGeneIDs()
    {
        Set<String> tsgDriverGeneIDs = Sets.newHashSet();
        Map<String, HmfTranscriptRegion> allGenes = HmfGenePanelSupplier.allGenesMap37();

        for (String gene : DndsDriverGeneLikelihoodSupplier.tsgLikelihood().keySet())
        {
            tsgDriverGeneIDs.add(allGenes.get(gene).geneID());
        }

        return tsgDriverGeneIDs;
    }

    @NotNull
    private static List<StructuralVariant> readFromVcf(@NotNull String vcfFileLocation) throws IOException
    {
        VariantContextFilter filter = variantContext -> {
            final Set<String> filters = Sets.newHashSet(variantContext.getFilters());
            filters.removeAll(ALLOWED_FILTERS);
            return variantContext.isNotFiltered() || filters.isEmpty();
        };

        final StructuralVariantFactory factory = new StructuralVariantFactory(filter);

        try (final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(vcfFileLocation,
                new VCFCodec(),
                false)) {
            reader.iterator().forEach(factory::addVariantContext);
        }

        return factory.results();
    }

    @NotNull
    private List<EnrichedStructuralVariant> enrichStructuralVariants(
            final String sampleId,
            @NotNull final IndexedFastaSequenceFile indexedFastaSequenceFile, @NotNull List<StructuralVariant> variants)
    {
        final PurityContext purityContext = mDbAccess.readPurityContext(sampleId);

        if (purityContext == null)
        {
            LOGGER.warn("sample({} unable to retrieve purple data, enrichment may be incomplete", sampleId);
        }

        final PurityAdjuster purityAdjuster = purityContext == null
                ? new PurityAdjuster(Gender.FEMALE, 1, 1)
                : new PurityAdjuster(purityContext.gender(), purityContext.bestFit().purity(), purityContext.bestFit().normFactor());

        final List<PurpleCopyNumber> copyNumberList = mDbAccess.readCopynumbers(sampleId);

        final Multimap<Chromosome, PurpleCopyNumber> copyNumbers =
                Multimaps.index(copyNumberList, x -> HumanChromosome.fromString(x.chromosome()));

        return new EnrichedStructuralVariantFactory(indexedFastaSequenceFile, purityAdjuster, copyNumbers).enrich(variants);
    }

    private final List<StructuralVariantAnnotation> createAnnotations(List<EnrichedStructuralVariant> enrichedVariants)
    {
        final Map<Integer, List<GeneAnnotation>> svIdGeneTranscriptsMap = mSvGeneTranscriptCollection.getSvIdGeneTranscriptsMap();

        List<StructuralVariantAnnotation> annotations = Lists.newArrayList();

        for(final EnrichedStructuralVariant var : enrichedVariants)
        {
            StructuralVariantAnnotation annotation = new StructuralVariantAnnotation(var);

            List<GeneAnnotation> genesList = svIdGeneTranscriptsMap.get(var.primaryKey());

            if(genesList != null)
            {
                for(GeneAnnotation geneAnnotation : genesList)
                {
                    geneAnnotation.setSvData(var);
                }

                annotation.annotations().addAll(genesList);
            }

            annotations.add(annotation);
        }

        return annotations;
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        if (cmd.hasOption(LOG_DEBUG))
            Configurator.setRootLevel(Level.DEBUG);

        StructuralVariantAnnotator svAnnotator = new StructuralVariantAnnotator(cmd);

        if(!svAnnotator.initialise())
            return;

        svAnnotator.run();

        LOGGER.info("run complete");
    }

    @NotNull
    private static Options createBasicOptions()
    {
        final Options options = new Options();
        options.addOption(VCF_FILE, true, "Path to the vcf file.");
        options.addOption(SAMPLE, true, "Tumor sample.");
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");
        options.addOption(ENSEMBL_DB, true, "Annotate structural variants using this Ensembl DB URI");
        options.addOption(ENSEMBL_DB_PASS, true, "Ensembl DB password if required");
        options.addOption(ENSEMBL_DB_USER, true, "Ensembl DB username if required");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");
        options.addOption(SV_PON_FILE, true, "PON file for SVs");
        options.addOption(REF_GENOME, true, "Path to the ref genome fasta file.");
        options.addOption(DATA_OUTPUT_DIR, true, "Path to persist annotations to file");

        SvFusionAnalyser.addCmdLineArgs(options);

        // testing options
        options.addOption(SOURCE_SVS_FROM_DB, false, "Skip annotations, including Ensemble DB data sync, for testing only)");
        options.addOption(LOAD_ANNOTATIONS_FROM_FILE, false, "Load existing annotations previously written to file");
        options.addOption(SKIP_DB_UPLOAD, false, "Skip uploading fusions and disruptions to database, off by default");
        options.addOption(TEST_SV_LIMIT, true, "Optional: only analyser X variants to save processing time");

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static DatabaseAccess databaseAccess(@NotNull final CommandLine cmd) throws SQLException
    {
        final String userName = cmd.getOptionValue(DB_USER);
        final String password = cmd.getOptionValue(DB_PASS);
        final String databaseUrl = cmd.getOptionValue(DB_URL);
        final String jdbcUrl = "jdbc:" + databaseUrl;
        return new DatabaseAccess(userName, password, jdbcUrl);
    }
}
