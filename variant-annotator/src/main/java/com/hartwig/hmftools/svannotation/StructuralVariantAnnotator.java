package com.hartwig.hmftools.svannotation;

import static com.hartwig.hmftools.svanalysis.annotators.SvPONAnnotator.PON_FILTER_PASS;
import static com.hartwig.hmftools.svanalysis.annotators.SvPONAnnotator.PON_FILTER_PON;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
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
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.svanalysis.annotators.SvPONAnnotator;
import com.hartwig.hmftools.svannotation.analysis.ImmutableStructuralVariantAnalysis;
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
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.VariantContextFilter;
import htsjdk.variant.vcf.VCFCodec;

public class StructuralVariantAnnotator
{

    // Let PON filtered SVs through since GRIDSS PON filtering is performed upstream
    private static Set<String> ALLOWED_FILTERS = Sets.newHashSet( "INFERRED", PON_FILTER_PON, PON_FILTER_PASS);
    private static final Logger LOGGER = LogManager.getLogger(StructuralVariantAnnotator.class);

    private static final String SAMPLE = "sample";
    private static final String VCF_FILE = "vcf_file";
    private static final String ENSEMBL_DB = "ensembl_db";
    private static final String ENSEMBL_DB_LOCAL = "local_ensembl";
    private static final String ENSEMBL_DB_USER = "ensembl_user";
    private static final String ENSEMBL_DB_PASS = "ensembl_pass";
    private static final String FUSION_PAIRS_CSV = "fusion_pairs_csv";
    private static final String PROMISCUOUS_FIVE_CSV = "promiscuous_five_csv";
    private static final String PROMISCUOUS_THREE_CSV = "promiscuous_three_csv";
    private static final String REF_GENOME = "ref_genome";
    private static final String DATA_OUTPUT_DIR = "data_output_dir";

    private static final String SOURCE_SVS_FROM_DB = "source_svs_from_db";
    private static final String LOAD_ANNOTATIONS_FROM_FILE = "load_annotations";
    private static final String LOG_DEBUG = "log_debug";
    private static final String SV_PON_FILE = "sv_pon_file";

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    private String mSampleId;

    private final CommandLine mCmdLineArgs;
    private DatabaseAccess mDbAccess;
    private boolean mSourceSvFromDB;

    public StructuralVariantAnnotator(final CommandLine cmd)
    {
        mSampleId = "";
        mDbAccess = null;
        mSourceSvFromDB = false;
        mCmdLineArgs = cmd;
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

        return true;
    }

    public boolean run()
    {
        LOGGER.info("annotating variants for sample({})", mSampleId);

        List<EnrichedStructuralVariant> enrichedVariants;

        if(mSourceSvFromDB)
        {
            // optionally load existing SVs from the database rather than from VCF
            enrichedVariants = mDbAccess.readStructuralVariants(mSampleId);

            if (enrichedVariants.isEmpty())
            {
                LOGGER.info("no SVs loaded from DB");
                return false;
            }

            LOGGER.debug("sample({}) loaded {} SVs from DB", mSampleId, enrichedVariants.size());
        }
        else
        {
            List<EnrichedStructuralVariant> svList = loadSVsFromFile();

            LOGGER.info("persisting {} SVs to database", svList.size());

            mDbAccess.writeStructuralVariants(mSampleId, svList);

            // re-read the data to get primaryId field as a foreign key for disruptions and fusions
            enrichedVariants = mDbAccess.readStructuralVariants(mSampleId);;
        }

        DatabaseAccess ensembleDBConn = null;
        VariantAnnotator annotator = null;
        try
        {
            if (mCmdLineArgs.hasOption(ENSEMBL_DB_LOCAL))
            {
                // CHSH: The same credential work for both hmf and ensembl DBs
                final String ensembleJdbcUrl = "jdbc:" + mCmdLineArgs.getOptionValue(ENSEMBL_DB);
                final String ensembleUser = mCmdLineArgs.getOptionValue(ENSEMBL_DB_USER);
                final String ensemblePassword = mCmdLineArgs.getOptionValue(ENSEMBL_DB_PASS);

                LOGGER.debug("connecting to local ensembl DB: {}", mCmdLineArgs.getOptionValue(ENSEMBL_DB));

                ensembleDBConn = new DatabaseAccess(ensembleUser, ensemblePassword, ensembleJdbcUrl);
            }

            annotator = ensembleDBConn != null
                    ? new MySQLAnnotator(ensembleDBConn.context())
                    : MySQLAnnotator.make("jdbc:" + mCmdLineArgs.getOptionValue(ENSEMBL_DB));
        }
        catch (SQLException e)
        {
            LOGGER.warn("Ensembl DB connection failed: {}", e.toString());
            return false;
        }

        LOGGER.debug("loading fusion data");
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

        // TODO (KODU): Add potentially actionable genes, see also instantation in patient reporter.
        final StructuralVariantAnalyzer analyzer = new StructuralVariantAnalyzer(annotator, tsgDriverGeneIDs(), knownFusionsModel);

        // final StructuralVariantAnalysis analysis = analyzer.run(enrichedVariants);

        List<StructuralVariantAnnotation> annotations;
        if(mCmdLineArgs.hasOption(LOAD_ANNOTATIONS_FROM_FILE) && mCmdLineArgs.hasOption(DATA_OUTPUT_DIR))
        {
            annotations = loadAnnotations(enrichedVariants);
        }
        else
        {
            annotations = analyzer.findAnnotations(enrichedVariants);

            // optionally persist to save having to look up this Ensembl data again
            if(mCmdLineArgs.hasOption(DATA_OUTPUT_DIR))
                writeAnnotations(annotations);
        }

        LOGGER.debug("finding disruptions and fusions");
        final List<GeneFusion> fusions = analyzer.findFusions(annotations);
        final List<GeneDisruption> disruptions = analyzer.findDisruptions(annotations);

        final StructuralVariantAnalysis analysis = ImmutableStructuralVariantAnalysis.of(annotations, fusions, disruptions);

        LOGGER.debug("persisting annotations to database");
        final StructuralVariantAnnotationDAO annotationDAO = new StructuralVariantAnnotationDAO(mDbAccess.context());
        annotationDAO.write(analysis);

        return true;
    }

    private void writeAnnotations(final List<StructuralVariantAnnotation> annotations)
    {
        String outputFilename = mCmdLineArgs.getOptionValue(DATA_OUTPUT_DIR);

        if(!outputFilename.endsWith("/"))
            outputFilename += "/";

        outputFilename += mSampleId + "_sv_ensembl_data.csv";

        try
        {
            Path outputFile = Paths.get(outputFilename);

            BufferedWriter writer = Files.newBufferedWriter(outputFile);

            for(final StructuralVariantAnnotation annotation : annotations)
            {
                for(final GeneAnnotation geneAnnotation : annotation.annotations())
                {
                    writer.write(annotation.variant().id());

                    // Gene info: isStart, geneName, geneStableId, geneStrand, synonyms, entrezIds, karyotypeBand
                    writer.write(
                            String.format(",%s,%s,%s,%d,%s,%s,%s",
                                    geneAnnotation.isStart(),
                                    geneAnnotation.geneName(),
                                    geneAnnotation.stableId(),
                                    geneAnnotation.strand(),
                                    "synonyms",
                                    "entrezIds",
                                    geneAnnotation.karyotypeBand()));



                    for(final Transcript transcript : geneAnnotation.transcripts())
                    {
                        // Transcript info: transcriptId,exonUpstream, exonUpstreamPhase, exonDownstream, exonDownstreamPhase, exonMax, canonical, codingStart, codingEnd
                        writer.write(
                                String.format(",%s,%d,%d,%d,%d,%d,%s,%d,%d",
                                        transcript.transcriptId(),
                                        transcript.exonUpstream(),
                                        transcript.exonUpstreamPhase(),
                                        transcript.exonDownstream(),
                                        transcript.exonDownstreamPhase(),
                                        transcript.exonMax(),
                                        transcript.isCanonical(),
                                        transcript.codingStart(),
                                        transcript.codingEnd()));

                        writer.newLine();
                    }
                }
            }
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing gene annotations");
        }
    }

    private final String getSampleGeneAnnotationsFilename()
    {
        String outputFilename = mCmdLineArgs.getOptionValue(DATA_OUTPUT_DIR);

        if(!outputFilename.endsWith("/"))
            outputFilename += "/";

        return outputFilename + mSampleId + "sv_gene_annotations.csv";
    }

    private static int VAR_ID_COL_INDEX = 0;
    private static int GENE_NAME_COL_INDEX = 2;
    private static int TRANSCRIPT_ID_COL_INDEX = 10;

    private final List<StructuralVariantAnnotation> loadAnnotations(List<EnrichedStructuralVariant> enrichedVariants)
    {
        List<StructuralVariantAnnotation> annotations = Lists.newArrayList();

        final String filename = getSampleGeneAnnotationsFilename();

        if (filename.isEmpty())
            return annotations;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            // skip field names
            String line = fileReader.readLine();

            if (line == null)
            {
                LOGGER.error("Empty copy number CSV file({})", filename);
                return annotations;
            }

            int varIndex = 0;
            EnrichedStructuralVariant currentVar = enrichedVariants.get(varIndex);
            StructuralVariantAnnotation annotation = new StructuralVariantAnnotation(currentVar);

            GeneAnnotation currentGene = null;
            Transcript currentTranscript = null;

            while ((line = fileReader.readLine()) != null)
            {
                // parse CSV data
                String[] items = line.split(",");

                int item = VAR_ID_COL_INDEX;

                // check if still on the same variant
                final String varId = items[item++];

                if(!varId.equals(currentVar.id()))
                {
                    // done with current variant
                    annotations.add(annotation);

                    ++varIndex;

                    if(varIndex >= enrichedVariants.size())
                        break;

                    currentVar = enrichedVariants.get(varIndex);
                    annotation = new StructuralVariantAnnotation(currentVar);
                }

                // isStart, geneName, geneStableId, geneStrand, synonyms, entrezIds, karyotypeBand

                final List<String> synonyms = Lists.newArrayList();
                final List<Integer> entrezIds = Lists.newArrayList();

                final String geneName = items[GENE_NAME_COL_INDEX];

                if(currentGene != null && !currentGene.geneName().equals(geneName))
                {
                    // add to annotation and prepare a new one
                    annotation.annotations().add(currentGene);
                }

                if(currentGene == null || !currentGene.geneName().equals(geneName))
                {
                    currentGene = new GeneAnnotation(
                            currentVar,
                            Boolean.parseBoolean(items[item++]),
                            items[item++],
                            items[item++],
                            Integer.parseInt(items[item++]),
                            synonyms,
                            entrezIds,
                            items[item++]);
                }

                final String transcriptId = items[TRANSCRIPT_ID_COL_INDEX];

                if(currentTranscript != null && !currentTranscript.transcriptId().equals(transcriptId))
                {
                    currentGene.transcripts().add(currentTranscript);
                }

                if(currentTranscript == null || !currentTranscript.transcriptId().equals(transcriptId))
                {
                    // transcriptId, exonUpstream, exonUpstreamPhase, exonDownstream, exonDownstreamPhase, exonMax, canonical, codingStart, codingEnd
                    currentTranscript = new Transcript(
                            currentGene,
                            items[item++],
                            Integer.parseInt(items[item++]),
                            Integer.parseInt(items[item++]),
                            Integer.parseInt(items[item++]),
                            Integer.parseInt(items[item++]),
                            Integer.parseInt(items[item++]),
                            Boolean.parseBoolean(items[item++]),
                            Long.parseLong(items[item++]),
                            Long.parseLong(items[item++]));
                }
            }
        }
        catch(IOException e)
        {
            LOGGER.error("failed to load sample gene annotations({}): {}", filename, e.toString());
        }

        return annotations;
    }

    private List<EnrichedStructuralVariant> loadSVsFromFile()
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
            svList = enrichStructuralVariants(indexedFastaSequenceFile, variants);
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
        Map<String, HmfTranscriptRegion> allGenes = HmfGenePanelSupplier.allGenesMap();

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
            @NotNull final IndexedFastaSequenceFile indexedFastaSequenceFile, @NotNull List<StructuralVariant> variants)
    {
        final PurityContext purityContext = mDbAccess.readPurityContext(mSampleId);

        if (purityContext == null)
        {
            LOGGER.warn("sample({} unable to retrieve purple data, enrichment may be incomplete", mSampleId);
        }

        final PurityAdjuster purityAdjuster = purityContext == null
                ? new PurityAdjuster(Gender.FEMALE, 1, 1)
                : new PurityAdjuster(purityContext.gender(), purityContext.bestFit().purity(), purityContext.bestFit().normFactor());

        final List<PurpleCopyNumber> copyNumberList = mDbAccess.readCopynumbers(mSampleId);

        final Multimap<Chromosome, PurpleCopyNumber> copyNumbers =
                Multimaps.index(copyNumberList, x -> HumanChromosome.fromString(x.chromosome()));

        return new EnrichedStructuralVariantFactory(indexedFastaSequenceFile, purityAdjuster, copyNumbers).enrich(variants);
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
        options.addOption(FUSION_PAIRS_CSV, true, "Path towards a CSV containing white-listed gene fusion pairs.");
        options.addOption(PROMISCUOUS_FIVE_CSV, true, "Path towards a CSV containing white-listed promiscuous 5' genes.");
        options.addOption(PROMISCUOUS_THREE_CSV, true, "Path towards a CSV containing white-listed promiscuous 3' genes.");
        options.addOption(ENSEMBL_DB, true, "Annotate structural variants using this Ensembl DB URI");
        options.addOption(ENSEMBL_DB_PASS, true, "Ensembl DB password if required");
        options.addOption(ENSEMBL_DB_USER, true, "Ensembl DB username if required");
        options.addOption(ENSEMBL_DB_LOCAL, false, "Connect to local Ensembl DB");
        options.addOption(SOURCE_SVS_FROM_DB, false, "Skip annotations, including Ensemble DB data sync, for testing only)");
        options.addOption(LOAD_ANNOTATIONS_FROM_FILE, false, "Load existing annotations previously written to file");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");
        options.addOption(SV_PON_FILE, true, "PON file for SVs");
        options.addOption(REF_GENOME, true, "Path to the ref genome fasta file.");
        options.addOption(DATA_OUTPUT_DIR, true, "Path to persist annotations to file");
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
