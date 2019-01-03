package com.hartwig.hmftools.svannotation;

import static com.hartwig.hmftools.common.variant.structural.annotation.SvGeneTranscriptCollection.getSampleGeneAnnotationsFilename;
import static com.hartwig.hmftools.common.variant.structural.annotation.SvPONAnnotator.PON_FILTER_PASS;
import static com.hartwig.hmftools.common.variant.structural.annotation.SvPONAnnotator.PON_FILTER_PON;
import static com.hartwig.hmftools.svannotation.analysis.SvFusionAnalyser.FUSION_PAIRS_CSV;
import static com.hartwig.hmftools.svannotation.analysis.SvFusionAnalyser.PROMISCUOUS_FIVE_CSV;
import static com.hartwig.hmftools.svannotation.analysis.SvFusionAnalyser.PROMISCUOUS_THREE_CSV;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.sql.SQLException;
import java.util.List;
import java.util.Map;
import java.util.Set;

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
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFileLoader;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.StructuralVariantAnalysis;
import com.hartwig.hmftools.common.variant.structural.annotation.StructuralVariantAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.SvGeneTranscriptCollection;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.patientdb.dao.StructuralVariantAnnotationDAO;
import com.hartwig.hmftools.svannotation.analysis.StructuralVariantAnalyzer;
import com.hartwig.hmftools.svannotation.analysis.SvFusionAnalyser;

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
import htsjdk.variant.variantcontext.filter.VariantContextFilter;

public class StructuralVariantAnnotator
{
    private static final String SAMPLE = "sample";
    private static final String VCF_FILE = "vcf_file";
    private static final String ENSEMBL_DB = "ensembl_db";
    private static final String ENSEMBL_DB_USER = "ensembl_user";
    private static final String ENSEMBL_DB_PASS = "ensembl_pass";
    private static final String REF_GENOME = "ref_genome";
    private static final String ENSEMBL_DATA_DIR = "ensembl_data_dir";
    private static final String DATA_OUTPUT_DIR = "data_output_dir";

    private static final String SOURCE_SVS_FROM_DB = "source_svs_from_db";
    private static final String LOAD_ANNOTATIONS_FROM_FILE = "load_annotations";
    private static final String SKIP_DB_UPLOAD = "skip_db_upload";
    private static final String LOG_DEBUG = "log_debug";
    private static final String SAMPLE_RNA_FILE = "sample_rna_file";
    private static final String TRANS_EXON_FILE = "trans_exon_file";
    private static final String GENE_DATA_FILE = "gene_data_file";

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    private String mSampleId;
    private String mOutputDir; // for writing fusion data
    private String mEnsemblDataDir; // for loading and writing cached Ensembl data

    private final CommandLine mCmdLineArgs;
    private DatabaseAccess mDbAccess;
    private boolean mSourceSvFromDB;
    private SvGeneTranscriptCollection mSvGeneTranscriptCollection;
    private boolean mUploadAnnotations;

    // Let PON filtered SVs through since GRIDSS PON filtering is performed upstream
    private static final Set<String> ALLOWED_FILTERS = Sets.newHashSet("INFERRED", PON_FILTER_PON, PON_FILTER_PASS);

    private static final Logger LOGGER = LogManager.getLogger(StructuralVariantAnnotator.class);

    private StructuralVariantAnnotator(final CommandLine cmd)
    {
        mSampleId = "";
        mOutputDir = "";
        mEnsemblDataDir = "";
        mDbAccess = null;
        mSourceSvFromDB = false;
        mCmdLineArgs = cmd;
        mSvGeneTranscriptCollection = new SvGeneTranscriptCollection();
    }

    private boolean initialise()
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

        mOutputDir = mCmdLineArgs.getOptionValue(DATA_OUTPUT_DIR, "");
        mEnsemblDataDir = mCmdLineArgs.getOptionValue(ENSEMBL_DATA_DIR, "");
        mSvGeneTranscriptCollection.setDataPath(mEnsemblDataDir);
        mUploadAnnotations = !mCmdLineArgs.hasOption(SKIP_DB_UPLOAD);

        return true;
    }

    public boolean run()
    {
        // create the annotator - typically a MySQL connection to the Ensembl database
        MySQLAnnotator annotator = null;

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
            }
            catch (SQLException e)
            {
                LOGGER.warn("Ensembl DB connection failed: {}", e.toString());
                return false;
            }
        }

        LOGGER.debug("loading known fusion data");
        KnownFusionsModel knownFusionsModel;

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

        final StructuralVariantAnalyzer svAnalyser = new StructuralVariantAnalyzer(annotator, tsgDriverGeneIDs(), knownFusionsModel);

        if(mCmdLineArgs.hasOption(SAMPLE_RNA_FILE))
        {
            final String rnaFile = mCmdLineArgs.getOptionValue(SAMPLE_RNA_FILE);
            svAnalyser.getFusionAnalyser().loadSampleRnaData(rnaFile);
        }

        if (mCmdLineArgs.hasOption(TRANS_EXON_FILE) && mCmdLineArgs.hasOption(GENE_DATA_FILE))
        {
            final String transExonFile = mCmdLineArgs.getOptionValue(TRANS_EXON_FILE);
            mSvGeneTranscriptCollection.loadTranscriptExonData(transExonFile);

            final String geneDataFile = mCmdLineArgs.getOptionValue(GENE_DATA_FILE);
            mSvGeneTranscriptCollection.loadEnsemblGeneData(geneDataFile);
        }

        List<String> samplesList = Lists.newArrayList();

        if(mSampleId.isEmpty() || mSampleId.equals("*"))
        {
            samplesList = mDbAccess.structuralVariantSampleList("");

            LOGGER.info("loaded {} samples from database", samplesList.size());
        }
        else
        {
            if(mSampleId.contains(","))
            {
                String[] sampleIds = mSampleId.split(",");
                for (int i = 0; i < sampleIds.length; ++i)
                {
                    samplesList.add(sampleIds[i]);
                }
            }
            else
            {
                samplesList.add(mSampleId);
            }
        }

        for(final String sampleId : samplesList)
        {
            runSample(sampleId, svAnalyser, annotator != null ? annotator.getEnsemblDAO() : null);
        }

        svAnalyser.getFusionAnalyser().onCompleted();

        return true;
    }

    private boolean outputFileExists(final String sampleId)
    {
        final String outputFilename = getSampleGeneAnnotationsFilename(mEnsemblDataDir, sampleId);

        Path outputFile = Paths.get(outputFilename);

        return Files.exists(outputFile);
    }

    private void runSample(final String sampleId, final StructuralVariantAnalyzer svAnalyser, EnsemblDAO ensemblDAO)
    {
        //if(!mOverwriteEnsembleFiles && outputFileExists(sampleId))
        //    return;

        LOGGER.info("annotating variants for sample({})", sampleId);

        List<EnrichedStructuralVariant> enrichedVariants;

        if(mSourceSvFromDB)
        {
            // optionally load existing SVs from the database rather than from VCF
            enrichedVariants = mDbAccess.readStructuralVariants(sampleId);

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
            enrichedVariants = mDbAccess.readStructuralVariants(sampleId);
        }

        List<StructuralVariantAnnotation> annotations;

        if(!mEnsemblDataDir.isEmpty() && mSvGeneTranscriptCollection.loadSampleGeneTranscripts(sampleId))
        {
            annotations = createAnnotations(sampleId, enrichedVariants);

            LOGGER.debug("loaded {} Ensembl annotations from file", annotations.size());
        }
        else if(!mSvGeneTranscriptCollection.getGeneExonDataMap().isEmpty())
        {
            annotations = createAnnotations(enrichedVariants);

            LOGGER.debug("loaded {} Ensembl annotations from file", annotations.size());
        }
        else
        {
            LOGGER.debug("sample({}) finding Ensembl annotations", sampleId);

            annotations = svAnalyser.annotateVariants(enrichedVariants);

            LOGGER.debug("sample({}) matched {} annotations from Ensembl database", sampleId, annotations.size());

            // optionally persist to save having to look up this Ensembl data again
            if(!mEnsemblDataDir.isEmpty())
                mSvGeneTranscriptCollection.writeAnnotations(sampleId, annotations);
        }

        LOGGER.debug("sample({}) finding disruptions and fusions", sampleId);
        final StructuralVariantAnalysis analysis = svAnalyser.runOnAnnotations(annotations);

        if(!mOutputDir.isEmpty())
        {
            String clusterInfo = ",,";
            svAnalyser.getFusionAnalyser().writeFusions(analysis.fusions(), mOutputDir, sampleId, clusterInfo);
            svAnalyser.getFusionAnalyser().writeRnaMatchData(sampleId, mOutputDir, analysis.fusions(), annotations, ensemblDAO);
        }

        if(mUploadAnnotations)
        {
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
            LOGGER.debug("loading indexed fasta reference file");
            final String fastaFileLocation = mCmdLineArgs.getOptionValue(REF_GENOME);
            final IndexedFastaSequenceFile indexedFastaSequenceFile = new IndexedFastaSequenceFile(new File(fastaFileLocation));

            LOGGER.debug("reading VCF File");
            final List<StructuralVariant> variants = readFromVcf(mCmdLineArgs.getOptionValue(VCF_FILE));

            LOGGER.debug("enriching structural variants based on purple data");
            svList = enrichStructuralVariants(sampleId, indexedFastaSequenceFile, variants);
        }
        catch(IOException e)
        {
            LOGGER.error("failed to read ref genome file");
            return svList;
        }

        return svList;
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

        return StructuralVariantFileLoader.fromFile(vcfFileLocation, filter);
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

    private List<StructuralVariantAnnotation> createAnnotations(List<EnrichedStructuralVariant> enrichedVariants)
    {
        List<StructuralVariantAnnotation> annotations = Lists.newArrayList();

        for (final EnrichedStructuralVariant var : enrichedVariants)
        {
            StructuralVariantAnnotation annotation = new StructuralVariantAnnotation(var);

            List<GeneAnnotation> genesList = mSvGeneTranscriptCollection.findGeneAnnotationsBySv(
                    var.primaryKey(), true, var.chromosome(true), var.position(true), var.orientation(true));

            if(var.end() != null)
            {
                genesList.addAll(mSvGeneTranscriptCollection.findGeneAnnotationsBySv(
                        var.primaryKey(), false, var.chromosome(false), var.position(false), var.orientation(false)));
            }

            if(genesList == null)
                continue;

            for(GeneAnnotation geneAnnotation : genesList)
            {
                geneAnnotation.setSvData(var);
            }

            annotation.annotations().addAll(genesList);
            annotations.add(annotation);
        }

        return annotations;
    }

    private List<StructuralVariantAnnotation> createAnnotations(final String sampleId, List<EnrichedStructuralVariant> enrichedVariants)
    {
        final Map<Integer, List<GeneAnnotation>> svIdGeneTranscriptsMap = mSvGeneTranscriptCollection.getSvIdGeneTranscriptsMap();

        List<StructuralVariantAnnotation> annotations = Lists.newArrayList();
        boolean idsUpdated = false;

        for(final EnrichedStructuralVariant var : enrichedVariants)
        {
            StructuralVariantAnnotation annotation = new StructuralVariantAnnotation(var);

            List<GeneAnnotation> genesList = svIdGeneTranscriptsMap.get(var.primaryKey());

            if(genesList == null)
                continue;

            for(GeneAnnotation geneAnnotation : genesList)
            {
                geneAnnotation.setSvData(var);
            }

            annotation.annotations().addAll(genesList);
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
        options.addOption(SAMPLE_RNA_FILE, true, "Sample RNA data to match");
        options.addOption(TRANS_EXON_FILE, true, "Ensembl transcript exon data");
        options.addOption(GENE_DATA_FILE, true, "Ensembl gene data");
        options.addOption(REF_GENOME, true, "Path to the ref genome fasta file.");
        options.addOption(DATA_OUTPUT_DIR, true, "Path to persist annotations to file");
        options.addOption(ENSEMBL_DATA_DIR, true, "Cached Ensembl data path");

        SvFusionAnalyser.addCmdLineArgs(options);

        // testing options
        options.addOption(SOURCE_SVS_FROM_DB, false, "Skip annotations, including Ensemble DB data sync, for testing only)");
        options.addOption(LOAD_ANNOTATIONS_FROM_FILE, false, "Load existing annotations previously written to file");
        options.addOption(SKIP_DB_UPLOAD, false, "Skip uploading fusions and disruptions to database, off by default");

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
