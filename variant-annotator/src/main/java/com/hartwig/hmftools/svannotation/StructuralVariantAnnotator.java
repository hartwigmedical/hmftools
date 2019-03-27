package com.hartwig.hmftools.svannotation;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory.INFERRED;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory.PON_FILTER_PON;
import static com.hartwig.hmftools.svannotation.SvGeneTranscriptCollection.PRE_GENE_PROMOTOR_DISTANCE;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariantFactory;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFileLoader;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneDisruption;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.ImmutableStructuralVariantAnalysis;
import com.hartwig.hmftools.common.variant.structural.annotation.StructuralVariantAnalysis;
import com.hartwig.hmftools.common.variant.structural.annotation.StructuralVariantAnnotation;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.patientdb.dao.StructuralVariantAnnotationDAO;
import com.hartwig.hmftools.svannotation.analysis.SvDisruptionAnalyser;
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
    private SvDisruptionAnalyser mDisruptionAnalyser;
    private SvFusionAnalyser mFusionAnalyser;

    private String mSampleId;
    private String mOutputDir;
    private String mEnsemblDataDir;

    private final CommandLine mCmdLineArgs;
    private DatabaseAccess mDbAccess;
    private boolean mSourceSvFromDB;
    private SvGeneTranscriptCollection mSvGeneTranscriptCollection;
    private boolean mUploadAnnotations;
    private boolean mWriteBreakends;

    private static final Logger LOGGER = LogManager.getLogger(StructuralVariantAnnotator.class);

    private StructuralVariantAnnotator(final CommandLine cmd)
    {
        mDisruptionAnalyser = null;
        mFusionAnalyser = null;

        mSampleId = "";
        mOutputDir = "";
        mEnsemblDataDir = "";
        mDbAccess = null;
        mSourceSvFromDB = false;
        mCmdLineArgs = cmd;
        mSvGeneTranscriptCollection = new SvGeneTranscriptCollection();
        mWriteBreakends = false;
    }

    private boolean initialise()
    {
        if(!mCmdLineArgs.hasOption(SAMPLE) || mCmdLineArgs.getOptionValue(SAMPLE).isEmpty())
        {
            LOGGER.error("required Sample config missing");
            return false;
        }

        mSampleId = mCmdLineArgs.getOptionValue(SAMPLE);

        try
        {
            mDbAccess = databaseAccess(mCmdLineArgs);
        }
        catch (SQLException e)
        {
            LOGGER.error("Database connection failed: {}", e.toString());
            return false;
        }

        mSourceSvFromDB = mCmdLineArgs.hasOption(SOURCE_SVS_FROM_DB);

        mOutputDir = mCmdLineArgs.getOptionValue(DATA_OUTPUT_DIR, "");

        if (!mOutputDir.endsWith(File.separator))
            mOutputDir += File.separator;

        mEnsemblDataDir = mCmdLineArgs.getOptionValue(ENSEMBL_DATA_DIR, "");
        mUploadAnnotations = !mCmdLineArgs.hasOption(SKIP_DB_UPLOAD);
        mWriteBreakends = mCmdLineArgs.hasOption(WRITE_BREAKENDS);

        mDisruptionAnalyser = new SvDisruptionAnalyser(mOutputDir);
        mFusionAnalyser = new SvFusionAnalyser(mCmdLineArgs, mSvGeneTranscriptCollection, mOutputDir);

        return true;
    }

    public boolean run()
    {
        List<String> samplesList = Lists.newArrayList();

        if (mEnsemblDataDir.isEmpty())
        {
            LOGGER.error("Ensembl data cache directory missing");
            return false;
        }

        mSvGeneTranscriptCollection.setDataPath(mEnsemblDataDir);
        mSvGeneTranscriptCollection.loadEnsemblData();

        if (samplesList.isEmpty())
        {
            if (mSampleId.equals("*"))
            {
                samplesList = mDbAccess.structuralVariantSampleList("");

                LOGGER.info("batch mode: loaded {} samples from database", samplesList.size());
            }
            else
            {
                samplesList.add(mSampleId);
            }
        }

        for (String sampleId : samplesList)
        {
            runSample(sampleId, samplesList.size() > 1);
        }

        mFusionAnalyser.onCompleted();
        mDisruptionAnalyser.onCompleted();
        mSvGeneTranscriptCollection.close();

        return true;
    }

    private void runSample(final String sampleId, boolean hasMultipleSamples)
    {
        LOGGER.info("annotating variants for sample({})", sampleId);

        List<EnrichedStructuralVariant> enrichedVariants;

        if (mSourceSvFromDB)
        {
            // Optionally load existing SVs from the database rather than from VCF
            enrichedVariants = mDbAccess.readStructuralVariants(sampleId);

            if (enrichedVariants.isEmpty())
            {
                LOGGER.debug("Sample({}) no SVs loaded from DB", sampleId);
                return;
            }

            LOGGER.debug("Sample({}) loaded {} SVs from DB", sampleId, enrichedVariants.size());
        }
        else
        {
            List<EnrichedStructuralVariant> svList = loadSVsFromVCF();

            LOGGER.info("Sample({}) persisting {} SVs to database", sampleId, svList.size());

            mDbAccess.writeStructuralVariants(sampleId, svList);

            // Re-read the data to get primaryId field as a foreign key for disruptions and fusions
            enrichedVariants = mDbAccess.readStructuralVariants(sampleId);
        }

        List<StructuralVariantAnnotation> annotations = Lists.newArrayList();

        if (mSvGeneTranscriptCollection.hasCachedEnsemblData())
        {
            annotations = createAnnotations(enrichedVariants);

            LOGGER.debug("Loaded {} Ensembl annotations from file", annotations.size());
        }

        List<GeneDisruption> disruptions = mDisruptionAnalyser.findDisruptions(annotations);
        List<GeneFusion> fusions = mFusionAnalyser.findFusions(annotations);

        LOGGER.debug("sample({}) found {} disruptions and {} fusions", sampleId, disruptions.size(), fusions.size());

        if (!mOutputDir.isEmpty())
        {
            String clusterInfo = ",,";
            mFusionAnalyser.writeFusions(fusions, sampleId, clusterInfo, hasMultipleSamples);
            mDisruptionAnalyser.writeDisruptions(disruptions, sampleId, hasMultipleSamples);

            if (mWriteBreakends)
            {
                mSvGeneTranscriptCollection.writeBreakendData(sampleId, annotations);
            }
        }

        if (mUploadAnnotations)
        {
            LOGGER.debug("persisting annotations to database");
            final StructuralVariantAnnotationDAO annotationDAO = new StructuralVariantAnnotationDAO(mDbAccess.context());

            annotationDAO.deleteAnnotationsForSample(sampleId);

            final StructuralVariantAnalysis analysis = ImmutableStructuralVariantAnalysis.of(annotations, fusions, disruptions);
            annotationDAO.write(analysis, sampleId);
        }
    }

    private List<EnrichedStructuralVariant> loadSVsFromVCF()
    {
        List<EnrichedStructuralVariant> svList = Lists.newArrayList();

        try
        {
            LOGGER.debug("loading indexed fasta reference file");
            final String fastaFileLocation = mCmdLineArgs.getOptionValue(REF_GENOME);
            final IndexedFastaSequenceFile indexedFastaSequenceFile = new IndexedFastaSequenceFile(new File(fastaFileLocation));

            final String vcfFile = mCmdLineArgs.getOptionValue(VCF_FILE);
            LOGGER.info("Reading VCF File: {}", vcfFile);
            final List<StructuralVariant> variants = readFromVcf(vcfFile);

            LOGGER.debug("enriching structural variants based on purple data");
            svList = enrichStructuralVariants(indexedFastaSequenceFile, variants);
        }
        catch (IOException e)
        {
            LOGGER.error("failed to read ref genome file");
            return svList;
        }

        return svList;
    }

    // Let PON filtered SVs through since GRIDSS PON filtering is performed upstream
    private static final Set<String> ALLOWED_FILTERS = Sets.newHashSet(INFERRED, PON_FILTER_PON, "PASS");

    @NotNull
    private static List<StructuralVariant> readFromVcf(@NotNull String vcfFileLocation) throws IOException
    {
        VariantContextFilter filter = variantContext ->
        {
            final Set<String> filters = Sets.newHashSet(variantContext.getFilters());
            filters.removeAll(ALLOWED_FILTERS);

            return variantContext.isNotFiltered() || filters.isEmpty();
        };

        return StructuralVariantFileLoader.fromFile(vcfFileLocation, filter);
    }

    @NotNull
    private List<EnrichedStructuralVariant> enrichStructuralVariants(@NotNull final IndexedFastaSequenceFile indexedFastaSequenceFile,
            @NotNull final List<StructuralVariant> variants)
    {
        return new EnrichedStructuralVariantFactory(indexedFastaSequenceFile).enrich(variants);
    }

    private List<StructuralVariantAnnotation> createAnnotations(List<EnrichedStructuralVariant> enrichedVariants)
    {
        List<StructuralVariantAnnotation> annotations = Lists.newArrayList();

        for (final EnrichedStructuralVariant var : enrichedVariants)
        {
            StructuralVariantAnnotation annotation = new StructuralVariantAnnotation(var);

            Integer primaryKey = var.primaryKey();
            assert primaryKey != null; // Not sure why this assert is valid here...

            List<GeneAnnotation> genesList = mSvGeneTranscriptCollection.findGeneAnnotationsBySv(primaryKey,
                    true,
                    var.chromosome(true),
                    var.position(true),
                    PRE_GENE_PROMOTOR_DISTANCE);

            if (var.end() != null)
            {
                genesList.addAll(mSvGeneTranscriptCollection.findGeneAnnotationsBySv(primaryKey,
                        false,
                        var.chromosome(false),
                        var.position(false),
                        PRE_GENE_PROMOTOR_DISTANCE));
            }

            if (genesList == null)
            {
                continue;
            }

            for (GeneAnnotation geneAnnotation : genesList)
            {
                geneAnnotation.setSvData(var);
            }

            annotation.annotations().addAll(genesList);
            annotations.add(annotation);
        }

        return annotations;
    }

    public static void writeEnsemblDataFiles(final CommandLine cmd)
    {
        EnsemblDAO ensemblData = new EnsemblDAO(cmd);
        ensemblData.writeDataCacheFiles(cmd.getOptionValue(DATA_OUTPUT_DIR));
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        if (cmd.hasOption(LOG_DEBUG))
        {
            Configurator.setRootLevel(Level.DEBUG);
        }

        if (cmd.hasOption(WRITE_ENSEMBL_CACHE))
        {
            writeEnsemblDataFiles(cmd);
            return;
        }

        StructuralVariantAnnotator svAnnotator = new StructuralVariantAnnotator(cmd);

        if (!svAnnotator.initialise())
        {
            return;
        }

        svAnnotator.run();

        LOGGER.info("run complete");
    }

    private static final String SOURCE_SVS_FROM_DB = "source_svs_from_db";
    private static final String LOAD_ANNOTATIONS_FROM_FILE = "load_annotations";
    private static final String WRITE_BREAKENDS = "write_breakends";

    private static final String SKIP_DB_UPLOAD = "skip_db_upload";

    // configuration
    private static final String SAMPLE = "sample";
    private static final String VCF_FILE = "vcf_file";
    private static final String REF_GENOME = "ref_genome";
    private static final String ENSEMBL_DATA_DIR = "ensembl_data_dir";
    private static final String DATA_OUTPUT_DIR = "data_output_dir";

    private static final String LOG_DEBUG = "log_debug";
    private static final String WRITE_ENSEMBL_CACHE = "write_ensembl_cache";

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    @NotNull
    private static Options createBasicOptions()
    {
        final Options options = new Options();
        options.addOption(VCF_FILE, true, "Path to the vcf file.");
        options.addOption(SAMPLE, true, "Tumor sample.");
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");
        options.addOption(REF_GENOME, true, "Path to the ref genome fasta file.");
        options.addOption(DATA_OUTPUT_DIR, true, "Path to persist annotations to file");
        options.addOption(ENSEMBL_DATA_DIR, true, "Cached Ensembl data path");
        options.addOption(WRITE_ENSEMBL_CACHE, false, "Write Ensembl cached data files and exit");

        SvFusionAnalyser.addCmdLineArgs(options);
        EnsemblDAO.addCmdLineArgs(options);

        // testing options
        options.addOption(SOURCE_SVS_FROM_DB, false, "Skip annotations, including Ensemble DB data sync, for testing only)");
        options.addOption(LOAD_ANNOTATIONS_FROM_FILE, false, "Load existing annotations previously written to file");
        options.addOption(WRITE_BREAKENDS, false, "Write breakend data to file");
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
