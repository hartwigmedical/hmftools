package com.hartwig.hmftools.svannotation;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory.INFERRED;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory.PON_FILTER_PON;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.patientdb.LoadStructuralVariants.convertSvData;
import static com.hartwig.hmftools.svanalysis.gene.SvGeneTranscriptCollection.PRE_GENE_PROMOTOR_DISTANCE;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariantFactory;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFile;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFileLoader;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneDisruption;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;
import com.hartwig.hmftools.common.variant.structural.annotation.StructuralVariantAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.patientdb.dao.StructuralVariantFusionDAO;
import com.hartwig.hmftools.svanalysis.fusion.SvDisruptionAnalyser;
import com.hartwig.hmftools.svanalysis.fusion.SvFusionAnalyser;
import com.hartwig.hmftools.svanalysis.gene.EnsemblDAO;
import com.hartwig.hmftools.svanalysis.gene.SvGeneTranscriptCollection;

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
    /* SV Annotator can run in either a single-sample or a batch mode:
        - When in batch mode it sources SVs from the database and runs fusions logic (soon to be replaced by this function in Linx)
        - Otherwise it reads SVs from a VCF, enriches them with purple data and writes the results to file and/or DB
    */

    private String mSampleId;
    private String mOutputDir;
    private String mEnsemblDataDir;
    private final CommandLine mCmdLineArgs;
    private DatabaseAccess mDbAccess;
    private boolean mBatchMode; // for bulk fusion analysis

    private boolean mRunFusions;
    private SvGeneTranscriptCollection mSvGeneTranscriptCollection;
    private SvDisruptionAnalyser mDisruptionAnalyser;
    private SvFusionAnalyser mFusionAnalyser;
    private boolean mUploadAnnotations;

    private static final Logger LOGGER = LogManager.getLogger(StructuralVariantAnnotator.class);

    private StructuralVariantAnnotator(final CommandLine cmd)
    {
        mSampleId = "";
        mOutputDir = "";
        mEnsemblDataDir = "";
        mDbAccess = null;
        mCmdLineArgs = cmd;
        mBatchMode = false;
        mRunFusions = false;
        mDisruptionAnalyser = null;
        mFusionAnalyser = null;
        mSvGeneTranscriptCollection = null;
    }

    private boolean initialise()
    {
        if(!mCmdLineArgs.hasOption(SAMPLE) || mCmdLineArgs.getOptionValue(SAMPLE).isEmpty())
        {
            LOGGER.error("required Sample config missing");
            return false;
        }

        mSampleId = mCmdLineArgs.getOptionValue(SAMPLE);

        if(mCmdLineArgs.hasOption(DB_URL))
        {
            try
            {
                mDbAccess = databaseAccess(mCmdLineArgs);
            }
            catch (SQLException e)
            {
                LOGGER.error("Database connection failed: {}", e.toString());
                return false;
            }
        }

        mOutputDir = mCmdLineArgs.getOptionValue(DATA_OUTPUT_DIR, "");

        if (!mOutputDir.endsWith(File.separator))
            mOutputDir += File.separator;

        if(mCmdLineArgs.hasOption(ENSEMBL_DATA_DIR))
        {
            mRunFusions = true;
            mSvGeneTranscriptCollection = new SvGeneTranscriptCollection();
            mEnsemblDataDir = mCmdLineArgs.getOptionValue(ENSEMBL_DATA_DIR, "");

            mSvGeneTranscriptCollection.setDataPath(mEnsemblDataDir);
            mSvGeneTranscriptCollection.loadEnsemblData(false);

            mUploadAnnotations = !mCmdLineArgs.hasOption(SKIP_DB_UPLOAD);

            mFusionAnalyser = new SvFusionAnalyser(mCmdLineArgs, mSvGeneTranscriptCollection, mOutputDir);

            mDisruptionAnalyser = new SvDisruptionAnalyser(mCmdLineArgs, mSvGeneTranscriptCollection);
            mDisruptionAnalyser.setOutputDir(mOutputDir);
        }

        return true;
    }

    public boolean run()
    {
        List<String> samplesList = Lists.newArrayList();

        if (mSampleId.equals("*"))
        {
            if(mDbAccess == null)
            {
                LOGGER.warn("DB connection required for batch mode");
                return false;
            }

            mBatchMode = true;
            samplesList = mDbAccess.structuralVariantSampleList("");

            LOGGER.info("batch mode: loaded {} samples from database", samplesList.size());
        }
        else
        {
            samplesList.add(mSampleId);
        }

        for (String sampleId : samplesList)
        {
            runSample(sampleId, samplesList.size() > 1);
        }

        if(mRunFusions)
        {
            mFusionAnalyser.onCompleted();
            mDisruptionAnalyser.onCompleted();
            mSvGeneTranscriptCollection.close();
        }

        return true;
    }

    private void runSample(final String sampleId, boolean hasMultipleSamples)
    {
        LOGGER.info("annotating variants for sample({})", sampleId);

        List<StructuralVariantData> svDataList = null;

        if (mBatchMode)
        {
            // load existing SVs from the database rather than from VCF
            svDataList = mDbAccess.readStructuralVariantData(sampleId);

            if (svDataList.isEmpty())
            {
                LOGGER.debug("sample({}) no SVs loaded from DB", sampleId);
                return;
            }

            LOGGER.debug("Sample({}) loaded {} SVs from DB", sampleId, svDataList.size());
        }
        else
        {
            List<EnrichedStructuralVariant> enrichedVariants = loadSVsFromVCF();

            // generate a unique ID for each record
            int svId = 0;
            svDataList = Lists.newArrayList();

            for(EnrichedStructuralVariant var : enrichedVariants)
            {
                svDataList.add(convertSvData(var, svId++));
            }

            if(mDbAccess != null)
            {
                LOGGER.info("sample({}) persisting {} SVs to database", sampleId, enrichedVariants.size());
                mDbAccess.writeStructuralVariants(sampleId, svDataList);

                // Re-read the data to get primaryId field as a foreign key for disruptions and fusions
                // enrichedVariants = mDbAccess.readStructuralVariants(sampleId);
                svDataList = mDbAccess.readStructuralVariantData(sampleId);
            }
            else
            {
            }

            // write data to file
            try
            {
                final String svFilename = StructuralVariantFile.generateFilename(mOutputDir, sampleId);
                StructuralVariantFile.write(svFilename, svDataList);
            }
            catch(IOException e)
            {
                LOGGER.error("failed to write SV data: {}", e.toString());
            }
        }

        if(mRunFusions)
        {
            List<StructuralVariantAnnotation> annotations = Lists.newArrayList();

            if (mSvGeneTranscriptCollection.hasCachedEnsemblData())
            {
                annotations = addGeneAnnotations(svDataList);

                LOGGER.debug("Loaded {} Ensembl annotations from file", annotations.size());
            }

            List<GeneDisruption> disruptions = mDisruptionAnalyser.findDisruptions(annotations);
            List<GeneFusion> fusions = mFusionAnalyser.findFusions(annotations);

            LOGGER.debug("sample({}) found {} disruptions and {} fusions", sampleId, disruptions.size(), fusions.size());

            if (!mOutputDir.isEmpty())
            {
                mFusionAnalyser.writeFusions(fusions, sampleId, hasMultipleSamples);
                mDisruptionAnalyser.writeDisruptions(disruptions, sampleId, hasMultipleSamples);
            }

            if (mUploadAnnotations)
            {
                LOGGER.debug("persisting annotations to database");
                final StructuralVariantFusionDAO annotationDAO = new StructuralVariantFusionDAO(mDbAccess.context());

                annotationDAO.deleteAnnotationsForSample(sampleId);

                List<Transcript> allTranscripts = Lists.newArrayList();

                for (StructuralVariantAnnotation annotation : annotations)
                {
                    for (GeneAnnotation geneAnnotation : annotation.annotations())
                    {
                        allTranscripts.addAll(geneAnnotation.transcripts());
                    }
                }

                annotationDAO.writeBreakendsAndFusions(sampleId, allTranscripts, fusions);
            }
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

    private List<StructuralVariantAnnotation> addGeneAnnotations(List<StructuralVariantData> svDataList)
    {
        List<StructuralVariantAnnotation> annotations = Lists.newArrayList();

        for (final StructuralVariantData var : svDataList)
        {
            StructuralVariantAnnotation annotation = new StructuralVariantAnnotation(var);

            List<GeneAnnotation> genesList = mSvGeneTranscriptCollection.findGeneAnnotationsBySv(
                    var.id(), true, var.startChromosome(), var.startPosition(), var.startOrientation(), PRE_GENE_PROMOTOR_DISTANCE);

            if (var.type() != SGL)
            {
                genesList.addAll(mSvGeneTranscriptCollection.findGeneAnnotationsBySv(
                        var.id(), false, var.endChromosome(), var.endPosition(), var.endOrientation(), PRE_GENE_PROMOTOR_DISTANCE));
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

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        if (cmd.hasOption(LOG_DEBUG))
        {
            Configurator.setRootLevel(Level.DEBUG);
        }

        StructuralVariantAnnotator svAnnotator = new StructuralVariantAnnotator(cmd);

        if (!svAnnotator.initialise())
        {
            return;
        }

        svAnnotator.run();

        LOGGER.info("run complete");
    }

    private static final String SKIP_DB_UPLOAD = "skip_db_upload";

    // configuration
    private static final String SAMPLE = "sample";
    private static final String VCF_FILE = "vcf_file";
    private static final String REF_GENOME = "ref_genome";
    private static final String DATA_OUTPUT_DIR = "data_output_dir";
    public static final String ENSEMBL_DATA_DIR = "ensembl_data_dir";

    private static final String LOG_DEBUG = "log_debug";

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    @NotNull
    private static Options createBasicOptions()
    {
        final Options options = new Options();
        options.addOption(VCF_FILE, true, "Path to the vcf file");
        options.addOption(SAMPLE, true, "Tumor sample");
        options.addOption(DB_USER, true, "Database user name");
        options.addOption(DB_PASS, true, "Database password");
        options.addOption(DB_URL, true, "Database url");
        options.addOption(REF_GENOME, true, "Path to the ref genome fasta file.");
        options.addOption(DATA_OUTPUT_DIR, true, "Sample data directory");
        options.addOption(ENSEMBL_DATA_DIR, true, "Cached Ensembl data path");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");

        SvFusionAnalyser.addCmdLineArgs(options);
        EnsemblDAO.addCmdLineArgs(options);

        // testing options
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
