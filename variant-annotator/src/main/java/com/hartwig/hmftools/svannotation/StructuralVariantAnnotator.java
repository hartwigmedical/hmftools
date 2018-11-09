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
import java.nio.file.StandardOpenOption;
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
    private static final String ENSEMBL_DB_USER = "ensembl_user";
    private static final String ENSEMBL_DB_PASS = "ensembl_pass";
    private static final String FUSION_PAIRS_CSV = "fusion_pairs_csv";
    private static final String PROMISCUOUS_FIVE_CSV = "promiscuous_five_csv";
    private static final String PROMISCUOUS_THREE_CSV = "promiscuous_three_csv";
    private static final String REF_GENOME = "ref_genome";
    private static final String DATA_OUTPUT_DIR = "data_output_dir";

    private static final String SOURCE_SVS_FROM_DB = "source_svs_from_db";
    private static final String LOAD_ANNOTATIONS_FROM_FILE = "load_annotations";
    private static final String SKIP_DB_UPLOAD = "skip_db_upload";
    private static final String WRITE_RESULTS_FILE = "write_results_file";
    private static final String LOG_DEBUG = "log_debug";
    private static final String SV_PON_FILE = "sv_pon_file";

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    private static final String TEST_SV_LIMIT = "test_limit_sv_count";

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
                    ensembleDBConn = new DatabaseAccess(ensembleUser, ensemblePassword, ensembleUrl);
                    annotator = new MySQLAnnotator(ensembleDBConn.context());
                }
                else
                {
                    MySQLAnnotator.make(ensembleUrl);
                }
            } catch (SQLException e)
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
        final StructuralVariantAnalyzer analyzer = new StructuralVariantAnalyzer(annotator, tsgDriverGeneIDs(), knownFusionsModel);

        // final StructuralVariantAnalysis analysis = analyzer.run(enrichedVariants);

        List<StructuralVariantAnnotation> annotations;
        if(mCmdLineArgs.hasOption(LOAD_ANNOTATIONS_FROM_FILE) && mCmdLineArgs.hasOption(DATA_OUTPUT_DIR))
        {
            annotations = loadAnnotations(enrichedVariants);

            LOGGER.debug("loaded {} Ensembl annotations from file", annotations.size());
        }
        else if(annotator != null)
        {
            int testSvLimit = Integer.parseInt(mCmdLineArgs.getOptionValue(TEST_SV_LIMIT, "0"));

            if(testSvLimit > 0)
            {
                List<EnrichedStructuralVariant> subset = Lists.newArrayList(enrichedVariants.subList(0, testSvLimit));
                annotations = analyzer.findAnnotations(subset);
            }
            else
            {
                annotations = analyzer.findAnnotations(enrichedVariants);
            }

            LOGGER.debug("matched {} annotations from Ensembl database", annotations.size());

            // optionally persist to save having to look up this Ensembl data again
            if(mCmdLineArgs.hasOption(DATA_OUTPUT_DIR))
                writeAnnotations(annotations);
        }
        else
        {
            LOGGER.error("Ensemble data not loaded from DB nor file");
            return false;
        }

        LOGGER.debug("finding disruptions and fusions");
        final List<GeneFusion> fusions = analyzer.findFusions(annotations);
        final List<GeneDisruption> disruptions = analyzer.findDisruptions(annotations);

        LOGGER.debug("found {} disruptions and {} fusions", disruptions.size(), fusions.size());

        final StructuralVariantAnalysis analysis = ImmutableStructuralVariantAnalysis.of(annotations, fusions, disruptions);

        if(mCmdLineArgs.hasOption(WRITE_RESULTS_FILE))
        {
            writeFusions(analysis.fusions());
            // writeDisruptions(analysis.disruptions());
        }

        if(!mCmdLineArgs.hasOption(SKIP_DB_UPLOAD))
        {
            LOGGER.debug("persisting annotations to database");
            final StructuralVariantAnnotationDAO annotationDAO = new StructuralVariantAnnotationDAO(mDbAccess.context());

            annotationDAO.deleteAnnotationsForSample(mSampleId);
            annotationDAO.write(analysis, mSampleId);

        }

        return true;
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

    private void writeAnnotations(final List<StructuralVariantAnnotation> annotations)
    {
        LOGGER.debug("writing {} annotations to file", annotations.size());

        String outputFilename = getSampleGeneAnnotationsFilename();

        try
        {
            Path outputFile = Paths.get(outputFilename);

            BufferedWriter writer = Files.newBufferedWriter(outputFile, StandardOpenOption.CREATE);

            for(final StructuralVariantAnnotation annotation : annotations)
            {
                if(annotation.annotations().isEmpty())
                {
                    // LOGGER.debug("SV({}) has no annotations", annotation.variant().primaryKey());
                    continue;
                }

                for(final GeneAnnotation geneAnnotation : annotation.annotations())
                {
                    String synonymnsStr = "";
                    for(final String syn : geneAnnotation.synonyms())
                    {
                        if(!synonymnsStr.isEmpty())
                            synonymnsStr += ";";

                        synonymnsStr += syn;
                    }

                    String entrezIdsStr = "";
                    for(final Integer eId : geneAnnotation.entrezIds())
                    {
                        if(!entrezIdsStr.isEmpty())
                            entrezIdsStr += ";";

                        entrezIdsStr += eId;
                    }

                    for(final Transcript transcript : geneAnnotation.transcripts())
                    {
                        writer.write(String.format("%d", annotation.variant().primaryKey()));

                        // Gene info: isStart, geneName, geneStableId, geneStrand, synonyms, entrezIds, karyotypeBand
                        writer.write(
                                String.format(",%s,%s,%s,%d,%s,%s,%s",
                                        geneAnnotation.isStart(),
                                        geneAnnotation.geneName(),
                                        geneAnnotation.stableId(),
                                        geneAnnotation.strand(),
                                        synonymnsStr,
                                        entrezIdsStr,
                                        geneAnnotation.karyotypeBand()));

                        // Transcript info: transcriptId,exonUpstream, exonUpstreamPhase, exonDownstream, exonDownstreamPhase, exonStart, exonEnd, exonMax, canonical, codingStart, codingEnd
                        writer.write(
                                String.format(",%s,%d,%d,%d,%d,%d,%d,%d,%s,%d,%d",
                                        transcript.transcriptId(),
                                        transcript.exonUpstream(),
                                        transcript.exonUpstreamPhase(),
                                        transcript.exonDownstream(),
                                        transcript.exonDownstreamPhase(),
                                        transcript.exonStart(),
                                        transcript.exonEnd(),
                                        transcript.exonMax(),
                                        transcript.isCanonical(),
                                        transcript.codingStart(),
                                        transcript.codingEnd()));

                        writer.newLine();
                    }
                }
            }

            writer.close();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing gene annotations");
        }
    }

    private void writeFusions(final List<GeneFusion> fusions)
    {
        if(fusions.isEmpty())
            return;

        String outputFilename = mCmdLineArgs.getOptionValue(DATA_OUTPUT_DIR);

        if(!outputFilename.endsWith("/"))
            outputFilename += "/";

        outputFilename += mSampleId + "_" + "sv_fusions.csv";

        try
        {
            Path outputFile = Paths.get(outputFilename);

            BufferedWriter writer = Files.newBufferedWriter(outputFile, StandardOpenOption.CREATE);

            writer.write("SampleId");
            writer.write(",StartSvId,StartChr,StartPos,StartOrient,StartType");
            writer.write(",StartGene,StartTranscript,StartRegionType,StartExon,StartPhase,StartExonStart,StartExonEnd,StartCodingStart,StartCodingEnd");
            writer.write(",EndSvId,EndChr,EndPos,EndOrient,EndType");
            writer.write(",EndGene,EndTranscript,EndRegionType,EndExon,EndPhase,EndExonStart,EndExonEnd,EndCodingStart,EndCodingEnd");
            writer.newLine();

            for(final GeneFusion fusion : fusions)
            {
                final Transcript startTrans = fusion.upstreamLinkedAnnotation();
                final Transcript endTrans = fusion.downstreamLinkedAnnotation();

                boolean startUsesStart = startTrans.parent().isStart();
                boolean endUsesStart = endTrans.parent().isStart();

                final EnrichedStructuralVariant startVar = startTrans.parent().variant();
                final EnrichedStructuralVariant endVar = endTrans.parent().variant();

                writer.write(String.format("%s", mSampleId));

                // write upstream SV, transcript and exon info
                writer.write(
                        String.format(",%d,%s,%d,%d,%s",
                                startVar.primaryKey(), startVar.chromosome(startUsesStart), startVar.position(startUsesStart), startVar.orientation(startUsesStart), startVar.type()));

                writer.write(
                        String.format(",%s,%s,%s,%d,%d,%d,%d,%d,%d",
                                startTrans.parent().geneName(), startTrans.transcriptId(), startTrans.getRegionType(), startTrans.exonUpstream(), startTrans.exonUpstreamPhase(),
                                startTrans.exonStart(), startTrans.exonEnd(),
                                startTrans.codingStart() != null ? startTrans.codingStart() : 0,
                                startTrans.codingEnd() != null ? startTrans.codingEnd() : 0));

                writer.write(
                        String.format(",%d,%s,%d,%d,%s",
                                endVar.primaryKey(), endVar.chromosome(endUsesStart), endVar.position(endUsesStart), endVar.orientation(endUsesStart), endVar.type()));

                writer.write(
                        String.format(",%s,%s,%s,%d,%d,%d,%d,%d,%d",
                                endTrans.parent().geneName(), endTrans.transcriptId(), endTrans.getRegionType(), endTrans.exonUpstream(), endTrans.exonUpstreamPhase(),
                                endTrans.exonStart(), endTrans.exonEnd(),
                                endTrans.codingStart() != null ? endTrans.codingStart() : 0,
                                endTrans.codingEnd() != null ? endTrans.codingEnd() : 0));

                writer.newLine();
            }

            writer.close();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing gene annotations");
        }
    }

    private void writeDisruptions(final List<GeneDisruption> disruptions)
    {
        if(disruptions.isEmpty())
            return;

        String outputFilename = mCmdLineArgs.getOptionValue(DATA_OUTPUT_DIR);

        if(!outputFilename.endsWith("/"))
            outputFilename += "/";

        outputFilename += mSampleId + "_" + "sv_disruption.csv";

        try
        {
            Path outputFile = Paths.get(outputFilename);

            BufferedWriter writer = Files.newBufferedWriter(outputFile, StandardOpenOption.CREATE);

            writer.close();
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

        return outputFilename + mSampleId + "_" + "sv_ensembl_data.csv";
    }

    private static int VAR_ID_COL_INDEX = 0;

    // gene data: isStart, geneName, geneStableId, geneStrand, synonyms, entrezIds, karyotypeBand
    private static int GENE_IS_START_COL_INDEX = 1;
    private static int GENE_NAME_COL_INDEX = 2;
    private static int GENE_STABLE_ID_COL_INDEX = 3;
    private static int GENE_STRAND_INDEX = 4;
    private static int GENE_SYNS_COL_INDEX = 5;
    private static int GENE_EIDS_COL_INDEX = 6;
    private static int GENE_KARYOTYPE_COL_INDEX = 7;

    // transcript data: transcriptId, exonUpstream, exonUpstreamPhase, exonDownstream, exonDownstreamPhase, exonStart, exonEnd, exonMax, canonical, codingStart, codingEnd
    private static int TRANSCRIPT_ID_COL_INDEX = 8;
    private static int TRANSCRIPT_EUS_COL_INDEX = 9;
    private static int TRANSCRIPT_EUP_COL_INDEX = 10;
    private static int TRANSCRIPT_EDS_COL_INDEX = 11;
    private static int TRANSCRIPT_EDP_COL_INDEX = 12;
    private static int TRANSCRIPT_ESTART_COL_INDEX = 13;
    private static int TRANSCRIPT_EEND_COL_INDEX = 14;
    private static int TRANSCRIPT_EMAX_COL_INDEX = 15;
    private static int TRANSCRIPT_CAN_COL_INDEX = 16;
    private static int TRANSCRIPT_CS_COL_INDEX = 17;
    private static int TRANSCRIPT_CE_COL_INDEX = 18;

    private final List<StructuralVariantAnnotation> loadAnnotations(List<EnrichedStructuralVariant> enrichedVariants)
    {
        List<StructuralVariantAnnotation> annotations = Lists.newArrayList();

        final String filename = getSampleGeneAnnotationsFilename();

        if (filename.isEmpty())
            return annotations;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine();

            if (line == null)
            {
                LOGGER.error("empty copy number CSV file({})", filename);
                return annotations;
            }

            int varIndex = 0;
            EnrichedStructuralVariant currentVar = enrichedVariants.get(varIndex);
            StructuralVariantAnnotation annotation = new StructuralVariantAnnotation(currentVar);

            GeneAnnotation currentGene = null;

            int fileIndex = 0;

            while (line != null)
            {
                // parse CSV data
                String[] items = line.split(",");

                // check if still on the same variant
                final String varId = items[VAR_ID_COL_INDEX];

                if(!varId.equals(currentVar.primaryKey().toString()))
                {
                    if(annotation != null)
                    {
                        // done with current variant and gene
                        if(currentGene != null)
                            annotation.annotations().add(currentGene);

                        annotations.add(annotation);

                        currentGene = null;
                    }

                    ++varIndex;

                    if(varIndex >= enrichedVariants.size())
                    {
                        LOGGER.warn("variants exhausted before file end, at index({}): fileVar({})", fileIndex, varId);
                        break;
                    }

                    currentVar = enrichedVariants.get(varIndex);
                    annotation = null;

                    if(!varId.equals(currentVar.primaryKey().toString()))
                    {
                        // it's possible that the SV has no ensembl data, in which it needs to be skipped over
                        // LOGGER.debug("variant mismatch at index({}): currentVar({}) vs fileVar({}), skipping SV", fileIndex, currentVar.primaryKey(), varId);
                        continue;
                    }

                    annotation = new StructuralVariantAnnotation(currentVar);
                }

                // isStart, geneName, geneStableId, geneStrand, synonyms, entrezIds, karyotypeBand
                final String geneName = items[GENE_NAME_COL_INDEX];

                if(currentGene == null || !currentGene.geneName().equals(geneName))
                {
                    if(currentGene != null)
                    {
                        // add to annotation and prepare a new one
                        annotation.annotations().add(currentGene);
                    }

                    String[] synonymsStr = items[GENE_SYNS_COL_INDEX].split(";");
                    final List<String> synonyms = Lists.newArrayList(synonymsStr);

                    String[] entrezIdStr = items[GENE_EIDS_COL_INDEX].split(";");

                    final List<Integer> entrezIds = Lists.newArrayList();

                    for (int i = 0; i < entrezIdStr.length; ++i)
                    {
                        if(!entrezIdStr[i].isEmpty())
                            entrezIds.add(Integer.parseInt(entrezIdStr[i]));
                    }

                    currentGene = new GeneAnnotation(
                            currentVar,
                            Boolean.parseBoolean(items[GENE_IS_START_COL_INDEX]),
                            geneName,
                            items[GENE_STABLE_ID_COL_INDEX],
                            Integer.parseInt(items[GENE_STRAND_INDEX]),
                            synonyms,
                            entrezIds,
                            items[GENE_KARYOTYPE_COL_INDEX]);
                }

                final String transcriptId = items[TRANSCRIPT_ID_COL_INDEX];

                // transcriptId, exonUpstream, exonUpstreamPhase, exonDownstream, exonDownstreamPhase, exonMax, canonical, codingStart, codingEnd
                Transcript transcript = new Transcript(
                        currentGene,
                        transcriptId,
                        Integer.parseInt(items[TRANSCRIPT_EUS_COL_INDEX]),
                        Integer.parseInt(items[TRANSCRIPT_EUP_COL_INDEX]),
                        Integer.parseInt(items[TRANSCRIPT_EDS_COL_INDEX]),
                        Integer.parseInt(items[TRANSCRIPT_EDP_COL_INDEX]),
                        Long.parseLong(items[TRANSCRIPT_ESTART_COL_INDEX]),
                        Long.parseLong(items[TRANSCRIPT_EEND_COL_INDEX]),
                        Integer.parseInt(items[TRANSCRIPT_EMAX_COL_INDEX]),
                        Boolean.parseBoolean(items[TRANSCRIPT_CAN_COL_INDEX]),
                        items[TRANSCRIPT_CS_COL_INDEX].equals("null") ? null : Long.parseLong(items[TRANSCRIPT_CS_COL_INDEX]),
                        items[TRANSCRIPT_CE_COL_INDEX].equals("null") ? null : Long.parseLong(items[TRANSCRIPT_CE_COL_INDEX]));

                currentGene.addTranscript(transcript);

                line = fileReader.readLine();

                if(line == null)
                {
                    // add the last annotation
                    annotations.add(annotation);
                    break;
                }

                ++fileIndex;
            }

        }
        catch(IOException e)
        {
            LOGGER.error("failed to load sample gene annotations({}): {}", filename, e.toString());
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
        options.addOption(FUSION_PAIRS_CSV, true, "Path towards a CSV containing white-listed gene fusion pairs.");
        options.addOption(PROMISCUOUS_FIVE_CSV, true, "Path towards a CSV containing white-listed promiscuous 5' genes.");
        options.addOption(PROMISCUOUS_THREE_CSV, true, "Path towards a CSV containing white-listed promiscuous 3' genes.");
        options.addOption(ENSEMBL_DB, true, "Annotate structural variants using this Ensembl DB URI");
        options.addOption(ENSEMBL_DB_PASS, true, "Ensembl DB password if required");
        options.addOption(ENSEMBL_DB_USER, true, "Ensembl DB username if required");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");
        options.addOption(SV_PON_FILE, true, "PON file for SVs");
        options.addOption(REF_GENOME, true, "Path to the ref genome fasta file.");
        options.addOption(DATA_OUTPUT_DIR, true, "Path to persist annotations to file");

        // testing options
        options.addOption(SOURCE_SVS_FROM_DB, false, "Skip annotations, including Ensemble DB data sync, for testing only)");
        options.addOption(LOAD_ANNOTATIONS_FROM_FILE, false, "Load existing annotations previously written to file");
        options.addOption(SKIP_DB_UPLOAD, false, "Skip uploading fusions and disruptions to database, off by default");
        options.addOption(WRITE_RESULTS_FILE, false, "Write fusions and disruptions to file, off by default");
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
