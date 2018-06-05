package com.hartwig.hmftools.bachelorpp;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.sql.SQLException;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multimaps;
import com.hartwig.hmftools.bachelorpp.types.BachelorGermlineVariant;
import com.hartwig.hmftools.common.gene.CanonicalTranscriptFactory;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.bed.BEDFileLoader;
import com.hartwig.hmftools.common.variant.ClonalityCutoffKernel;
import com.hartwig.hmftools.common.variant.ClonalityFactory;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariantFactory;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariant;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariantFactory;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

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
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class BachelorPP {

    private static final Logger LOGGER = LogManager.getLogger(BachelorPP.class);

    private static final String SAMPLE = "sample";
    private static final String WRITE_VCF_FILE = "write_vcf_file";
    private static final String BACH_INPUT_FILE = "bach_input_file";
    private static final String VCF_INPUT_FILE = "vcf_input_file";
    private static final String REF_GENOME = "ref_genome";
    private static final String HIGH_CONFIDENCE_BED = "high_confidence_bed";
    private static final String MPILEUP_DIR = "mpileup_dir";
    private static final String WRITE_TO_DB = "write_to_db";
    private static final String DATA_OUTPUT_PATH = "data_output_path";
    private static final String LOG_DEBUG = "log_debug";

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    private static final String VCF_FILE_SUFFIX = "_germline_variants.vcf";

    public static void main(@NotNull final String[] args) throws ParseException, SQLException {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        if (cmd.hasOption(LOG_DEBUG)) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        DatabaseAccess dbAccess = null;

        try {
            dbAccess = databaseAccess(cmd);
        }
        catch(SQLException e)
        {
            LOGGER.error("DB connection failed");
            return;
        }

        String sampleId = cmd.getOptionValue(SAMPLE);

        if(sampleId == null || sampleId.equals("*"))
            sampleId = "";

        AlleleDepthLoader adLoader = new AlleleDepthLoader();

        adLoader.setSampleId(sampleId);
        adLoader.loadBachelorMatchData(cmd.getOptionValue(BACH_INPUT_FILE));

        if(cmd.hasOption(MPILEUP_DIR))
        {
            adLoader.loadMiniPileupData(cmd.getOptionValue(MPILEUP_DIR));
        }

        final List<BachelorGermlineVariant> bachRecords = adLoader.getBachelorVariants();

        final String outputDir = cmd.getOptionValue(DATA_OUTPUT_PATH);

        // create variant objects for VCF file writing and enrichment, and cache aginst bachelor record
        buildVariants(sampleId, bachRecords);

        if(cmd.hasOption(WRITE_VCF_FILE)) {

            writeVcfFile(sampleId, bachRecords, outputDir);
        }

        annotateRecords(sampleId, bachRecords, cmd, dbAccess);

        if(cmd.hasOption(WRITE_TO_DB))
            writeToDatabase(sampleId, bachRecords, dbAccess);

        // re-enable once fields are known and if VCF file or DB GermlineVariant data isn't sufficient
        // writeToFile(sampleId, bachRecords, outputDir);

        LOGGER.info("run complete");
    }

    private static void buildVariants(final String sampleId, List<BachelorGermlineVariant> bachRecords)
    {
        for(final BachelorGermlineVariant bachRecord : bachRecords) {
            VariantContextBuilder builder = new VariantContextBuilder();
            builder.id(bachRecord.variantId());
            builder.loc(bachRecord.chromosome(), bachRecord.position(), bachRecord.position());

            List<String> alleleStrings = Lists.newArrayList();
            alleleStrings.add(bachRecord.ref());
            alleleStrings.add(bachRecord.alts());
            builder.alleles(alleleStrings);

            List<Genotype> genoTypes = Lists.newArrayList();

            GenotypeBuilder gBuilder = new GenotypeBuilder(sampleId, builder.getAlleles());
            int[] adCounts = { bachRecord.getRefCount(), bachRecord.getAltCount() };
            gBuilder.AD(adCounts);
            genoTypes.add(gBuilder.make());

            builder.genotypes(genoTypes);
            VariantContext variantContext = builder.make();

            variantContext.getCommonInfo().putAttribute("AD", bachRecord.getAltCount(), true);
            variantContext.getCommonInfo().addFilter("PASS");


            Optional<SomaticVariant> somVariant = SomaticVariantFactory.unfilteredInstance().createVariant(sampleId, variantContext);

            if(somVariant.isPresent())
                bachRecord.setSomaticVariant(somVariant.get());
        }
    }

    private static void writeVcfFile(final String sampleId, final List<BachelorGermlineVariant> bachRecords, final String outputDir)
    {
        String outputFileName = outputDir;
        if (!outputFileName.endsWith("/")) {
            outputFileName += "/";
        }

        outputFileName += sampleId + VCF_FILE_SUFFIX;

        final File sampleVcfFile = new File("/Users/charlesshale/data/runs/sample_somatic.vcf");
        final VCFFileReader vcfReader = new VCFFileReader(sampleVcfFile, false);
        final VCFHeader sampleHeader = vcfReader.getFileHeader();

        final VCFHeader header = sampleHeader; // new VCFHeader();
        header.addMetaDataLine(new VCFInfoHeaderLine("AD", 1, VCFHeaderLineType.Integer, "Allele Depth for germline variant"));
        header.addMetaDataLine(new VCFHeaderLine("BachelorPP", BachelorPP.class.getPackage().getImplementationVersion()));

        // Filter.UpdateVCFHeader(header);
        //  AlleleFrequency.UpdateVCFHeader(header);

        // setup VCF
        final VariantContextWriter writer = new VariantContextWriterBuilder().setReferenceDictionary(header.getSequenceDictionary())
                .setOutputFile(outputFileName)
                .build();

        writer.writeHeader(header);

        // variants.sort(new VariantContextComparator(header.getSequenceDictionary()));

        // write variants
        for(final BachelorGermlineVariant bachRecord : bachRecords)
        {
            writer.add(bachRecord.getVariantContext());
        }

        writer.close();
    }

    private static void annotateRecords(
            final String sampleId, List<BachelorGermlineVariant> bachRecords,
            final CommandLine cmd, final DatabaseAccess dbAccess)
    {
        final String highConfidenceBed = cmd.getOptionValue(HIGH_CONFIDENCE_BED);
        final String fastaFileLocation = cmd.getOptionValue(REF_GENOME);

        List<SomaticVariant> variants = bachRecords.stream().map(BachelorGermlineVariant::getSomaticVariant).collect(Collectors.toList());

        Multimap<String, GenomeRegion> highConfidenceRegions = null;
        IndexedFastaSequenceFile indexedFastaSequenceFile = null;

        try {

            // LOGGER.info("Reading somatic VCF File");
            // final List<SomaticVariant> variants = SomaticVariantFactory.filteredInstance(filter).fromVCFFile(sample, vcfFileLocation);
            // variants = SomaticVariantFactory.unfilteredInstance().fromVCFFile(sampleId, vcfFile);

             LOGGER.debug("reading high confidence bed file");
             highConfidenceRegions = BEDFileLoader.fromBedFile(highConfidenceBed);

             LOGGER.debug("loading indexed fasta reference file");
             indexedFastaSequenceFile = new IndexedFastaSequenceFile(new File(fastaFileLocation));
        }
        catch(IOException e)
        {
            LOGGER.error("Reference file loading failed");
            return;
        }

        LOGGER.debug("querying purple database");
        final PurityContext purityContext = dbAccess.readPurityContext(sampleId);

        if (purityContext == null) {
            LOGGER.warn("Unable to retrieve purple data. Enrichment may be incomplete.");
        }

        final PurityAdjuster purityAdjuster = purityContext == null
                ? new PurityAdjuster(Gender.FEMALE, 1, 1)
                : new PurityAdjuster(purityContext.gender(), purityContext.bestFit().purity(), purityContext.bestFit().normFactor());

        final Multimap<String, PurpleCopyNumber> copyNumbers =
                Multimaps.index(dbAccess.readCopynumbers(sampleId), PurpleCopyNumber::chromosome);

        final Multimap<String, FittedRegion> copyNumberRegions =
                Multimaps.index(dbAccess.readCopyNumberRegions(sampleId), FittedRegion::chromosome);

        LOGGER.debug("incorporating purple purity");
        final PurityAdjustedSomaticVariantFactory purityAdjustmentFactory =
                new PurityAdjustedSomaticVariantFactory(purityAdjuster, copyNumbers, copyNumberRegions);
        final List<PurityAdjustedSomaticVariant> purityAdjustedVariants = purityAdjustmentFactory.create(variants);

        final double clonalPloidy = ClonalityCutoffKernel.clonalCutoff(purityAdjustedVariants);

        LOGGER.debug("enriching variants");
        final EnrichedSomaticVariantFactory enrichedSomaticVariantFactory = new EnrichedSomaticVariantFactory(highConfidenceRegions,
                indexedFastaSequenceFile,
                new ClonalityFactory(purityAdjuster, clonalPloidy),
                CanonicalTranscriptFactory.create(HmfGenePanelSupplier.allGeneList()));

        final List<EnrichedSomaticVariant> enrichedVariants = enrichedSomaticVariantFactory.enrich(purityAdjustedVariants);

        for(final EnrichedSomaticVariant var : enrichedVariants)
        {
            for(BachelorGermlineVariant bachRecord : bachRecords)
            {
                if(bachRecord.chromosome().equals(var.chromosome()) && bachRecord.position() == var.position())
                {
                    bachRecord.setEnrichedVariant(var);
                    break;
                }
            }
        }
    }

    public static void writeToDatabase(final String sampleId, final List<BachelorGermlineVariant> bachRecords, final DatabaseAccess dbAccess)
    {
        LOGGER.info("writing germline reports to database");
        final GermlineVariantDAO germlineDAO = new GermlineVariantDAO(dbAccess.context());
        germlineDAO.write(sampleId, bachRecords);
    }

    public static void writeToFile(final String sampleId, final List<BachelorGermlineVariant> bachRecords, final String outputDir)
    {
        LOGGER.debug("writing germline reports to file");

        String outputFileName = outputDir;
        if (!outputFileName.endsWith("/")) {
            outputFileName += "/";
        }

        outputFileName += sampleId + "_germline_variants.csv";

        Path outputFile = Paths.get(outputFileName);

        try
        {
            BufferedWriter writer = Files.newBufferedWriter(outputFile);

            // writer.write("SampleId,SvId,Type,ChrStart,PosStart,OrientStart,ChrEnd,PosEnd,OrientEnd,Filters\n");

            for(final BachelorGermlineVariant bachRecord : bachRecords)
            {
//                writer.write(
//                        String.format("%s,%s,%s,%s,%d,%d,%s,%d,%d,%s",
//                                sampleId, bachRecord..id(), var.type(),
//                                var.chromosome(true), var.position(true), var.orientation(true),
//                                var.chromosome(false), var.position(false), var.orientation(false), filtersStr));
//
//                writer.newLine();
            }
        }
            catch (final IOException e) {
            LOGGER.error("error writing to outputFile");
        }
    }

    @NotNull
    private static Options createBasicOptions() {
        final Options options = new Options();
        options.addOption(SAMPLE, true, "Tumor sample.");

        options.addOption(BACH_INPUT_FILE, true, "Name of input Bachelor file");
        options.addOption(VCF_INPUT_FILE, true, "Optional input VCF file if not using generated one");
        options.addOption(MPILEUP_DIR, true, "Name of input Mini Pileup file");
        options.addOption(REF_GENOME, true, "Path to the ref genome fasta file");
        options.addOption(HIGH_CONFIDENCE_BED, true, "Path to the high confidence bed file");

        options.addOption(DB_USER, true, "Database user name");
        options.addOption(DB_PASS, true, "Database password");
        options.addOption(DB_URL, true, "Database url");

        options.addOption(DATA_OUTPUT_PATH, true, "CSV output directory");
        options.addOption(WRITE_VCF_FILE, false, "Whether to output a VCF file");
        options.addOption(WRITE_TO_DB, false, "Whether to upload records to DB");

        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, default level is Info");

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
