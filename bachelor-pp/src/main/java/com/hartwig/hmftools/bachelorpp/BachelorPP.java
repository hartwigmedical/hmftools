package com.hartwig.hmftools.bachelorpp;

import static com.hartwig.hmftools.bachelorpp.BachelorDataCollection.loadBachelorFilters;
import static com.hartwig.hmftools.common.variant.VariantConsequence.FRAMESHIFT_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.MISSENSE_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.SPLICE_ACCEPTOR_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.SPLICE_DONOR_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.STOP_GAINED;
import static com.hartwig.hmftools.common.variant.VariantConsequence.SYNONYMOUS_VARIANT;
import static com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationFactory.SNPEFF_IDENTIFIER;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.sql.SQLException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.bachelorpp.types.BachelorGermlineVariant;
import com.hartwig.hmftools.bachelorpp.types.BachelorRecordFilter;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.collect.Multimaps;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.FittedRegionFile;
import com.hartwig.hmftools.common.region.BEDFileLoader;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.ClonalityCutoffKernel;
import com.hartwig.hmftools.common.variant.ClonalityFactory;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariantFactory;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariant;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariantFactory;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
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
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class BachelorPP {

    private static final Logger LOGGER = LogManager.getLogger(BachelorPP.class);

    // config items
    private static final String SAMPLE = "sample";
    private static final String REF_GENOME = "ref_genome";
    private static final String HIGH_CONFIDENCE_BED = "high_confidence_bed";
    private static final String WRITE_TO_DB = "write_to_db";
    private static final String LOG_DEBUG = "log_debug";
    private static final String SAMPLE_PATH = "sample_path";
    private static final String PURPLE_DATA_DIRECTORY = "purple_data_dir";
    private static final String BACH_INPUT_FILE = "bachelor_file";

    // file locations
    private static final String BACHELOR_SUB_DIRECTORY = "bachelor";
    private static final String DEFAULT_BACH_INPUT_FILE = "bachelor_output.csv";
    private static final String WHITELIST_FILE = "whitelist_file";
    private static final String BLACKLIST_FILE = "blacklist_file";

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    private DatabaseAccess mDbAccess;
    private Multimap<String, GenomeRegion> mHighConfidenceRegions;
    private IndexedFastaSequenceFile mIndexedFastaSequenceFile;
    private boolean mIsBatchMode;
    private boolean mApplyFilters;
    private BufferedWriter mWriter;

    public BachelorPP()
    {
        mDbAccess = null;
        mHighConfidenceRegions = null;
        mIndexedFastaSequenceFile = null;
        mIsBatchMode = false;
        mApplyFilters = false;
        mWriter = null;
    }

    public boolean initialise(final CommandLine cmd)
    {
        try
        {
            mDbAccess = databaseAccess(cmd);
        }
        catch (SQLException e)
        {
            LOGGER.error("DB connection failed: {}", e.toString());
            return false;
        }

        if(cmd.hasOption(HIGH_CONFIDENCE_BED) && cmd.hasOption(REF_GENOME))
        {
            final String highConfidenceBed = cmd.getOptionValue(HIGH_CONFIDENCE_BED);
            final String fastaFileLocation = cmd.getOptionValue(REF_GENOME);

            try
            {
                LOGGER.debug("reading high confidence bed file");
                mHighConfidenceRegions = BEDFileLoader.fromBedFile(highConfidenceBed);

                LOGGER.debug("loading indexed fasta reference file");
                mIndexedFastaSequenceFile = new IndexedFastaSequenceFile(new File(fastaFileLocation));
            }
            catch (IOException e)
            {
                LOGGER.error("reference file loading failed");
                return false;
            }
        }

        return true;
    }

    public boolean run(final CommandLine cmd)
    {
        String sampleId = cmd.getOptionValue(SAMPLE);

        if (sampleId == null || sampleId.equals("*"))
        {
            sampleId = "";
            mIsBatchMode = true;
        }

        String sampleDirectory = cmd.getOptionValue(SAMPLE_PATH);

        if(!sampleDirectory.endsWith("/"))
            sampleDirectory += "/";

        String bachelorDirectory = sampleDirectory;

        if(!mIsBatchMode)
            bachelorDirectory += BACHELOR_SUB_DIRECTORY;

        String bachelorInputFile;

        if(cmd.hasOption(BACH_INPUT_FILE))
            bachelorInputFile = cmd.getOptionValue(BACH_INPUT_FILE);
        else
            bachelorInputFile = bachelorDirectory + "/" + DEFAULT_BACH_INPUT_FILE;

        BachelorDataCollection dataCollection = new BachelorDataCollection();
        dataCollection.setSampleId(sampleId);

        if(!dataCollection.loadBachelorData(bachelorInputFile))
            return false;

        if (dataCollection.getBachelorVariants().isEmpty())
        {
            LOGGER.debug("sample({}) has no records to process", sampleId);
            return false;
        }

        List<BachelorGermlineVariant> bachRecords = dataCollection.getBachelorVariants();

        mApplyFilters = cmd.hasOption(WHITELIST_FILE) || cmd.hasOption(BLACKLIST_FILE);

        if(mApplyFilters)
        {
            LOGGER.debug("applying white and black list filters");
            bachRecords = applyFilters(bachRecords, cmd);

            if (bachRecords.isEmpty())
            {
                LOGGER.info("all records were filtered out");
                return true;
            }

            LOGGER.debug("{} records after filtering", bachRecords.size());
        }

        if(!mIsBatchMode)
        {
            AlleleDepthLoader adLoader = new AlleleDepthLoader();
            adLoader.setSampleId(sampleId);

            if (!adLoader.loadMiniPileupData(bachelorDirectory) || !adLoader.applyPileupData(bachRecords))
                return false;
        }

        Map<String, List<BachelorGermlineVariant>> sampleRecordsMap = new HashMap();

        if(mIsBatchMode)
        {
            for (final BachelorGermlineVariant bachRecord : bachRecords)
            {
                List<BachelorGermlineVariant> sampleRecords = sampleRecordsMap.get(bachRecord.sampleId());

                if(sampleRecords == null)
                {
                    sampleRecords = Lists.newArrayList();
                    sampleRecordsMap.put(bachRecord.sampleId(), sampleRecords);
                }

                sampleRecords.add(bachRecord);
            }

            if(sampleRecordsMap.size() > 1)
            {
                LOGGER.debug("{} unique samples after filtering", sampleRecordsMap.size());
            }
        }
        else
        {
            sampleRecordsMap.put(sampleId, bachRecords);
        }

        for(final Map.Entry<String, List<BachelorGermlineVariant>> entry : sampleRecordsMap.entrySet())
        {
            final String specificSample = entry.getKey();
            List<BachelorGermlineVariant> sampleRecords = entry.getValue();

            // create variant objects for VCF file writing and enrichment, and cache aginst bachelor record
            buildVariants(specificSample, sampleRecords);

            annotateRecords(specificSample, sampleRecords, cmd, sampleDirectory);

            int validRecordCount = 0;
            for (final BachelorGermlineVariant bachRecord : sampleRecords)
            {
                if (bachRecord.isValid())
                    ++validRecordCount;
            }

            if(validRecordCount == 0)
            {
                LOGGER.info("sample({}) has no valid germline reports", specificSample, validRecordCount);
                continue;
            }

            if (cmd.hasOption(WRITE_TO_DB))
            {
                LOGGER.info("sample({}) writing {} germline reports to database", specificSample, validRecordCount);
                writeToDatabase(specificSample, sampleRecords);
            }

            writeToFile(specificSample, sampleRecords, bachelorDirectory);
        }

        try
        {
            if (mWriter != null)
                mWriter.close();
        }
        catch (IOException e)
        {
            LOGGER.error("error closing output file: {}", e.toString());
        }

        return true;
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        if (cmd.hasOption(LOG_DEBUG)) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        BachelorPP bachelorPP = new BachelorPP();

        if(!bachelorPP.initialise(cmd))
            return;

        if(!bachelorPP.run(cmd))
            return;

        LOGGER.info("run complete");
    }

    private void buildVariants(final String sampleId, List<BachelorGermlineVariant> bachRecords)
    {
        for (final BachelorGermlineVariant bachRecord : bachRecords)
        {
            VariantContextBuilder builder = new VariantContextBuilder();
            builder.id(bachRecord.variantId());
            builder.loc(bachRecord.chromosome(), bachRecord.position(), bachRecord.position() + bachRecord.ref().length() - 1);

            List<Allele> alleles = Lists.newArrayList();
            alleles.add(Allele.create(bachRecord.ref(), true));
            alleles.add(Allele.create(bachRecord.alts(), false));
            builder.alleles(alleles);

            List<Genotype> genoTypes = Lists.newArrayList();

            GenotypeBuilder gBuilder = new GenotypeBuilder(sampleId, builder.getAlleles());
            int[] adCounts = { bachRecord.getRefCount(), bachRecord.getAltCount() };
            gBuilder.AD(adCounts);
            genoTypes.add(gBuilder.make());

            builder.genotypes(genoTypes);
            VariantContext variantContext = builder.make();

            variantContext.getCommonInfo().putAttribute("AD", bachRecord.getAltCount(), true);
            variantContext.getCommonInfo().addFilter("PASS");
            variantContext.getCommonInfo().putAttribute(SNPEFF_IDENTIFIER, bachRecord.annotations());

            bachRecord.setVariantContext(variantContext);

            Optional<SomaticVariant> somVariant = SomaticVariantFactory.unfilteredInstance().createVariant(sampleId, variantContext);

            if (somVariant.isPresent())
            {
                bachRecord.setSomaticVariant(somVariant.get());
            }
            else
            {
                LOGGER.error("sample({}) var({}:{}) somatic variant creation failed", sampleId, bachRecord.variantId(), bachRecord.position());
            }
        }
    }

    private void annotateRecords(final String sampleId, List<BachelorGermlineVariant> bachRecords, final CommandLine cmd,
            final String sampleDirectory)
    {
        final PurityContext purityContext;
        final Multimap<Chromosome, PurpleCopyNumber> copyNumbers;
        final Multimap<Chromosome, FittedRegion> copyNumberRegions;

        if(cmd.hasOption(PURPLE_DATA_DIRECTORY))
        {
            final String purplePath = sampleDirectory + cmd.getOptionValue(PURPLE_DATA_DIRECTORY);

            LOGGER.debug("sample({}) loading purple data from file", sampleId);

            try
            {
                purityContext = FittedPurityFile.read(purplePath, sampleId);

                List<PurpleCopyNumber> copyNumberData = PurpleCopyNumberFile.read(PurpleCopyNumberFile.generateFilename(purplePath, sampleId));

                copyNumbers = Multimaps.fromRegions(copyNumberData);

                List<FittedRegion> fittedRegionData = FittedRegionFile.read(FittedRegionFile.generateFilename(purplePath, sampleId));

                copyNumberRegions = Multimaps.fromRegions(fittedRegionData);
            }
            catch (IOException e)
            {
                LOGGER.error("failed to read purple data from {}: {}", purplePath, e.toString());
                return;
            }
        }
        else
        {
            LOGGER.debug("sample({}) loading purple data from database", sampleId);

            purityContext = mDbAccess.readPurityContext(sampleId);

            if (purityContext == null)
            {
                LOGGER.warn("failed to read purity data");
            }

            copyNumbers = Multimaps.fromRegions(mDbAccess.readCopynumbers(sampleId));

            copyNumberRegions = Multimaps.fromRegions(mDbAccess.readCopyNumberRegions(sampleId));
        }

        List<SomaticVariant> variants = bachRecords.stream()
                .filter(x -> x.getSomaticVariant() != null)
                .map(x -> x.getSomaticVariant())
                .collect(Collectors.toList());

        final PurityAdjuster purityAdjuster = purityContext == null
                ? new PurityAdjuster(Gender.FEMALE, 1, 1)
                : new PurityAdjuster(purityContext.gender(), purityContext.bestFit().purity(), purityContext.bestFit().normFactor());

        final PurityAdjustedSomaticVariantFactory purityAdjustmentFactory =
                new PurityAdjustedSomaticVariantFactory(purityAdjuster, copyNumbers, copyNumberRegions);

        final List<PurityAdjustedSomaticVariant> purityAdjustedVariants = purityAdjustmentFactory.create(variants);

        final double clonalPloidy = ClonalityCutoffKernel.clonalCutoff(purityAdjustedVariants);

        for(PurityAdjustedSomaticVariant var : purityAdjustedVariants)
        {
            for (BachelorGermlineVariant bachRecord : bachRecords)
            {
                if (bachRecord.chromosome().equals(var.chromosome()) && bachRecord.position() == var.position())
                {
                    double adjVaf = 0;

                    if(bachRecord.isHomozygous())
                        adjVaf = purityAdjuster.purityAdjustedVAFWithHomozygousNormal(var.chromosome(), var.adjustedCopyNumber(), var.alleleFrequency());
                    else
                        adjVaf = purityAdjuster.purityAdjustedVAFWithHeterozygousNormal(var.chromosome(), var.adjustedCopyNumber(), var.alleleFrequency());

                    if(Double.isNaN(adjVaf) || Double.isInfinite(adjVaf))
                        adjVaf = 0;

                    bachRecord.setAdjustedVaf(adjVaf);
                    break;
                }
            }
        }

        LOGGER.debug("sample({}) enriching variants", sampleId);

        final EnrichedSomaticVariantFactory enrichedSomaticVariantFactory = new EnrichedSomaticVariantFactory(
                mHighConfidenceRegions,
                mIndexedFastaSequenceFile,
                new ClonalityFactory(purityAdjuster, clonalPloidy));

        final List<EnrichedSomaticVariant> enrichedVariants = enrichedSomaticVariantFactory.enrich(purityAdjustedVariants);

        for (BachelorGermlineVariant bachRecord : bachRecords)
        {
            boolean matched = false;

            for (final EnrichedSomaticVariant var : enrichedVariants)
            {
                if (bachRecord.chromosome().equals(var.chromosome()) && bachRecord.position() == var.position())
                {
                    bachRecord.setEnrichedVariant(var);
                    matched = true;
                    break;
                }
            }

            if (!matched)
            {
                LOGGER.debug("sample({}) enriched variant not found: var({}) gene({}) transcript({}) chr({}) position({})",
                        sampleId, bachRecord.variantId(), bachRecord.gene(), bachRecord.transcriptId(),
                        bachRecord.chromosome(), bachRecord.position());
            }
        }
    }

    private static List<BachelorGermlineVariant> applyFilters(final List<BachelorGermlineVariant> bachRecords, CommandLine cmdLineArgs)
    {
        List<BachelorRecordFilter> whitelistFilters = loadBachelorFilters(cmdLineArgs.getOptionValue(WHITELIST_FILE));
        List<BachelorRecordFilter> blacklistFilters = loadBachelorFilters(cmdLineArgs.getOptionValue(BLACKLIST_FILE));

        if (whitelistFilters.isEmpty() && blacklistFilters.isEmpty())
            return bachRecords;

        List<BachelorGermlineVariant> filteredRecords = Lists.newArrayList();

        for (int index = 0; index < bachRecords.size(); ++index)
        {
            final BachelorGermlineVariant bachRecord = bachRecords.get(index);

            if (index > 0 && (index % 1000) == 0)
            {
                LOGGER.info("processed {} records", index);
            }

            boolean keepRecord = false;
            BachelorRecordFilter matchedFilter = null;

            if (bachRecord.hasEffect(FRAMESHIFT_VARIANT) || bachRecord.hasEffect(STOP_GAINED)
            || bachRecord.hasEffect(SPLICE_ACCEPTOR_VARIANT) || bachRecord.hasEffect(SPLICE_DONOR_VARIANT))
            {
                keepRecord = true;

                for (final BachelorRecordFilter filter : blacklistFilters)
                {
                    if (filter.matches(bachRecord))
                    {
                        matchedFilter = filter;
                        break;
                    }
                }
            }
            else if (bachRecord.hasEffect(MISSENSE_VARIANT) || bachRecord.hasEffect(SYNONYMOUS_VARIANT))
            {
                keepRecord = false;

                for (final BachelorRecordFilter filter : whitelistFilters)
                {
                    if (filter.matches(bachRecord))
                    {
                        keepRecord = true;
                        matchedFilter = filter;
                        break;
                    }
                }
            }

            if(matchedFilter != null)
            {
                bachRecord.setDiagnosis(matchedFilter.Diagnosis);
                bachRecord.setSignificance(matchedFilter.Significance);
            }

            if(keepRecord)
                filteredRecords.add(bachRecord);
        }

        return filteredRecords;
    }

    private void writeToDatabase(final String sampleId, final List<BachelorGermlineVariant> bachRecords)
    {
        final GermlineVariantDAO germlineDAO = new GermlineVariantDAO(mDbAccess.context());
        germlineDAO.write(sampleId, bachRecords);
    }

    public void writeToFile(final String sampleId, final List<BachelorGermlineVariant> bachRecords, final String outputDir)
    {
        String outputFileName = outputDir;
        if (!outputFileName.endsWith("/"))
        {
            outputFileName += "/";
        }

        if(sampleId.equals("*"))
            outputFileName += sampleId + "_germline_variants.csv";
        else
            outputFileName += "bachelor_germline_variants.csv";

        Path outputFile = Paths.get(outputFileName);

        try
        {
            BufferedWriter writer = null;

            if(mWriter == null)
            {
                mWriter = Files.newBufferedWriter(outputFile, StandardOpenOption.CREATE);

                writer = mWriter;

                writer.write("SampleId,Program,Source,Chromosome,Position");

                writer.write(",Type,Ref,Alt,Gene,TranscriptId,DbsnpId,CosmicId,Effects,WorstCodingEffect,AltCount,RefCount");

                writer.write(",AdjCopyNumber,AdjustedVaf,HighConfidenceRegion,TrinucleotideContext,Microhomology,RepeatSequence,RepeatCount");

                writer.write(",HgvsProtein,HgvsCoding,Biallelic,Hotspot,Mappability,GermlineStatus,MinorAllelePloidy,Effects,Filter,CodonInfo");

                if(mApplyFilters)
                {
                    writer.write(",ClinvarDiagnosis,ClinvarSignificance");
                }

                writer.newLine();
            }
            else
            {
                writer = mWriter;
            }

            for (final BachelorGermlineVariant bachRecord : bachRecords)
            {
                if (!bachRecord.isValid())
                    continue;

                writer.write(
                        String.format("%s,%s,%s,%s,%d",
                        bachRecord.sampleId(),
                        bachRecord.program(),
                        bachRecord.source(),
                        bachRecord.chromosome(),
                        bachRecord.position()));

                final EnrichedSomaticVariant region = bachRecord.getEnrichedVariant();

                writer.write(
                        String.format(",%s,%s,%s,%s,%s,%s,%s,%s,%s,%d,%d",
                        region.type(),
                        region.ref(),
                        region.alt(),
                        bachRecord.gene(),
                        bachRecord.transcriptId(),
                        region.dbsnpID() == null ? "" : region.dbsnpID(),
                        region.canonicalCosmicID() == null ? "" : region.canonicalCosmicID(),
                        bachRecord.effects(),
                        region.worstCodingEffect(),
                        bachRecord.getAltCount(),
                        bachRecord.getRefCount()));

                writer.write(
                        String.format(",%.2f,%.2f,%s,%s,%s,%s,%d",
                        region.adjustedCopyNumber(),
                        bachRecord.getAdjustedVaf(),
                        region.highConfidenceRegion(),
                        region.trinucleotideContext(),
                        region.microhomology(),
                        region.repeatSequence(),
                        region.repeatCount()));

                writer.write(
                        String.format(",%s,%s,%s,%s,%s,%s,%.2f,%s,%s",
                                bachRecord.hgvsProtein(),
                                bachRecord.hgvsCoding(),
                                bachRecord.isBiallelic(),
                                region.hotspot(),
                                region.mappability(),
                                bachRecord.isHomozygous() ? "HOM" : "HET",
                                region.minorAllelePloidy(),
                                bachRecord.isLowScore() ? "ARTEFACT" : "PASS",
                                bachRecord.codonInfo()));

                if(mApplyFilters)
                {
                    writer.write(String.format(",%s,%s", bachRecord.getDiagnosis(), bachRecord.getSignificance()));
                }

                writer.newLine();
            }
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to outputFile: {}", e.toString());
        }
    }

    @NotNull
    private static Options createBasicOptions()
    {
        final Options options = new Options();
        options.addOption(SAMPLE, true, "Tumor sample.");

        options.addOption(SAMPLE_PATH, true, "Sample directory with a 'bachelor' sub-directory expected");
        options.addOption(BACH_INPUT_FILE, true, "Specific bachelor input file, if left out then assumes in sample path");
        options.addOption(SAMPLE_PATH, true, "Name of input Mini Pileup file");
        options.addOption(REF_GENOME, true, "Path to the ref genome fasta file");
        options.addOption(HIGH_CONFIDENCE_BED, true, "Path to the high confidence bed file");
        options.addOption(PURPLE_DATA_DIRECTORY, true, "Sub-directory with sample path for purple data");

        options.addOption(DB_USER, true, "Database user name");
        options.addOption(DB_PASS, true, "Database password");
        options.addOption(DB_URL, true, "Database url");

        options.addOption(WRITE_TO_DB, false, "Whether to upload records to DB");
        options.addOption(WHITELIST_FILE, true, "Additional whitelist checks on bachelor input file");
        options.addOption(BLACKLIST_FILE, true, "Additional blacklist checks on bachelor input file");

        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, default level is Info");

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
