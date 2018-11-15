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
import java.util.List;
import java.util.Optional;

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
    private static final String VCF_HEADER_FILE = "vcf_header_file";
    private static final String REF_GENOME = "ref_genome";
    private static final String HIGH_CONFIDENCE_BED = "high_confidence_bed";
    private static final String WRITE_TO_DB = "write_to_db";
    private static final String WRITE_VCF_FILE = "write_vcf_file";
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

    private static final String VCF_FILE_SUFFIX = "_germline_variants.vcf";

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        if (cmd.hasOption(LOG_DEBUG)) {
            Configurator.setRootLevel(Level.DEBUG);
        }

        String sampleId = cmd.getOptionValue(SAMPLE);

        if (sampleId == null || sampleId.equals("*"))
        {
            sampleId = "";
        }

        String sampleDirectory = cmd.getOptionValue(SAMPLE_PATH);

        if(!sampleDirectory.endsWith("/"))
            sampleDirectory += "/";

        final String bachelorDirectory = sampleDirectory + BACHELOR_SUB_DIRECTORY;

        String bachelorInputFile;

        if(cmd.hasOption(BACH_INPUT_FILE))
            bachelorInputFile = cmd.getOptionValue(BACH_INPUT_FILE);
        else
            bachelorInputFile = bachelorDirectory + "/" + BACH_INPUT_FILE;

        BachelorDataCollection dataCollection = new BachelorDataCollection();
        dataCollection.setSampleId(sampleId);

        if(!dataCollection.loadBachelorData(bachelorInputFile))
            return;

        if (dataCollection.getBachelorVariants().isEmpty())
        {
            LOGGER.debug("sample({}) has no records to process", sampleId);
            return;
        }

        final List<BachelorGermlineVariant> bachRecords = dataCollection.getBachelorVariants();

        if(cmd.hasOption(WHITELIST_FILE) || cmd.hasOption(BLACKLIST_FILE))
        {
            rewriteFilteredBachelorRecords(bachRecords, cmd);
            LOGGER.info("filtering complete");
            return;
        }

        AlleleDepthLoader adLoader = new AlleleDepthLoader();
        adLoader.setSampleId(sampleId);

        if (!adLoader.loadMiniPileupData(bachelorDirectory) || !adLoader.applyPileupData(bachRecords))
            return;

        DatabaseAccess dbAccess;
        try
        {
            dbAccess = databaseAccess(cmd);
        }
        catch (SQLException e)
        {
            LOGGER.error("DB connection failed");
            return;
        }

        // create variant objects for VCF file writing and enrichment, and cache aginst bachelor record
        buildVariants(sampleId, bachRecords);

        annotateRecords(sampleId, bachRecords, cmd, dbAccess, sampleDirectory);

        int validRecordCount = 0;
        for (final BachelorGermlineVariant bachRecord : bachRecords)
        {
            if (bachRecord.isValid())
                ++validRecordCount;
        }

        if(validRecordCount == 0)
        {
            LOGGER.info("sample({})has no valid germline reports", sampleId, validRecordCount);
            return;
        }

        LOGGER.info("sample({}) writing {} germline reports to database", sampleId, validRecordCount);

        if (cmd.hasOption(WRITE_TO_DB))
        {
            writeToDatabase(sampleId, bachRecords, dbAccess);
        }

        if (cmd.hasOption(WRITE_VCF_FILE))
        {
            writeVcfFile(sampleId, bachRecords, bachelorDirectory, cmd.getOptionValue(VCF_HEADER_FILE));
        }

        writeToFile(sampleId, bachRecords, bachelorDirectory);

        LOGGER.info("run complete");
    }

    private static void rewriteFilteredBachelorRecords(final List<BachelorGermlineVariant> bachRecords, CommandLine cmdLineArgs)
    {
        String outputDir = cmdLineArgs.getOptionValue(SAMPLE_PATH);

        List<BachelorRecordFilter> whitelistFilters = loadBachelorFilters(cmdLineArgs.getOptionValue(WHITELIST_FILE));
        List<BachelorRecordFilter> blacklistFilters = loadBachelorFilters(cmdLineArgs.getOptionValue(BLACKLIST_FILE));

        if(whitelistFilters.isEmpty() && blacklistFilters.isEmpty())
            return;

        String outputFileName = outputDir;
        if (!outputFileName.endsWith("/"))
        {
            outputFileName += "/";
        }

        outputFileName += "bachelor_filtered_output.csv";

        Path outputFile = Paths.get(outputFileName);

        try
        {
            BufferedWriter writer = Files.newBufferedWriter(outputFile, StandardOpenOption.CREATE);

            writer.write("SAMPLEID,SOURCE,PROGRAM,ID,GENE,TRANSCRIPT_ID,CHROM,POS,REF,ALTS,EFFECTS,ANNOTATIONS,HGVS_PROTEIN,IS_HOMOZYGOUS,PHRED_SCORE,HGVS_CODING,MATCH_TYPE,CLNDN,CLNSIG");
            writer.newLine();

            for(int index = 0; index < bachRecords.size(); ++index)
            {
                final BachelorGermlineVariant bachRecord = bachRecords.get(index);

                if(index > 0 && (index % 1000) == 0)
                {
                    LOGGER.info("processed {} records", index);
                }

                boolean keepRecord = false;
                BachelorRecordFilter matchedFilter = null;

                if(FRAMESHIFT_VARIANT.isParentTypeOf(bachRecord.effects())
                || STOP_GAINED.isParentTypeOf(bachRecord.effects())
                || SPLICE_ACCEPTOR_VARIANT.isParentTypeOf(bachRecord.effects())
                || SPLICE_DONOR_VARIANT.isParentTypeOf(bachRecord.effects()))
                {
                    keepRecord = true;

                    for(final BachelorRecordFilter filter : blacklistFilters)
                    {
                        if(filter.matches(bachRecord))
                        {
                            keepRecord = false;
                            matchedFilter = filter;
                            break;
                        }
                    }
                }
                else if(MISSENSE_VARIANT.isParentTypeOf(bachRecord.effects())
                || SYNONYMOUS_VARIANT.isParentTypeOf(bachRecord.effects()))
                {
                    keepRecord = false;

                    for(final BachelorRecordFilter filter : whitelistFilters)
                    {
                        if(filter.matches(bachRecord))
                        {
                            keepRecord = true;
                            matchedFilter = filter;
                            break;
                        }
                    }
                }

                if(keepRecord || matchedFilter != null)
                {
                    writer.write(
                            String.format("%s,%s,%s,%s",
                                    bachRecord.patient(),
                                    bachRecord.source(),
                                    bachRecord.program(),
                                    bachRecord.variantId()));

                    writer.write(
                            String.format(",%s,%s,%s,%d,%s,%s",
                                    bachRecord.gene(),
                                    bachRecord.transcriptId(),
                                    bachRecord.chromosome(),
                                    bachRecord.position(),
                                    bachRecord.ref(),
                                    bachRecord.alts()));

                    writer.write(
                            String.format(",%s,%s,%s,%s,%d,%s,%s",
                                    bachRecord.effects(),
                                    bachRecord.annotations(),
                                    bachRecord.hgvsProtein(),
                                    bachRecord.isHomozygous(),
                                    bachRecord.phredScore(),
                                    bachRecord.hgvsCoding(),
                                    bachRecord.matchType()));

                    writer.write(String.format(",%s,%s",
                        matchedFilter != null ? matchedFilter.Diagnosis : "",
                        matchedFilter != null ? matchedFilter.Significance : ""));

                    writer.newLine();
                }
            }

            writer.close();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to outputFile");
        }
    }

    private static void buildVariants(final String sampleId, List<BachelorGermlineVariant> bachRecords)
    {
        for (final BachelorGermlineVariant bachRecord : bachRecords)
        {
            if (bachRecord.getRefCount() == 0 || bachRecord.getAltCount() == 0)
            {
                continue;
            }

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
                LOGGER.error("somatic variant creation failed");
            }
        }
    }

    private static void annotateRecords(final String sampleId, List<BachelorGermlineVariant> bachRecords, final CommandLine cmd,
            final DatabaseAccess dbAccess, final String sampleDirectory)
    {
        final String highConfidenceBed = cmd.getOptionValue(HIGH_CONFIDENCE_BED);
        final String fastaFileLocation = cmd.getOptionValue(REF_GENOME);

        List<SomaticVariant> variants = Lists.newArrayList();
        for (final BachelorGermlineVariant bachRecord : bachRecords)
        {
            if (bachRecord.getSomaticVariant() != null)
            {
                variants.add(bachRecord.getSomaticVariant());
            }
        }

        Multimap<String, GenomeRegion> highConfidenceRegions;
        IndexedFastaSequenceFile indexedFastaSequenceFile;

        try
        {
            LOGGER.debug("reading high confidence bed file");
            highConfidenceRegions = BEDFileLoader.fromBedFile(highConfidenceBed);

            LOGGER.debug("loading indexed fasta reference file");
            indexedFastaSequenceFile = new IndexedFastaSequenceFile(new File(fastaFileLocation));
        }
        catch (IOException e)
        {
            LOGGER.error("reference file loading failed");
            return;
        }

        final PurityContext purityContext;
        final Multimap<Chromosome, PurpleCopyNumber> copyNumbers;
        final Multimap<Chromosome, FittedRegion> copyNumberRegions;

        if(cmd.hasOption(PURPLE_DATA_DIRECTORY))
        {
            final String purplePath = sampleDirectory + cmd.getOptionValue(PURPLE_DATA_DIRECTORY);

            LOGGER.debug("loading purple data from file");

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
            LOGGER.debug("loading purple data from database");

            purityContext = dbAccess.readPurityContext(sampleId);

            if (purityContext == null)
            {
                LOGGER.warn("failed to read purity data");
            }

            copyNumbers = Multimaps.fromRegions(dbAccess.readCopynumbers(sampleId));

            copyNumberRegions = Multimaps.fromRegions(dbAccess.readCopyNumberRegions(sampleId));
        }


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

        LOGGER.debug("enriching variants");
        final EnrichedSomaticVariantFactory enrichedSomaticVariantFactory = new EnrichedSomaticVariantFactory(highConfidenceRegions,
                indexedFastaSequenceFile,
                new ClonalityFactory(purityAdjuster, clonalPloidy));

        final List<EnrichedSomaticVariant> enrichedVariants = enrichedSomaticVariantFactory.enrich(purityAdjustedVariants);

        for (BachelorGermlineVariant bachRecord : bachRecords)
        {
            boolean matched = false;

            for (final EnrichedSomaticVariant var : enrichedVariants)
            {
                if (bachRecord.chromosome().equals(var.chromosome()) && bachRecord.position() == var.position()) {
                    bachRecord.setEnrichedVariant(var);
                    matched = true;
                    break;
                }
            }

            if (!matched)
            {
                LOGGER.info("sample({}) enriched variant not found: var({}) gene({}) transcript({}) chr({}) position({})",
                        sampleId, bachRecord.variantId(), bachRecord.gene(), bachRecord.transcriptId(),
                        bachRecord.chromosome(), bachRecord.position());
            }
        }
    }

    private static void writeToDatabase(final String sampleId, final List<BachelorGermlineVariant> bachRecords, final DatabaseAccess dbAccess)
    {
        final GermlineVariantDAO germlineDAO = new GermlineVariantDAO(dbAccess.context());
        germlineDAO.write(sampleId, bachRecords);
    }

    public static void writeToFile(final String sampleId, final List<BachelorGermlineVariant> bachRecords, final String outputDir)
    {
        String outputFileName = outputDir;
        if (!outputFileName.endsWith("/"))
        {
            outputFileName += "/";
        }

        outputFileName += sampleId + "_germline_variants.csv";

        Path outputFile = Paths.get(outputFileName);

        try {
            BufferedWriter writer = Files.newBufferedWriter(outputFile);

            writer.write("SampleId,Program,Source,Chromosome,Position");

            writer.write(",Type,Ref,Alt,Gene,TranscriptId,DbsnpId,CosmicId,Effects,WorstCodingEffect,AltCount,RefCount");

            writer.write(",AdjCopyNumber,AdjustedVaf,HighConfidenceRegion,TrinucleotideContext,Microhomology,RepeatSequence,RepeatCount");

            writer.write(",HgvsProtein,HgvsCoding,Biallelic,Hotspot,Mappability,GermlineStatus,MinorAllelePloidy,Effects");

            writer.newLine();

            for (final BachelorGermlineVariant bachRecord : bachRecords)
            {
                if (!bachRecord.isValid())
                    continue;

                writer.write(
                        String.format("%s,%s,%s,%s,%d",
                        sampleId,
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
                        String.format(",%s,%s,%s,%s,%s,%s,%.2f,%s",
                                bachRecord.hgvsProtein(),
                                bachRecord.hgvsCoding(),
                                bachRecord.isBiallelic(),
                                region.hotspot(),
                                region.mappability(),
                                bachRecord.isHomozygous() ? "HOM" : "HET",
                                region.minorAllelePloidy(),
                                bachRecord.isLowScore() ? "ARTEFACT" : "PASS"));

                writer.newLine();
            }

            writer.close();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to outputFile");
        }
    }

    public static void rewriteFilteredRecordsToFile(final List<BachelorGermlineVariant> bachRecords, final String outputDir)
    {
        String outputFileName = outputDir;
        if (!outputFileName.endsWith("/"))
        {
            outputFileName += "/";
        }

        outputFileName += "bachelor_filtered_output.csv";

        Path outputFile = Paths.get(outputFileName);

        try
        {
            BufferedWriter writer = Files.newBufferedWriter(outputFile, StandardOpenOption.CREATE);

            writer.write("SAMPLEID,SOURCE,PROGRAM,ID,GENE,TRANSCRIPT_ID,CHROM,POS,REF,ALTS,EFFECTS,ANNOTATIONS,HGVS_PROTEIN,IS_HOMOZYGOUS,PHRED_SCORE,HGVS_CODING,MATCH_TYPE");
            writer.newLine();

            for (final BachelorGermlineVariant bachRecord : bachRecords)
            {
                writer.write(
                        String.format("%s,%s,%s,%s",
                                bachRecord.patient(),
                                bachRecord.source(),
                                bachRecord.program(),
                                bachRecord.variantId()));

                writer.write(
                        String.format(",%s,%s,%s,%d,%s,%s",
                                bachRecord.gene(),
                                bachRecord.transcriptId(),
                                bachRecord.chromosome(),
                                bachRecord.position(),
                                bachRecord.ref(),
                                bachRecord.alts()));

                writer.write(
                        String.format(",%s,%s,%s,%s,%d,%s,%s",
                                bachRecord.effects(),
                                bachRecord.annotations(),
                                bachRecord.hgvsProtein(),
                                bachRecord.isHomozygous(),
                                bachRecord.phredScore(),
                                bachRecord.hgvsCoding(),
                                bachRecord.matchType()));

                writer.newLine();
            }

            writer.close();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to outputFile");
        }
    }

    private static void writeVcfFile(final String sampleId, final List<BachelorGermlineVariant> bachRecords, final String outputDir, final String vcfHeaderFile)
    {
        LOGGER.debug("writing germline VCF");

        String outputFileName = outputDir;
        if (!outputFileName.endsWith("/"))
        {
            outputFileName += "/";
        }

        outputFileName += sampleId + VCF_FILE_SUFFIX;

        final File sampleVcfFile = new File(vcfHeaderFile);
        final VCFFileReader vcfReader = new VCFFileReader(sampleVcfFile, false);
        final VCFHeader sampleHeader = vcfReader.getFileHeader();

        final VCFHeader header = sampleHeader; // new VCFHeader();
        header.addMetaDataLine(new VCFInfoHeaderLine("AD", 1, VCFHeaderLineType.Integer, "Allele Depth for germline variant"));
        // header.addMetaDataLine(new VCFHeaderLine("BachelorPP", BachelorPP.class.getPackage().getImplementationVersion()));
        // header.addMetaDataLine(new VCFHeaderLine("CHROM", "POS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCPCT02040008T"));

        // Effects.UpdateVCFHeader(header);
        //  AlleleFrequency.UpdateVCFHeader(header);

        // setup VCF
        final VariantContextWriter writer = new VariantContextWriterBuilder().setReferenceDictionary(header.getSequenceDictionary())
                .setOutputFile(outputFileName)
                .build();

        writer.writeHeader(header);

        // variants.sort(new VariantContextComparator(header.getSequenceDictionary()));

        // write variants
        try
        {
            for (final BachelorGermlineVariant bachRecord : bachRecords)
            {
                writer.add(bachRecord.getVariantContext());
            }
        }
        catch (final NullPointerException e)
        {
            LOGGER.error("failed to write VCF record: {}", e.toString());
        }

        writer.close();
    }

    @NotNull
    private static Options createBasicOptions()
    {
        final Options options = new Options();
        options.addOption(SAMPLE, true, "Tumor sample.");

        options.addOption(SAMPLE_PATH, true, "Sample directory with a 'bachelor' sub-directory expected");
        options.addOption(BACH_INPUT_FILE, true, "Specific bachelor input file, if left out then assumes in sample path");
        options.addOption(VCF_HEADER_FILE, true, "Temp: VCF file header example");
        options.addOption(SAMPLE_PATH, true, "Name of input Mini Pileup file");
        options.addOption(REF_GENOME, true, "Path to the ref genome fasta file");
        options.addOption(HIGH_CONFIDENCE_BED, true, "Path to the high confidence bed file");
        options.addOption(PURPLE_DATA_DIRECTORY, true, "Sub-directory with sample path for purple data");

        options.addOption(DB_USER, true, "Database user name");
        options.addOption(DB_PASS, true, "Database password");
        options.addOption(DB_URL, true, "Database url");

        options.addOption(WRITE_VCF_FILE, false, "Whether to output a VCF file");
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
