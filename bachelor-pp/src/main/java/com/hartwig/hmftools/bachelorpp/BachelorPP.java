package com.hartwig.hmftools.bachelorpp;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationFactory.SNPEFF_IDENTIFIER;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.bachelorpp.types.BachelorGermlineVariant;
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

public class BachelorPP
{

    private static final Logger LOGGER = LogManager.getLogger(BachelorPP.class);

    // config items
    private static final String SAMPLE = "sample";
    private static final String REF_GENOME = "ref_genome";
    private static final String HIGH_CONFIDENCE_BED = "high_confidence_bed";
    private static final String WRITE_TO_DB = "write_to_db";
    private static final String LOG_DEBUG = "log_debug";
    private static final String SAMPLE_PATH = "sample_path"; // path to the run directory for the sample
    private static final String PURPLE_DATA_DIRECTORY = "purple_data_dir"; // purple data directory within the sample fir
    private static final String BACH_DIRECTORY = "bachelor_dir"; // usually defaults to the 'bachelor' subdirectory of the sample dir
    private static final String BACH_INPUT_FILE = "bachelor_file"; // full path
    private static final String MAX_READ_COUNT = "max_read_count";

    // file locations
    private static final String SAMPLE_LIST_FILE = "sample_list_file";

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    private static final String DEFAULT_BACH_DIRECTORY = "bachelor";
    private static final String DEFAULT_BACH_INPUT_FILE = "bachelor_output.csv";

    private DatabaseAccess mDbAccess;
    private Multimap<String, GenomeRegion> mHighConfidenceRegions;
    private IndexedFastaSequenceFile mIndexedFastaSequenceFile;
    private boolean mIsBatchMode;
    private boolean mUploadRecordsToDB;
    private BufferedWriter mWriter;
    private AlleleDepthLoader mAllelDepthLoader;
    private List<String> mLimitedSampleList;

    private BachelorPP()
    {
        mDbAccess = null;
        mHighConfidenceRegions = null;
        mIndexedFastaSequenceFile = null;
        mIsBatchMode = false;
        mUploadRecordsToDB = false;
        mAllelDepthLoader = null;
        mWriter = null;

        mLimitedSampleList = Lists.newArrayList();
    }

    private boolean initialise(final CommandLine cmd)
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

        if (cmd.hasOption(HIGH_CONFIDENCE_BED) && cmd.hasOption(REF_GENOME))
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

        mUploadRecordsToDB = cmd.hasOption(WRITE_TO_DB);

        return true;
    }

    private boolean run(final CommandLine cmd)
    {
        String sampleId = cmd.getOptionValue(SAMPLE);

        if (sampleId == null || sampleId.equals("*"))
        {
            LOGGER.info("running in batch mode");
            sampleId = "";
            mIsBatchMode = true;

            if(cmd.hasOption(SAMPLE_LIST_FILE))
                loadSampleListFile(cmd.getOptionValue(SAMPLE_LIST_FILE));
        }

        String sampleDirectory = cmd.getOptionValue(SAMPLE_PATH);

        if (!sampleDirectory.endsWith(File.separator))
        {
            sampleDirectory += File.separator;
        }

        String bachelorDataDir;
        if (cmd.hasOption(BACH_DIRECTORY))
        {
            bachelorDataDir = sampleDirectory + cmd.getOptionValue(BACH_DIRECTORY);
            if (!bachelorDataDir.endsWith(File.separator))
            {
                bachelorDataDir += File.separator;
            }

            LOGGER.debug("using configured bachelor directory: {}", bachelorDataDir);
        }
        else
        {
            bachelorDataDir = sampleDirectory + DEFAULT_BACH_DIRECTORY + File.separator;
            LOGGER.debug("using default bachelor data dir: {}", bachelorDataDir);
        }

        String bachelorInputFile;
        if (cmd.hasOption(BACH_INPUT_FILE))
        {
            bachelorInputFile = cmd.getOptionValue(BACH_INPUT_FILE);
            LOGGER.info("loading specific input file: {}", bachelorInputFile);
        }
        else
        {
            bachelorInputFile = bachelorDataDir + DEFAULT_BACH_INPUT_FILE;
            LOGGER.info("loading sample default file: {}", bachelorInputFile);
        }

        if (!mIsBatchMode)
        {
            mAllelDepthLoader = new AlleleDepthLoader();
            mAllelDepthLoader.setSampleId(sampleId);

            if (!mAllelDepthLoader.loadMiniPileupData(bachelorDataDir))
            {
                return false;
            }
        }

        BachelorDataCollection dataCollection = new BachelorDataCollection();
        dataCollection.setSampleId(sampleId, mLimitedSampleList);
        dataCollection.setMaxReadCount(Integer.parseInt(cmd.getOptionValue(MAX_READ_COUNT, "0")));

        if (!dataCollection.loadBachelorData(bachelorInputFile))
        {
            return false;
        }

        while (dataCollection.processBachelorData())
        {
            processCurrentRecords(cmd, dataCollection.getBachelorVariants(), sampleId, sampleDirectory, bachelorDataDir);
        }

        closeBufferedWriter(mWriter);

        return true;
    }

    private void processCurrentRecords(final CommandLine cmd, List<BachelorGermlineVariant> bachRecords,
            final String sampleId, final String sampleDirectory, final String bachelorDataDir)
    {
        if (!mIsBatchMode)
        {
            if (!mAllelDepthLoader.applyPileupData(bachRecords))
            {
                return;
            }
        }

        Map<String, List<BachelorGermlineVariant>> sampleRecordsMap = Maps.newHashMap();

        if (mIsBatchMode)
        {
            String currentSample = "";
            List<BachelorGermlineVariant> sampleRecords = null;

            for (final BachelorGermlineVariant bachRecord : bachRecords)
            {
                if (currentSample.isEmpty() || !currentSample.equals(bachRecord.SampleId))
                {
                    currentSample = bachRecord.SampleId;
                    sampleRecords = sampleRecordsMap.get(bachRecord.SampleId);

                    if (sampleRecords == null)
                    {
                        sampleRecords = Lists.newArrayList();
                        sampleRecordsMap.put(bachRecord.SampleId, sampleRecords);
                    }
                }

                sampleRecords.add(bachRecord);
            }

            if (sampleRecordsMap.size() > 1)
            {
                LOGGER.info("{} unique samples after filtering", sampleRecordsMap.size());
            }
        }
        else
        {
            sampleRecordsMap.put(sampleId, bachRecords);
        }

        for (final Map.Entry<String, List<BachelorGermlineVariant>> entry : sampleRecordsMap.entrySet())
        {
            final String specificSample = entry.getKey();
            List<BachelorGermlineVariant> sampleRecords = entry.getValue();

            LOGGER.info("sample({}) processing {} germline reports", specificSample, sampleRecords.size());

            // sort by chromosome and position
            Collections.sort(sampleRecords);

            // create variant objects for VCF file writing and enrichment, and cache against bachelor record
            buildVariants(specificSample, sampleRecords);

            annotateRecords(specificSample, sampleRecords, cmd, sampleDirectory);

            int validRecordCount = 0;
            for (final BachelorGermlineVariant bachRecord : sampleRecords)
            {
                if (bachRecord.isValid())
                {
                    ++validRecordCount;
                }
            }

            if (validRecordCount == 0)
            {
                LOGGER.info("sample({}) has no valid germline reports", specificSample, validRecordCount);
                continue;
            }

            if (mUploadRecordsToDB)
            {
                LOGGER.info("sample({}) writing {} germline reports to database", specificSample, validRecordCount);
                writeToDatabase(specificSample, sampleRecords);
            }

            writeToFile(specificSample, sampleRecords, bachelorDataDir);
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        if (cmd.hasOption(LOG_DEBUG))
        {
            Configurator.setRootLevel(Level.DEBUG);
        }

        BachelorPP bachelorPP = new BachelorPP();

        if (!bachelorPP.initialise(cmd))
        {
            return;
        }

        if (!bachelorPP.run(cmd))
        {
            return;
        }

        LOGGER.info("run complete");
    }

    private void buildVariants(final String sampleId, List<BachelorGermlineVariant> bachRecords)
    {
        for (final BachelorGermlineVariant bachRecord : bachRecords)
        {
            VariantContextBuilder builder = new VariantContextBuilder();
            builder.id(bachRecord.VariantId);
            builder.loc(bachRecord.Chromosome, bachRecord.Position, bachRecord.Position + bachRecord.Ref.length() - 1);

            List<Allele> alleles = Lists.newArrayList();
            alleles.add(Allele.create(bachRecord.Ref, true));
            alleles.add(Allele.create(bachRecord.Alts, false));
            builder.alleles(alleles);

            List<Genotype> genoTypes = Lists.newArrayList();

            GenotypeBuilder gBuilder = new GenotypeBuilder(sampleId, builder.getAlleles());
            int[] adCounts = { bachRecord.getTumorRefCount(), bachRecord.getTumorAltCount() };
            gBuilder.AD(adCounts);
            gBuilder.DP(bachRecord.getGermlineReadDepth());
            genoTypes.add(gBuilder.make());

            builder.genotypes(genoTypes);
            VariantContext variantContext = builder.make();

            variantContext.getCommonInfo().addFilter("PASS");
            variantContext.getCommonInfo().putAttribute(SNPEFF_IDENTIFIER, bachRecord.Annotations);

            bachRecord.setVariantContext(variantContext);

            SomaticVariant somVariant = SomaticVariantFactory.unfilteredInstance().createSomaticVariant(sampleId, variantContext);
            bachRecord.setSomaticVariant(somVariant);
        }
    }

    private void annotateRecords(final String sampleId, List<BachelorGermlineVariant> bachRecords, final CommandLine cmd,
            final String sampleDirectory)
    {
        final PurityContext purityContext;
        final Multimap<Chromosome, PurpleCopyNumber> copyNumbers;
        final Multimap<Chromosome, FittedRegion> copyNumberRegions;

        if (cmd.hasOption(PURPLE_DATA_DIRECTORY))
        {
            final String purplePath = sampleDirectory + cmd.getOptionValue(PURPLE_DATA_DIRECTORY);

            LOGGER.debug("sample({}) loading purple data from file using path {}", sampleId, purplePath);

            try
            {
                purityContext = FittedPurityFile.read(purplePath, sampleId);

                List<PurpleCopyNumber> copyNumberData =
                        PurpleCopyNumberFile.read(PurpleCopyNumberFile.generateFilename(purplePath, sampleId));

                copyNumbers = Multimaps.fromRegions(copyNumberData);

                List<FittedRegion> fittedRegionData = FittedRegionFile.read(FittedRegionFile.generateFilename(purplePath, sampleId));

                copyNumberRegions = Multimaps.fromRegions(fittedRegionData);
            } catch (IOException e)
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
                .map(BachelorGermlineVariant::getSomaticVariant)
                .collect(Collectors.toList());

        final PurityAdjuster purityAdjuster = purityContext == null
                ? new PurityAdjuster(Gender.FEMALE, 1, 1)
                : new PurityAdjuster(purityContext.gender(), purityContext.bestFit().purity(), purityContext.bestFit().normFactor());

        final PurityAdjustedSomaticVariantFactory purityAdjustmentFactory =
                new PurityAdjustedSomaticVariantFactory(purityAdjuster, copyNumbers, copyNumberRegions);

        final List<PurityAdjustedSomaticVariant> purityAdjustedVariants = purityAdjustmentFactory.create(variants);

        final List<PurityAdjustedSomaticVariant> validVariants =
                purityAdjustedVariants.stream().filter(x -> !Double.isNaN(x.ploidy())).collect(Collectors.toList());

        final double clonalPloidy = ClonalityCutoffKernel.clonalCutoff(validVariants);

        for (PurityAdjustedSomaticVariant var : purityAdjustedVariants)
        {
            for (BachelorGermlineVariant bachRecord : bachRecords)
            {
                if (bachRecord.Chromosome.equals(var.chromosome()) && bachRecord.Position == var.position())
                {
                    double adjVaf;

                    if (bachRecord.IsHomozygous)
                    {
                        adjVaf = purityAdjuster.purityAdjustedVAFWithHomozygousNormal(
                                var.chromosome(), var.adjustedCopyNumber(), var.alleleFrequency());
                    }
                    else
                    {
                        adjVaf = purityAdjuster.purityAdjustedVAFWithHeterozygousNormal(
                                var.chromosome(), var.adjustedCopyNumber(), var.alleleFrequency());
                    }

                    if (Double.isNaN(adjVaf) || Double.isInfinite(adjVaf))
                    {
                        adjVaf = 0;
                    }

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
                if (bachRecord.Chromosome.equals(var.chromosome()) && bachRecord.Position == var.position())
                {
                    bachRecord.setEnrichedVariant(var);
                    matched = true;
                    break;
                }
            }

            if (!matched)
            {
                LOGGER.debug("sample({}) enriched variant not found: var({}) gene({}) transcript({}) chr({}) position({})",
                        sampleId, bachRecord.VariantId, bachRecord.Gene, bachRecord.TranscriptId,
                        bachRecord.Chromosome, bachRecord.Position);
            }
        }
    }

    private void writeToDatabase(final String sampleId, final List<BachelorGermlineVariant> bachRecords)
    {
        final GermlineVariantDAO germlineDAO = new GermlineVariantDAO(mDbAccess.context());
        germlineDAO.write(sampleId, bachRecords);
    }

    private void writeToFile(final String sampleId, final List<BachelorGermlineVariant> bachRecords, final String outputDir)
    {
        try
        {
            if (mWriter == null)
            {
                String outputFileName = outputDir;
                if (!outputFileName.endsWith(File.separator))
                {
                    outputFileName += File.separator;
                }

                if (!mIsBatchMode)
                {
                    outputFileName += sampleId + "_germline_variants.csv";
                }
                else
                {
                    outputFileName += "bachelor_germline_variants.csv";
                }

                mWriter = createBufferedWriter(outputFileName, false);

                mWriter.write("SampleId,Program,Chromosome,Position,Type,Ref,Alt,Gene,TranscriptId,DbsnpId,CosmicId");
                mWriter.write(",Effects,CodingEffect,VcfReadData,GermlineAltCount,GermlineReadDepth,TumorAltCount,TumorReadDepth");
                mWriter.write(",AdjCopyNumber,AdjustedVaf,HighConfidenceRegion,TrinucleotideContext,Microhomology,RepeatSequence,RepeatCount");
                mWriter.write(",HgvsProtein,HgvsCoding,Biallelic,Hotspot,Mappability,GermlineStatus,MinorAllelePloidy,Filter,CodonInfo");
                mWriter.newLine();
            }

            BufferedWriter writer = mWriter;

            for (final BachelorGermlineVariant bachRecord : bachRecords)
            {
                if (!bachRecord.isValid())
                {
                    continue;
                }

                writer.write(
                        String.format("%s,%s,%s,%d",
                                bachRecord.SampleId,
                                bachRecord.Program,
                                bachRecord.Chromosome,
                                bachRecord.Position));

                final EnrichedSomaticVariant region = bachRecord.getEnrichedVariant();

                writer.write(
                        String.format(",%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%d,%d,%d,%d",
                                region.type(),
                                region.ref(),
                                region.alt(),
                                bachRecord.Gene,
                                bachRecord.TranscriptId,
                                region.dbsnpID() == null ? "" : region.dbsnpID(),
                                region.canonicalCosmicID() == null ? "" : region.canonicalCosmicID(),
                                bachRecord.Effects,
                                bachRecord.CodingEffect,
                                bachRecord.isReadDataSet(), bachRecord.getGermlineAltCount(), bachRecord.getGermlineReadDepth(),
                                bachRecord.getTumorAltCount(), bachRecord.getTumorReadDepth()));

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
                                bachRecord.HgvsProtein,
                                bachRecord.HgvsCoding,
                                bachRecord.isBiallelic(),
                                region.hotspot(),
                                region.mappability(),
                                bachRecord.IsHomozygous ? "HOM" : "HET",
                                region.minorAllelePloidy(),
                                bachRecord.isLowScore() ? "ARTEFACT" : "PASS",
                                bachRecord.CodonInfo));

                writer.newLine();
            }
        } catch (final IOException e)
        {
            LOGGER.error("error writing to outputFile: {}", e.toString());
        }
    }

    private void loadSampleListFile(final String filename)
    {
        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine(); // skip header

            while ((line = fileReader.readLine()) != null)
            {
                String[] items = line.split(",");

                final String sampleId = items[0];

                mLimitedSampleList.add(items[0]);
            }

            LOGGER.info("loaded {} specific sample ids", mLimitedSampleList.size());

        }
        catch (IOException exception)
        {
            LOGGER.error("Failed to read sample list input CSV file({}): {}", filename, exception.toString());
        }
    }


    @NotNull
    private static Options createBasicOptions()
    {
        final Options options = new Options();
        options.addOption(SAMPLE, true, "Tumor sample");
        options.addOption(SAMPLE_PATH, true, "Typically the sample run directory and then assumes a 'bachelor' sub-directory for bachelor and mini-pileup input files");
        options.addOption(BACH_DIRECTORY, true, "Override for specific bachelor input dir, if left out then assumes in sample path & 'bachelor' sub-directory");
        options.addOption(BACH_INPUT_FILE, true, "Override for specific bachelor input file, if left out then assumes in bachelor_dir & *germline_variants.csv");
        options.addOption(REF_GENOME, true, "Path to the ref genome fasta file");
        options.addOption(HIGH_CONFIDENCE_BED, true, "Path to the high confidence bed file");
        options.addOption(PURPLE_DATA_DIRECTORY, true, "Sub-directory with sample path for purple data");
        options.addOption(MAX_READ_COUNT, true, "Optional - for buffered input file reading");
        options.addOption(SAMPLE_LIST_FILE, true, "Optional: limiting list of sample IDs to process");

        options.addOption(DB_USER, true, "Database user name");
        options.addOption(DB_PASS, true, "Database password");
        options.addOption(DB_URL, true, "Database url");

        options.addOption(WRITE_TO_DB, false, "Whether to upload records to DB");

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
