package com.hartwig.hmftools.bachelor;

import static com.hartwig.hmftools.bachelor.BachelorApplication.DEFAULT_BACH_DIRECTORY;
import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONE;
import static com.hartwig.hmftools.common.variant.CodingEffect.SPLICE;
import static com.hartwig.hmftools.common.variant.VariantType.INDEL;
import static com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationFactory.SNPEFF_IDENTIFIER;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.bachelor.types.BachelorDataCollection;
import com.hartwig.hmftools.bachelor.types.BachelorGermlineVariant;
import com.hartwig.hmftools.bachelor.types.RunDirectory;
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
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class BachelorPostProcess
{
    private DatabaseAccess mDbAccess;
    private Multimap<String, GenomeRegion> mHighConfidenceRegions;
    private IndexedFastaSequenceFile mIndexedFastaSeqFile;
    private BufferedWriter mWriter;
    private AlleleDepthLoader mAllelDepthLoader;
    private BamCountReader mBamCountReader;
    private boolean mReadBamsDirect;

    private CommandLine mCmdLineArgs;
    private String mBatchDataDir;
    private boolean mIsBatchMode;
    private boolean mUploadRecordsToDB;
    private List<BachelorGermlineVariant> mBachRecords;
    private String mBachDataDir;
    private String mSampleDataDir;
    private boolean mUsingBatchOutput;

    // config items
    public static final String REF_GENOME = "ref_genome";
    private static final String HIGH_CONFIDENCE_BED = "high_confidence_bed";
    private static final String WRITE_TO_DB = "write_to_db";
    private static final String PURPLE_DATA_DIRECTORY = "purple_data_dir"; // purple data directory within the sample fir
    private static final String BACH_DIRECTORY = "bachelor_dir"; // usually defaults to the 'bachelor' subdirectory of the sample dir
    private static final String BACH_INPUT_FILE = "bachelor_file"; // full path
    private static final String READ_BAMS_DIRECT = "bam_direct"; // skip BAM slicing and use of Mini-Pileup file reading

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    private static final String DEFAULT_BACH_INPUT_FILE = "bachelor_output.csv";


    private static final Logger LOGGER = LogManager.getLogger(BachelorPostProcess.class);

    public BachelorPostProcess()
    {
        mCmdLineArgs = null;
        mBachRecords = Lists.newArrayList();
        mBatchDataDir = "";
        mBachDataDir = "";
        mSampleDataDir = "";
        mIsBatchMode = false;

        mDbAccess = null;
        mHighConfidenceRegions = null;
        mIndexedFastaSeqFile = null;
        mIsBatchMode = false;
        mUploadRecordsToDB = false;
        mAllelDepthLoader = null;
        mBamCountReader = null;
        mAllelDepthLoader = null;
        mReadBamsDirect = false;
        mWriter = null;
    }

    public boolean initialise(final CommandLine cmd, boolean isBatchMode, final String batchOutputDir)
    {
        mCmdLineArgs = cmd;
        mIsBatchMode = isBatchMode;

        try
        {
            mDbAccess = databaseAccess(cmd);
        }
        catch (SQLException e)
        {
            LOGGER.error("DB connection failed: {}", e.toString());
            return false;
        }

        mBatchDataDir = batchOutputDir;
        mUsingBatchOutput = mIsBatchMode && !mBatchDataDir.isEmpty();

        mReadBamsDirect = cmd.hasOption(READ_BAMS_DIRECT);

        if(mReadBamsDirect)
        {
            mBamCountReader = new BamCountReader();
        }
        else
        {
            mAllelDepthLoader = new AlleleDepthLoader();
        }

        if (cmd.hasOption(HIGH_CONFIDENCE_BED) && cmd.hasOption(REF_GENOME))
        {
            final String highConfidenceBed = cmd.getOptionValue(HIGH_CONFIDENCE_BED);
            final String refGenomeFile = cmd.getOptionValue(REF_GENOME);

            try
            {
                LOGGER.debug("reading high confidence bed file");
                mHighConfidenceRegions = BEDFileLoader.fromBedFile(highConfidenceBed);

                LOGGER.debug("loading indexed fasta reference file");
                mIndexedFastaSeqFile = new IndexedFastaSequenceFile(new File(refGenomeFile));

                mBamCountReader.initialise(cmd, mIndexedFastaSeqFile);
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

    private void loadBachelorRecords()
    {
        String bachelorInputFile;
        if (mCmdLineArgs.hasOption(BACH_INPUT_FILE))
        {
            bachelorInputFile = mCmdLineArgs.getOptionValue(BACH_INPUT_FILE);
            LOGGER.info("loading specific input file: {}", bachelorInputFile);
        }
        else
        {
            bachelorInputFile = mBachDataDir + DEFAULT_BACH_INPUT_FILE;
            LOGGER.info("loading sample default file: {}", bachelorInputFile);
        }

        BachelorDataCollection dataCollection = new BachelorDataCollection();

        if (!dataCollection.loadBachelorData(bachelorInputFile))
        {
            return;
        }

        mBachRecords.addAll(dataCollection.getBachelorVariants());
    }

    public void run(RunDirectory runDir, final String sampleId, @Nullable  List<BachelorGermlineVariant> bachRecords)
    {
        if(mUsingBatchOutput)
        {
            mBachDataDir = mBatchDataDir;
            mSampleDataDir = mBachDataDir;
        }
        else
        {
            mSampleDataDir = runDir.sampleDir().toString() + File.separator;
            mBachDataDir = mSampleDataDir + DEFAULT_BACH_DIRECTORY + File.separator;
        }

        if(bachRecords == null)
        {
            loadBachelorRecords();
        }
        else
        {
            mBachRecords = bachRecords;
        }

        processCurrentRecords(mBachRecords);

        if(!mUsingBatchOutput)
        {
            closeBufferedWriter(mWriter);
        }
    }

    private void processCurrentRecords(List<BachelorGermlineVariant> bachRecords)
    {
        if(mReadBamsDirect)
        {
            mBamCountReader.readBamCounts(bachRecords, mBachDataDir);
        }
        else
        {

            if (!mAllelDepthLoader.loadMiniPileupData(mBachDataDir))
                return;

            if (!mAllelDepthLoader.applyPileupData(bachRecords))
                return;
        }

        Map<String, List<BachelorGermlineVariant>> sampleRecordsMap = Maps.newHashMap();

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

        for (final Map.Entry<String, List<BachelorGermlineVariant>> entry : sampleRecordsMap.entrySet())
        {
            final String specificSample = entry.getKey();
            sampleRecords = entry.getValue();

            LOGGER.info("sample({}) processing {} germline reports", specificSample, sampleRecords.size());

            // sort by chromosome and position
            Collections.sort(sampleRecords);

            // create variant objects for VCF file writing and enrichment, and cache against bachelor record
            buildVariants(specificSample, sampleRecords);

            annotateRecords(specificSample, sampleRecords);

            filterRecords(sampleRecords);

            if (sampleRecords.isEmpty())
            {
                LOGGER.info("sample({}) has no valid germline reports", specificSample);
                continue;
            }

            if (mUploadRecordsToDB)
            {
                LOGGER.info("sample({}) writing {} germline reports to database", specificSample, sampleRecords.size());
                writeToDatabase(specificSample, sampleRecords);
            }

            writeToFile(specificSample, sampleRecords);
        }
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

    private void annotateRecords(final String sampleId, List<BachelorGermlineVariant> bachRecords)
    {
        final PurityContext purityContext;
        final Multimap<Chromosome, PurpleCopyNumber> copyNumbers;
        final Multimap<Chromosome, FittedRegion> copyNumberRegions;

        if (mCmdLineArgs.hasOption(PURPLE_DATA_DIRECTORY))
        {
            final String purplePath = mSampleDataDir + mCmdLineArgs.getOptionValue(PURPLE_DATA_DIRECTORY);

            LOGGER.debug("sample({}) loading purple data from file using path {}", sampleId, purplePath);

            try
            {
                purityContext = FittedPurityFile.read(purplePath, sampleId);

                List<PurpleCopyNumber> copyNumberData =
                        PurpleCopyNumberFile.read(PurpleCopyNumberFile.generateFilename(purplePath, sampleId));

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
                mIndexedFastaSeqFile,
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

    private static int INDEL_REPEAT_LIMIT = 8;

    private void filterRecords(final List<BachelorGermlineVariant> bachRecords)
    {
        // currently the only filter is on INDELs with either microhomology matching the gain or loss, or with high repeat count
        int index = 0;
        while(index < bachRecords.size())
        {
            BachelorGermlineVariant bachRecord = bachRecords.get(index);

            if (!bachRecord.isValid())
            {
                bachRecords.remove(index);
                continue;
            }

            final EnrichedSomaticVariant enrichedVariant = bachRecord.getEnrichedVariant();

            if(enrichedVariant.type() == INDEL && (bachRecord.CodingEffect == SPLICE || bachRecord.CodingEffect == NONE))
            {
                int repeatCount = enrichedVariant.repeatCount();

                if(repeatCount > INDEL_REPEAT_LIMIT)
                {
                    LOGGER.debug("filtered var({}) indel {} with high repeatCount({})",
                            bachRecord.asString(), bachRecord.CodingEffect, repeatCount);
                    bachRecords.remove(index);
                    continue;
                }

                final String microhomology = enrichedVariant.microhomology();
                final String ref = bachRecord.Ref;
                final String alt = bachRecord.Alts;

                String mergeStr1;
                String mergeStr2;
                String compareStr;
                if(alt.length() > ref.length())
                {
                    mergeStr1 = ref + microhomology;
                    mergeStr2 = microhomology + ref;
                    compareStr = alt;
                }
                else
                {
                    mergeStr1 = alt + microhomology;
                    mergeStr2 = microhomology + alt;
                    compareStr = ref;
                }

                if(compareStr.equals(mergeStr1) || compareStr.equals(mergeStr2))
                {
                    LOGGER.debug("filtered var({}) indel {} with ref, alt and microHom equal",
                            bachRecord.asString(), bachRecord.CodingEffect, repeatCount);
                    bachRecords.remove(index);
                    continue;
                }
            }

            ++index;
        }
    }

    private void writeToDatabase(final String sampleId, final List<BachelorGermlineVariant> bachRecords)
    {
        final GermlineVariantDAO germlineDAO = new GermlineVariantDAO(mDbAccess.context());
        germlineDAO.write(sampleId, bachRecords);
    }

    private void writeToFile(final String sampleId, final List<BachelorGermlineVariant> bachRecords)
    {
        try
        {
            if (mWriter == null)
            {
                String outputFileName = mBachDataDir;
                if (!outputFileName.endsWith(File.separator))
                {
                    outputFileName += File.separator;
                }

                if (!mIsBatchMode || mBatchDataDir.isEmpty())
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
                mWriter.write(",ClinvarMatch,ClinvarSignificance,ClinvarSigInfo");
                mWriter.newLine();
            }

            BufferedWriter writer = mWriter;

            for (final BachelorGermlineVariant bachRecord : bachRecords)
            {
                if (!bachRecord.isValid())
                {
                    continue;
                }

                writer.write(String.format("%s,%s,%s,%d",
                        bachRecord.SampleId, bachRecord.Program, bachRecord.Chromosome, bachRecord.Position));

                final EnrichedSomaticVariant enrichedVariant = bachRecord.getEnrichedVariant();

                writer.write(String.format(",%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%d,%d,%d,%d",
                        enrichedVariant.type(), enrichedVariant.ref(), enrichedVariant.alt(),
                        bachRecord.Gene, bachRecord.TranscriptId,
                        enrichedVariant.dbsnpID() == null ? "" : enrichedVariant.dbsnpID(),
                        enrichedVariant.canonicalCosmicID() == null ? "" : enrichedVariant.canonicalCosmicID(),
                        bachRecord.Effects, bachRecord.CodingEffect,
                        bachRecord.isReadDataSet(), bachRecord.getGermlineAltCount(), bachRecord.getGermlineReadDepth(),
                        bachRecord.getTumorAltCount(), bachRecord.getTumorReadDepth()));

                writer.write(String.format(",%.2f,%.2f,%s,%s,%s,%s,%d",
                        enrichedVariant.adjustedCopyNumber(), bachRecord.getAdjustedVaf(), enrichedVariant.highConfidenceRegion(),
                        enrichedVariant.trinucleotideContext(), enrichedVariant.microhomology(), enrichedVariant.repeatSequence(),
                        enrichedVariant.repeatCount()));

                writer.write(String.format(",%s,%s,%s,%s,%s,%s,%.2f,%s,%s,%s,%s,%s",
                        bachRecord.HgvsProtein, bachRecord.HgvsCoding, bachRecord.isBiallelic(), enrichedVariant.hotspot(),
                        enrichedVariant.mappability(), bachRecord.IsHomozygous ? "HOM" : "HET", enrichedVariant.minorAllelePloidy(),
                        bachRecord.isLowScore() ? "ARTEFACT" : "PASS", bachRecord.CodonInfo,
                        bachRecord.getClinvarMatch(), bachRecord.getClinvarSig(), bachRecord.getClinvarSigInfo()));

                writer.newLine();
            }
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to outputFile: {}", e.toString());
        }
    }

    public void close()
    {
        if(mUsingBatchOutput)
        {
            closeBufferedWriter(mWriter);
        }
    }

    public static void addCmdLineOptions(Options options)
    {
        options.addOption(BACH_DIRECTORY, true, "Override for specific bachelor input dir, if left out then assumes in sample path & 'bachelor' sub-directory");
        options.addOption(BACH_INPUT_FILE, true, "Override for specific bachelor input file, if left out then assumes in bachelor_dir & *germline_variants.csv");
        options.addOption(REF_GENOME, true, "Path to the ref genome fasta file");
        options.addOption(HIGH_CONFIDENCE_BED, true, "Path to the high confidence bed file");
        options.addOption(PURPLE_DATA_DIRECTORY, true, "Sub-directory with sample path for purple data");
        options.addOption(READ_BAMS_DIRECT, false, "Read tumor alt and read depth from available BAM file");

        options.addOption(DB_USER, true, "Database user name");
        options.addOption(DB_PASS, true, "Database password");
        options.addOption(DB_URL, true, "Database url");

        options.addOption(WRITE_TO_DB, false, "Whether to upload records to DB");
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
