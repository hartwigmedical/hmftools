package com.hartwig.hmftools.bachelor;

import static com.hartwig.hmftools.bachelor.LoadGermlineVariants.DB_URL;
import static com.hartwig.hmftools.bachelor.types.BachelorConfig.BATCH_FILE;
import static com.hartwig.hmftools.bachelor.types.BachelorConfig.INTERIM_FILENAME;
import static com.hartwig.hmftools.bachelor.types.BachelorConfig.READ_BAMS_DIRECT;
import static com.hartwig.hmftools.bachelor.types.BachelorConfig.REF_GENOME;
import static com.hartwig.hmftools.bachelor.types.BachelorConfig.databaseAccess;
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
import com.hartwig.hmftools.bachelor.types.BachelorConfig;
import com.hartwig.hmftools.bachelor.types.BachelorDataCollection;
import com.hartwig.hmftools.bachelor.types.BachelorGermlineVariant;
import com.hartwig.hmftools.bachelor.types.GermlineVariant;
import com.hartwig.hmftools.bachelor.types.GermlineVariantFile;
import com.hartwig.hmftools.bachelor.types.ImmutableGermlineVariant;
import com.hartwig.hmftools.bachelor.types.RunDirectory;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariantFactory;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariant;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariantFactory;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class BachelorPostProcess
{
    private final DatabaseAccess mDbAccess;
    private final BachelorConfig mConfig;
    private IndexedFastaSequenceFile mIndexedFastaSeqFile;
    private BufferedWriter mWriter;
    private final AlleleDepthLoader mAllelDepthLoader;
    private BamCountReader mBamCountReader;
    private final boolean mReadBamsDirect;

    private List<BachelorGermlineVariant> mBachRecords;

    private static final Logger LOGGER = LogManager.getLogger(BachelorPostProcess.class);

    BachelorPostProcess(final BachelorConfig config, final CommandLine cmd)
    {
        mConfig = config;
        mBachRecords = Lists.newArrayList();

        mReadBamsDirect = cmd.hasOption(READ_BAMS_DIRECT);

        mBamCountReader = mReadBamsDirect ? new BamCountReader() : null;
        mAllelDepthLoader = !mReadBamsDirect ? new AlleleDepthLoader() : null;

        mDbAccess = cmd.hasOption(DB_URL) ? databaseAccess(cmd) : null;

        mIndexedFastaSeqFile = null;

        if (cmd.hasOption(REF_GENOME))
        {
            final String refGenomeFile = cmd.getOptionValue(REF_GENOME);

            try
            {
                LOGGER.debug("Loading indexed fasta reference file");
                mIndexedFastaSeqFile = new IndexedFastaSequenceFile(new File(refGenomeFile));

                if(mBamCountReader != null)
                    mBamCountReader.initialise(cmd, mIndexedFastaSeqFile);
            }
            catch (IOException e)
            {
                LOGGER.error("Reference file loading failed");
                return;
            }
        }

        mWriter = null;
    }

    private void loadBachelorRecords()
    {
        String bachelorInputFile = mConfig.SampleDataDir;

        if(mConfig.IsBatchMode)
        {
            bachelorInputFile += BATCH_FILE + INTERIM_FILENAME;
        }
        else
        {
            bachelorInputFile += mConfig.SampleId + INTERIM_FILENAME;
        }

        LOGGER.info("Loading bachelor interim file: {}", bachelorInputFile);

        BachelorDataCollection dataCollection = new BachelorDataCollection();

        if (!dataCollection.loadBachelorData(bachelorInputFile))
        {
            return;
        }

        mBachRecords.addAll(dataCollection.getBachelorVariants());
    }

    public void run(@Nullable List<BachelorGermlineVariant> bachRecords)
    {
        if(bachRecords == null)
        {
            loadBachelorRecords();
        }
        else
        {
            mBachRecords = bachRecords;
        }

        processCurrentRecords(mBachRecords);

        if(!mConfig.IsBatchMode)
        {
            closeBufferedWriter(mWriter);
            mWriter = null;
        }
    }

    private void processCurrentRecords(List<BachelorGermlineVariant> bachRecords)
    {
        if(bachRecords.isEmpty() && !mConfig.IsBatchMode)
        {
            writeToFile(mConfig.SampleId, Lists.newArrayList());
            return;
        }

        long recordsWithTumorData = bachRecords.stream().filter(x -> x.isReadDataSet()).count();

        if(recordsWithTumorData < bachRecords.size())
        {
            if (mReadBamsDirect)
            {
                mBamCountReader.readBamCounts(bachRecords, mConfig.SampleDataDir);
            }
            else
            {
                if (!mAllelDepthLoader.loadMiniPileupData(mConfig.SampleDataDir))
                    return;

                if (!mAllelDepthLoader.applyPileupData(bachRecords))
                    return;
            }
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

            assert sampleRecords != null;
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
                writeToFile(specificSample, Lists.newArrayList());
                continue;
            }

            final List<GermlineVariant> germlineVariants = convert(sampleRecords);

            if (mDbAccess != null)
            {
                LOGGER.info("sample({}) writing {} germline reports to database", specificSample, sampleRecords.size());
                writeToDatabase(specificSample, germlineVariants);
            }

            writeToFile(specificSample, germlineVariants);
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
        final List<PurpleCopyNumber> copyNumbers;

        if (!mConfig.PurpleDataDir.isEmpty())
        {
            LOGGER.debug("sample({}) loading purple data from file using path {}", sampleId, mConfig.PurpleDataDir);

            try
            {
                purityContext = FittedPurityFile.read(mConfig.PurpleDataDir, sampleId);
                copyNumbers = PurpleCopyNumberFile.read(PurpleCopyNumberFile.generateFilenameForReading(mConfig.PurpleDataDir, sampleId));
            }
            catch (IOException e)
            {
                LOGGER.error("Failed to read purple data from {}: {}", mConfig.PurpleDataDir, e.toString());
                return;
            }
        }
        else
        {
            LOGGER.debug("sample({}) loading purple data from database", sampleId);

            purityContext = mDbAccess.readPurityContext(sampleId);

            if (purityContext == null)
            {
                LOGGER.warn("Failed to read purity data");
            }

            copyNumbers = mDbAccess.readCopynumbers(sampleId);
        }

        List<SomaticVariant> variants = bachRecords.stream()
                .filter(x -> x.getSomaticVariant() != null)
                .map(BachelorGermlineVariant::getSomaticVariant)
                .collect(Collectors.toList());

        final PurityAdjuster purityAdjuster = purityContext == null
                ? new PurityAdjuster(Gender.FEMALE, 1, 1)
                : new PurityAdjuster(purityContext.gender(), purityContext.bestFit().purity(), purityContext.bestFit().normFactor());

        final PurityAdjustedSomaticVariantFactory purityAdjustmentFactory =
                new PurityAdjustedSomaticVariantFactory(sampleId, purityAdjuster, copyNumbers);

        final List<PurityAdjustedSomaticVariant> purityAdjustedVariants = purityAdjustmentFactory.create(variants);

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

        final EnrichedSomaticVariantFactory enrichedSomaticVariantFactory = new EnrichedSomaticVariantFactory(mIndexedFastaSeqFile);

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

                if(repeatCount > 8)
                {
                    LOGGER.debug("Filtered var({}) indel {} with high repeatCount({})",
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
                    LOGGER.debug("Filtered var({}) indel {} with ref, alt and microHom equal",
                            bachRecord.asString(), bachRecord.CodingEffect, repeatCount);
                    bachRecords.remove(index);
                    continue;
                }
            }

            ++index;
        }
    }

    private void writeToDatabase(final String sampleId, final List<GermlineVariant> germlineVariants)
    {
        final GermlineVariantDAO germlineDAO = new GermlineVariantDAO(mDbAccess.context());
        germlineDAO.write(sampleId, germlineVariants);
    }

    private final List<GermlineVariant> convert(final List<BachelorGermlineVariant> bachRecords)
    {
        final List<GermlineVariant> germlineVariants = Lists.newArrayList();

        for(final BachelorGermlineVariant bachRecord : bachRecords)
        {
            if (!bachRecord.isValid())
                continue;

            final EnrichedSomaticVariant enrichedVariant = bachRecord.getEnrichedVariant();

            germlineVariants.add(ImmutableGermlineVariant.builder()
                    .chromosome(bachRecord.Chromosome)
                    .position(bachRecord.Position)
                    .filter(bachRecord.isLowScore() ? "ARTEFACT" : "PASS")
                    .type(enrichedVariant.type().toString())
                    .ref(enrichedVariant.ref())
                    .alts(enrichedVariant.alt())
                    .gene(bachRecord.Gene)
                    .cosmicId(enrichedVariant.canonicalCosmicID() == null ? "" : enrichedVariant.canonicalCosmicID())
                    .dbsnpId(enrichedVariant.dbsnpID() == null ? "" : enrichedVariant.dbsnpID())
                    .effects(bachRecord.Effects)
                    .codingEffect(bachRecord.CodingEffect)
                    .transcriptId(bachRecord.TranscriptId)
                    .alleleReadCount(bachRecord.getTumorAltCount())
                    .totalReadCount(bachRecord.getTumorReadDepth())
                    .adjustedCopyNumber(enrichedVariant.adjustedCopyNumber())
                    .adjustedVaf(bachRecord.getAdjustedVaf())
                    .trinucleotideContext(enrichedVariant.trinucleotideContext())
                    .microhomology(enrichedVariant.microhomology())
                    .repeatSequence(enrichedVariant.repeatSequence())
                    .repeatCount(enrichedVariant.repeatCount())
                    .hgvsProtein(bachRecord.HgvsProtein)
                    .hgvsCoding(bachRecord.HgvsCoding)
                    .biallelic(bachRecord.isBiallelic())
                    .hotspot(enrichedVariant.isHotspot())
                    .mappability(enrichedVariant.mappability())
                    .minorAllelePloidy(enrichedVariant.minorAllelePloidy())
                    .program(bachRecord.Program)
                    .variantId(bachRecord.VariantId)
                    .annotations(bachRecord.Annotations)
                    .phredScore(bachRecord.PhredScore)
                    .isHomozygous(bachRecord.IsHomozygous)
                    .matchType(bachRecord.MatchType)
                    .codonInfo(bachRecord.CodonInfo)
                    .clinvarMatch(bachRecord.getClinvarMatch())
                    .clinvarSignificance(bachRecord.getClinvarSig())
                    .clinvarSignificanceInfo(bachRecord.getClinvarSigInfo())
                    .build());

            /*
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
             */


        }

        return germlineVariants;
    }

    private void writeToFile(final String sampleId, final List<GermlineVariant> germlineVariants)
    {
        try
        {
            if(mConfig.IsBatchMode)
            {
                if (mWriter == null)
                {
                    final String filename = GermlineVariantFile.generateFilename(mConfig.OutputDir, BATCH_FILE);
                    mWriter = createBufferedWriter(filename, false);
                }

                for(GermlineVariant variant : germlineVariants)
                {
                    mWriter.write(GermlineVariantFile.toString(variant));
                    mWriter.newLine();
                }
            }
            else
            {
                final String filename = GermlineVariantFile.generateFilename(mConfig.OutputDir, sampleId);
                GermlineVariantFile.write(filename, germlineVariants);
            }
        }
        catch (final IOException e)
        {
            LOGGER.error("Error writing to outputFile: {}", e.toString());
        }
    }

    public void close()
    {
        if(mConfig.IsBatchMode)
        {
            closeBufferedWriter(mWriter);
        }
    }

}
