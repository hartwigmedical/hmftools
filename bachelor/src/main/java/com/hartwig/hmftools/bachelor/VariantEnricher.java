package com.hartwig.hmftools.bachelor;

import static com.hartwig.hmftools.bachelor.types.BachelorConfig.BACH_LOGGER;
import static com.hartwig.hmftools.bachelor.types.BachelorConfig.BATCH_FILE;
import static com.hartwig.hmftools.bachelor.types.BachelorConfig.DB_URL;
import static com.hartwig.hmftools.bachelor.types.BachelorConfig.REF_GENOME;
import static com.hartwig.hmftools.bachelor.types.BachelorConfig.databaseAccess;
import static com.hartwig.hmftools.bachelor.types.FilterType.ARTEFACT;
import static com.hartwig.hmftools.bachelor.types.FilterType.PASS;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONE;
import static com.hartwig.hmftools.common.variant.CodingEffect.SPLICE;
import static com.hartwig.hmftools.common.variant.VariantType.INDEL;
import static com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationFactory.SNPEFF_IDENTIFIER;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.bachelor.types.BachelorConfig;
import com.hartwig.hmftools.bachelor.types.BachelorGermlineVariant;
import com.hartwig.hmftools.bachelor.types.EnrichedSomaticVariant;
import com.hartwig.hmftools.bachelor.types.EnrichedSomaticVariantFactory;
import com.hartwig.hmftools.bachelor.types.GermlineVariant;
import com.hartwig.hmftools.bachelor.types.GermlineVariantFile;
import com.hartwig.hmftools.bachelor.types.ImmutableGermlineVariant;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariant;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariantFactory;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

class VariantEnricher
{
    private final DatabaseAccess mDbAccess;
    private final BachelorConfig mConfig;
    private IndexedFastaSequenceFile mIndexedFastaSeqFile;
    private BufferedWriter mWriter;
    private BamCountReader mBamCountReader;

    private List<BachelorGermlineVariant> mBachRecords;

    VariantEnricher(final BachelorConfig config, final CommandLine cmd)
    {
        mConfig = config;
        mBachRecords = Lists.newArrayList();

        mDbAccess = cmd.hasOption(DB_URL) ? databaseAccess(cmd) : null;

        mIndexedFastaSeqFile = null;
        mBamCountReader = null;

        if (cmd.hasOption(REF_GENOME))
        {
            final String refGenomeFile = cmd.getOptionValue(REF_GENOME);

            try
            {
                BACH_LOGGER.debug("Loading indexed fasta reference file");
                mIndexedFastaSeqFile = new IndexedFastaSequenceFile(new File(refGenomeFile));

                mBamCountReader = new BamCountReader();
                mBamCountReader.initialise(config.RefGenomeFile, mIndexedFastaSeqFile);
            }
            catch (IOException e)
            {
                BACH_LOGGER.error("Reference file loading failed");
                return;
            }
        }

        mWriter = null;
    }

    void run(@Nullable List<BachelorGermlineVariant> bachRecords)
    {
        mBachRecords = bachRecords;

        processCurrentRecords(mBachRecords);

        if(!mConfig.IsBatchMode)
        {
            closeBufferedWriter(mWriter);
            mWriter = null;
        }
    }

    private void processCurrentRecords(@Nullable List<BachelorGermlineVariant> bachRecords)
    {
        if(bachRecords.isEmpty() && !mConfig.IsBatchMode)
        {
            writeToFile(mConfig.SampleId, Lists.newArrayList());
            return;
        }

        long recordsWithTumorData = bachRecords.stream().filter(x -> x.isReadDataSet()).count();

        if(mBamCountReader != null && recordsWithTumorData < bachRecords.size())
        {
            mBamCountReader.readBamCounts(mConfig.BamFile, bachRecords);
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

            BACH_LOGGER.info("sample({}) processing {} germline reports", specificSample, sampleRecords.size());

            // sort by chromosome and position
            Collections.sort(sampleRecords);

            // create variant objects for VCF file writing and enrichment, and cache against bachelor record
            buildVariants(specificSample, sampleRecords);

            if(!mConfig.SkipEnrichment)
            {
                annotateRecords(specificSample, sampleRecords);
            }

            filterRecords(sampleRecords);

            if (sampleRecords.isEmpty())
            {
                BACH_LOGGER.info("sample({}) has no valid germline reports", specificSample);
                writeToFile(specificSample, Lists.newArrayList());
                continue;
            }

            final List<GermlineVariant> germlineVariants = convert(sampleRecords);

            if (mDbAccess != null)
            {
                BACH_LOGGER.info("sample({}) writing {} germline reports to database", specificSample, sampleRecords.size());
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
            BACH_LOGGER.debug("sample({}) loading purple data from file using path {}", sampleId, mConfig.PurpleDataDir);

            try
            {
                purityContext = FittedPurityFile.read(mConfig.PurpleDataDir, sampleId);
                copyNumbers = PurpleCopyNumberFile.read(PurpleCopyNumberFile.generateFilenameForReading(mConfig.PurpleDataDir, sampleId));
            }
            catch (IOException e)
            {
                BACH_LOGGER.error("Failed to read purple data from {}: {}", mConfig.PurpleDataDir, e.toString());
                return;
            }
        }
        else
        {
            BACH_LOGGER.debug("sample({}) loading purple data from database", sampleId);

            purityContext = mDbAccess.readPurityContext(sampleId);

            if (purityContext == null)
            {
                BACH_LOGGER.warn("Failed to read purity data");
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

        BACH_LOGGER.debug("sample({}) enriching variants", sampleId);

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
                BACH_LOGGER.debug("sample({}) enriched variant not found: var({}) gene({}) transcript({}) chr({}) position({})",
                        sampleId, bachRecord.VariantId, bachRecord.Gene, bachRecord.TranscriptId,
                        bachRecord.Chromosome, bachRecord.Position);
            }
        }
    }

    private static final int HIGH_INDEL_REPEAT_COUNT = 8;

    private void filterRecords(final List<BachelorGermlineVariant> bachRecords)
    {
        // currently the only filter is on INDELs with either microhomology matching the gain or loss, or with high repeat count
        int index = 0;
        while(index < bachRecords.size())
        {
            BachelorGermlineVariant bachRecord = bachRecords.get(index);

            if (!bachRecord.isValid(false))
            {
                bachRecords.remove(index);
                continue;
            }

            final SomaticVariant somaticVariant = bachRecord.getSomaticVariant();

            if(somaticVariant.type() == INDEL && (bachRecord.CodingEffect == SPLICE || bachRecord.CodingEffect == NONE))
            {
                int repeatCount = somaticVariant.repeatCount();

                if(repeatCount > HIGH_INDEL_REPEAT_COUNT)
                {
                    BACH_LOGGER.debug("Filtered var({}) indel {} with high repeatCount({})",
                            bachRecord.asString(), bachRecord.CodingEffect, repeatCount);
                    bachRecords.remove(index);
                    continue;
                }

                final String microhomology = somaticVariant.microhomology();
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
                    BACH_LOGGER.debug("Filtered var({}) indel {} with ref, alt and microHom equal",
                            bachRecord.asString(), bachRecord.CodingEffect);
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

    private List<GermlineVariant> convert(final List<BachelorGermlineVariant> bachRecords)
    {
        final List<GermlineVariant> germlineVariants = Lists.newArrayList();

        for(final BachelorGermlineVariant bachRecord : bachRecords)
        {
            if (!bachRecord.isValid(!mConfig.SkipEnrichment))
                continue;

            final SomaticVariant somaticVariant = bachRecord.getSomaticVariant();
            final EnrichedSomaticVariant enrichedVariant = bachRecord.getEnrichedVariant();

            germlineVariants.add(ImmutableGermlineVariant.builder()
                    .chromosome(bachRecord.Chromosome)
                    .position(bachRecord.Position)
                    .filter(bachRecord.filterType())
                    .type(somaticVariant.type().toString())
                    .ref(somaticVariant.ref())
                    .alts(somaticVariant.alt())
                    .gene(bachRecord.Gene)
                    .effects(bachRecord.Effects)
                    .codingEffect(bachRecord.CodingEffect)
                    .transcriptId(bachRecord.TranscriptId)
                    .alleleReadCount(bachRecord.getTumorAltCount())
                    .totalReadCount(bachRecord.getTumorReadDepth())
                    .adjustedCopyNumber(enrichedVariant != null ? enrichedVariant.adjustedCopyNumber() : somaticVariant.adjustedCopyNumber())
                    .adjustedVaf(bachRecord.getAdjustedVaf())
                    .trinucleotideContext(somaticVariant.trinucleotideContext())
                    .microhomology(somaticVariant.microhomology())
                    .repeatSequence(somaticVariant.repeatSequence())
                    .repeatCount(somaticVariant.repeatCount())
                    .hgvsProtein(bachRecord.HgvsProtein)
                    .hgvsCoding(bachRecord.HgvsCoding)
                    .biallelic(bachRecord.isBiallelic())
                    .minorAlleleJcn(enrichedVariant != null ? enrichedVariant.minorAlleleCopyNumber() : somaticVariant.minorAlleleCopyNumber())
                    .program(bachRecord.Program)
                    .variantId(bachRecord.VariantId)
                    .annotations(bachRecord.Annotations)
                    .phredScore(bachRecord.PhredScore)
                    .isHomozygous(bachRecord.IsHomozygous)
                    .matchType(bachRecord.getMatchType())
                    .codonInfo(bachRecord.CodonInfo)
                    .clinvarMatch(bachRecord.getClinvarMatch())
                    .clinvarSignificance(bachRecord.getClinvarSig())
                    .clinvarSignificanceInfo(bachRecord.getClinvarSigInfo())
                    .build());
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
            BACH_LOGGER.error("Error writing to outputFile: {}", e.toString());
        }
    }

    void close()
    {
        if(mConfig.IsBatchMode)
        {
            closeBufferedWriter(mWriter);
        }
    }
}
