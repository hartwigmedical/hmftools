package com.hartwig.hmftools.bachelor;

import static com.hartwig.hmftools.bachelor.types.BachelorConfig.BACH_LOGGER;
import static com.hartwig.hmftools.bachelor.types.BachelorConfig.REF_GENOME;
import static com.hartwig.hmftools.bachelor.types.PathogenicType.BLACK_LIST;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONE;
import static com.hartwig.hmftools.common.variant.CodingEffect.SPLICE;
import static com.hartwig.hmftools.common.variant.VariantHeader.PURPLE_CN_INFO;
import static com.hartwig.hmftools.common.variant.VariantHeader.PURPLE_MINOR_ALLELE_CN_INFO;
import static com.hartwig.hmftools.common.variant.VariantType.INDEL;
import static com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment.MICROHOMOLOGY_FLAG;
import static com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment.REPEAT_COUNT_FLAG;
import static com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment.REPEAT_SEQUENCE_FLAG;
import static com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment.TRINUCLEOTIDE_FLAG;
import static com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationFactory.SNPEFF_IDENTIFIER;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.bachelor.types.BachelorConfig;
import com.hartwig.hmftools.bachelor.types.BachelorGermlineVariant;
import com.hartwig.hmftools.bachelor.types.GermlineVariant;
import com.hartwig.hmftools.bachelor.types.GermlineVariantFile;
import com.hartwig.hmftools.bachelor.types.ImmutableGermlineVariant;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurityAdjusterTypicalChromosome;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.purity.PurityContextFile;
import com.hartwig.hmftools.common.variant.PurityAdjustedSomaticVariantFactory;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.germline.ImmutableReportableGermlineVariant;
import com.hartwig.hmftools.common.variant.germline.ReportableGermlineVariant;
import com.hartwig.hmftools.common.variant.germline.ReportableGermlineVariantFile;
import com.hartwig.hmftools.common.variant.repeat.RepeatContext;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.CommonInfo;
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
    private final RefGenomeEnrichment mRefGenomeEnrichment;

    private List<BachelorGermlineVariant> mBachRecords;

    VariantEnricher(final BachelorConfig config, final CommandLine cmd)
    {
        mConfig = config;
        mBachRecords = Lists.newArrayList();

        mDbAccess = createDatabaseAccess(cmd);

        mIndexedFastaSeqFile = null;
        mBamCountReader = null;

        if (cmd.hasOption(REF_GENOME))
        {
            final String refGenomeFile = cmd.getOptionValue(REF_GENOME);

            try
            {
                BACH_LOGGER.debug("Loading indexed fasta reference file");
                mIndexedFastaSeqFile = new IndexedFastaSequenceFile(new File(refGenomeFile));

                if (mConfig.BamFile != null)
                {
                    mBamCountReader = new BamCountReader();
                    mBamCountReader.initialise(config.RefGenomeFile, mIndexedFastaSeqFile);
                }
            }
            catch (IOException e)
            {
                BACH_LOGGER.error("Reference file loading failed");
            }
        }

        mRefGenomeEnrichment = new RefGenomeEnrichment(mIndexedFastaSeqFile);

        mWriter = null;
    }

    void run(@Nullable List<BachelorGermlineVariant> bachRecords)
    {
        mBachRecords = bachRecords;

        processCurrentRecords(mBachRecords);
        closeBufferedWriter(mWriter);
        mWriter = null;
    }

    private void processCurrentRecords(@Nullable List<BachelorGermlineVariant> bachRecords)
    {
        if(bachRecords.isEmpty())
        {
            // write empty output files even if no variants are found
            writeToFile(mConfig.SampleId, Lists.newArrayList());
            return;
        }

        if(mBamCountReader != null)
            mBamCountReader.readBamCounts(mConfig.BamFile, bachRecords);

        final String sampleId = mConfig.SampleId;

        try
        {
            // sort by chromosome and position
            Collections.sort(bachRecords);
        }
        catch(Exception e)
        {
            BACH_LOGGER.error("sorting error: {}", e.toString());
        }

        BACH_LOGGER.info("Sample({}) processing {} germline reports", sampleId, bachRecords.size());

        buildVariants(sampleId, bachRecords);

        if(!mConfig.PurpleDataDir.isEmpty())
        {
            addPurpleData(sampleId, bachRecords);
        }

        applyIndelFilter(bachRecords);

        final List<GermlineVariant> germlineVariants = convert(bachRecords);

        if (mDbAccess != null)
        {
            BACH_LOGGER.info("Sample({}) writing {} germline reports to database", sampleId, bachRecords.size());
            writeToDatabase(sampleId, germlineVariants);
        }

        writeToFile(sampleId, germlineVariants);
    }

    private void buildVariants(final String sampleId, List<BachelorGermlineVariant> bachRecords)
    {
        final SomaticVariantFactory variantFactory = new SomaticVariantFactory();

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

            final CommonInfo variantInfo = variantContext.getCommonInfo();
            variantInfo.addFilter(bachRecord.filterType().toString());
            variantInfo.putAttribute(SNPEFF_IDENTIFIER, bachRecord.Annotations);

            // enrich with data from the ref genome
            if(mRefGenomeEnrichment.isValid())
            {
                final VariantType variantType = VariantType.type(variantContext);

                final String sequence = mRefGenomeEnrichment.getSequence(bachRecord.Chromosome, bachRecord.Position, bachRecord.Ref);

                final Optional<RepeatContext> repeatContext = mRefGenomeEnrichment.getRepeatContext(
                        variantType, sequence, bachRecord.Position, bachRecord.Alts);

                final String trinucleotideContext =
                        mRefGenomeEnrichment.getTrinucleotideContext(bachRecord.Chromosome, bachRecord.Position);

                final String microhomology = mRefGenomeEnrichment.getMicrohomology(
                        variantType, bachRecord.Ref, bachRecord.Alts, bachRecord.Position, sequence);

                variantInfo.putAttribute(TRINUCLEOTIDE_FLAG, trinucleotideContext);
                variantInfo.putAttribute(MICROHOMOLOGY_FLAG, microhomology);

                if(repeatContext.isPresent())
                {
                    variantInfo.putAttribute(REPEAT_COUNT_FLAG, repeatContext.get().count());
                    variantInfo.putAttribute(REPEAT_SEQUENCE_FLAG, repeatContext.get().sequence());
                }
            }

            bachRecord.setVariantContext(variantContext);

            SomaticVariant somVariant = variantFactory.createSomaticVariant(sampleId, variantContext);
            bachRecord.setSomaticVariant(somVariant);
        }
    }

    private void addPurpleData(final String sampleId, final List<BachelorGermlineVariant> bachRecords)
    {
        final PurityContext purityContext;
        final List<PurpleCopyNumber> copyNumbers;

        if (!mConfig.PurpleDataDir.isEmpty())
        {
            BACH_LOGGER.debug("sample({}) loading purple data from file using path {}", sampleId, mConfig.PurpleDataDir);

            try
            {
                purityContext = PurityContextFile.read(mConfig.PurpleDataDir, sampleId);
                copyNumbers = PurpleCopyNumberFile.read(PurpleCopyNumberFile.generateFilenameForReading(mConfig.PurpleDataDir, sampleId));
            }
            catch (IOException e)
            {
                BACH_LOGGER.error("failed to read purple data from {}: {}", mConfig.PurpleDataDir, e.toString());
                return;
            }
        }
        else
        {
            BACH_LOGGER.debug("sample({}) loading purple data from database", sampleId);

            purityContext = mDbAccess.readPurityContext(sampleId);

            if (purityContext == null)
            {
                BACH_LOGGER.warn("sample({}) failed to read purity data", sampleId);
                return;
            }

            copyNumbers = mDbAccess.readCopynumbers(sampleId);

            if(copyNumbers.isEmpty())
            {
                BACH_LOGGER.warn("sample({}) failed to read copy number data", sampleId);
                return;
            }
        }

        final PurityAdjuster purityAdjuster =
                new PurityAdjusterTypicalChromosome(purityContext.gender(), purityContext.bestFit().purity(), purityContext.bestFit().normFactor());

        final PurityAdjustedSomaticVariantFactory purityAdjustmentFactory =
                new PurityAdjustedSomaticVariantFactory(sampleId, purityAdjuster, copyNumbers);

        for (BachelorGermlineVariant bachRecord : bachRecords)
        {
            VariantContext purpleVarContext = purityAdjustmentFactory.enrich(bachRecord.getVariantContext());

            double adjustedCopyNumber = purpleVarContext.getAttributeAsDouble(PURPLE_CN_INFO, 0);
            double alleleFrequency = bachRecord.getSomaticVariant().alleleFrequency();
            double minorAlleleCopyNumber = purpleVarContext.getAttributeAsDouble(PURPLE_MINOR_ALLELE_CN_INFO, 0);

            double adjVaf = bachRecord.IsHomozygous ?
                    purityAdjuster.purityAdjustedVAFWithHomozygousNormal(bachRecord.Chromosome, adjustedCopyNumber, alleleFrequency)
                    : purityAdjuster.purityAdjustedVAFWithHeterozygousNormal(bachRecord.Chromosome, adjustedCopyNumber, alleleFrequency);

            if (Double.isNaN(adjVaf) || Double.isInfinite(adjVaf))
                adjVaf = 0;

            bachRecord.setEnrichmentData(adjVaf, minorAlleleCopyNumber, adjustedCopyNumber);
        }
    }

    private static final int HIGH_INDEL_REPEAT_COUNT = 8;

    private void applyIndelFilter(final List<BachelorGermlineVariant> bachRecords)
    {
        // currently the only filter is on INDELs with either microhomology matching the gain or loss, or with high repeat count
        final List<BachelorGermlineVariant> invalidRecords = bachRecords.stream().filter(x -> !x.isValid(false)).collect(Collectors.toList());
        invalidRecords.forEach(x -> bachRecords.remove(x));

        for(BachelorGermlineVariant bachRecord : bachRecords)
        {
            final SomaticVariant somaticVariant = bachRecord.getSomaticVariant();

            if(somaticVariant.type() == INDEL && (bachRecord.CodingEffect == SPLICE || bachRecord.CodingEffect == NONE))
            {
                int repeatCount = somaticVariant.repeatCount();

                if(repeatCount > HIGH_INDEL_REPEAT_COUNT)
                {
                    BACH_LOGGER.debug("Var({}) indel {} with high repeatCount({}) black-listed",
                            bachRecord.asString(), bachRecord.CodingEffect, repeatCount);

                    bachRecord.setPathogenicType(BLACK_LIST);
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
                    BACH_LOGGER.debug("Var({}) indel {} with ref, alt and microHom equal black-listed",
                            bachRecord.asString(), bachRecord.CodingEffect);
                    bachRecord.setPathogenicType(BLACK_LIST);
                }
            }
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

        boolean requirePurpleEnrichmentData = !mConfig.PurpleDataDir.isEmpty();

        for(final BachelorGermlineVariant bachRecord : bachRecords)
        {
            if (!bachRecord.isValid(requirePurpleEnrichmentData))
                continue;

            final SomaticVariant somaticVariant = bachRecord.getSomaticVariant();

            germlineVariants.add(ImmutableGermlineVariant.builder()
                    .chromosome(bachRecord.Chromosome)
                    .position(bachRecord.Position)
                    .filter(bachRecord.filterType())
                    .reported(bachRecord.isReportable())
                    .pathogenic(bachRecord.pathogenicType())
                    .type(somaticVariant.type().toString())
                    .ref(somaticVariant.ref())
                    .alts(somaticVariant.alt())
                    .gene(bachRecord.Gene)
                    .effects(bachRecord.Effects)
                    .codingEffect(bachRecord.CodingEffect)
                    .transcriptId(bachRecord.TranscriptId)
                    .alleleReadCount(bachRecord.getTumorAltCount())
                    .totalReadCount(bachRecord.getTumorReadDepth())
                    .adjustedCopyNumber(bachRecord.getAdjustedCopyNumber())
                    .adjustedVaf(bachRecord.getAdjustedVaf())
                    .trinucleotideContext(somaticVariant.trinucleotideContext())
                    .microhomology(somaticVariant.microhomology())
                    .repeatSequence(somaticVariant.repeatSequence())
                    .repeatCount(somaticVariant.repeatCount())
                    .hgvsProtein(bachRecord.HgvsProtein)
                    .hgvsCoding(bachRecord.HgvsCoding)
                    .biallelic(bachRecord.isBiallelic())
                    .minorAlleleJcn(bachRecord.getMinorAlleleCopyNumber())
                    .refStatus(bachRecord.refStatus())
                    .clinvarInfo(bachRecord.getClinvarConsolidatedInfo())
                    .build());
        }

        return germlineVariants;
    }

    private void writeToFile(final String sampleId, final List<GermlineVariant> germlineVariants)
    {
        try
        {
            final String filename = GermlineVariantFile.generateFilename(mConfig.OutputDir, sampleId);
            GermlineVariantFile.write(filename, germlineVariants);

            final String reportableVariantsFilename = ReportableGermlineVariantFile.generateFilename(mConfig.OutputDir, sampleId);

            final List<ReportableGermlineVariant> reportableVariants = Lists.newArrayList();

            for(final GermlineVariant germlineVariant : germlineVariants)
            {
                if(!germlineVariant.reported())
                    continue;

                reportableVariants.add(ImmutableReportableGermlineVariant.builder()
                        .gene(germlineVariant.gene())
                        .chromosome(germlineVariant.chromosome())
                        .position(germlineVariant.position())
                        .ref(germlineVariant.ref())
                        .alt(germlineVariant.alts())
                        .codingEffect(germlineVariant.codingEffect())
                        .hgvsCoding(germlineVariant.hgvsCoding())
                        .hgvsProtein(germlineVariant.hgvsProtein())
                        .alleleReadCount(germlineVariant.alleleReadCount())
                        .totalReadCount(germlineVariant.totalReadCount())
                        .adjustedVaf(germlineVariant.adjustedVaf())
                        .adjustedCopyNumber(germlineVariant.adjustedCopyNumber())
                        .biallelic(germlineVariant.biallelic())
                        .build());
            }

            ReportableGermlineVariantFile.write(reportableVariantsFilename, reportableVariants);
        }
        catch (final IOException e)
        {
            BACH_LOGGER.error("Error writing to outputFile: {}", e.toString());
        }
    }
}
