package com.hartwig.hmftools.purple.germline;

import static com.hartwig.hmftools.common.purple.GermlineStatus.AMPLIFICATION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HET_DELETION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HOM_DELETION;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.ALLELE_FRACTION;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.REFERENCE_BREAKEND_READPAIR_COVERAGE;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.REFERENCE_BREAKEND_READ_COVERAGE;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.VARIANT_FRAGMENT_BREAKEND_COVERAGE;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.VARIANT_FRAGMENT_BREAKPOINT_COVERAGE;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsDouble;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsInt;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.config.PurpleConstants.WINDOW_SIZE;
import static com.hartwig.hmftools.purple.sv.SomaticSvCache.addEnrichedVariantContexts;

import java.io.File;
import java.util.Iterator;
import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.sv.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.sv.ImmutableEnrichedStructuralVariant;
import com.hartwig.hmftools.common.sv.ImmutableEnrichedStructuralVariantLeg;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantHeader;
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;
import com.hartwig.hmftools.purple.config.ReferenceData;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.purple.sv.StructuralRefContextEnrichment;
import com.hartwig.hmftools.purple.sv.VariantContextCollection;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

public class GermlineSvCache
{
    private final String mOutputVcfFilename;
    private final Optional<VCFHeader> mVcfHeader;
    private final VariantContextCollection mVariantCollection;
    private final IndexedFastaSequenceFile mRefGenomeFile;

    private final PurityContext mPurityContext;
    private final List<ObservedRegion> mFittedRegions;
    private final List<PurpleCopyNumber> mCopyNumbers;

    public GermlineSvCache()
    {
        mVcfHeader = Optional.empty();
        mOutputVcfFilename = Strings.EMPTY;
        mVariantCollection = new VariantContextCollection(null);
        mRefGenomeFile = null;
        mPurityContext = null;
        mFittedRegions = null;
        mCopyNumbers = null;
    }

    public GermlineSvCache(
            final String version, final String inputVcf, final String outputVcf, final ReferenceData referenceData,
            final List<ObservedRegion> fittedRegions, final List<PurpleCopyNumber> copyNumbers, final PurityContext purityContext)
    {
        mPurityContext = purityContext;
        mFittedRegions = fittedRegions;
        mCopyNumbers = copyNumbers;

        final VCFFileReader vcfReader = new VCFFileReader(new File(inputVcf), false);
        mOutputVcfFilename = outputVcf;
        mVcfHeader = Optional.of(generateOutputHeader(version, vcfReader.getFileHeader()));
        mVariantCollection = new VariantContextCollection(mVcfHeader.get());
        mRefGenomeFile = referenceData.RefGenome;

        for(VariantContext context : vcfReader)
        {
            mVariantCollection.add(context);
        }

        PPL_LOGGER.info("loaded {} germline SVs from {}", mVariantCollection.variants().size(), inputVcf);

        vcfReader.close();
    }

    public List<StructuralVariant> variants() { return mVariantCollection.variants(); }

    public void write()
    {
        if(!mVcfHeader.isPresent() || mOutputVcfFilename.isEmpty())
            return;

        try
        {
            final VariantContextWriter writer = new VariantContextWriterBuilder()
                    .setOutputFile(mOutputVcfFilename)
                    .setReferenceDictionary(mVcfHeader.get().getSequenceDictionary())
                    .setIndexCreator(new TabixIndexCreator(mVcfHeader.get().getSequenceDictionary(), new TabixFormat()))
                    .setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF)
                    .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                    .build();

            final StructuralRefContextEnrichment refEnricher = new StructuralRefContextEnrichment(mRefGenomeFile, writer::add);

            writer.writeHeader(refEnricher.enrichHeader(mVcfHeader.get()));

            // may be no reason to use the enriched collection, unsure if it adds any value
            for(StructuralVariant variant : mVariantCollection.variants())
            {
                annotateVariant(variant);
            }

            // now write the variants
            Iterator<VariantContext> variantIter = mVariantCollection.iterator();

            while(variantIter.hasNext())
            {
                VariantContext variant = variantIter.next();

                if(variant.isFiltered())
                    continue;

                refEnricher.accept(variant);
            }

            refEnricher.flush();

            writer.close();
        }
        catch(Exception e)
        {
            PPL_LOGGER.error("failed to write germline SV VCF({}): {}", mOutputVcfFilename, e.toString());
        }
    }

    private void annotateVariant(
            final StructuralVariant variant, final ObservedRegion[] fittedRegions, final PurpleCopyNumber[] copyNumbers)
    {
        ImmutableEnrichedStructuralVariant.Builder builder = ImmutableEnrichedStructuralVariant.builder().from(variant);

        double junctionCopyNumberTotal = 0;
        int legCount = 0;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            StructuralVariantLeg leg = se == SE_START ? variant.start() : variant.end();

            if(leg == null)
                continue;

            VariantContext context = se == SE_START ? variant.startContext() : variant.endContext();

            ImmutableEnrichedStructuralVariantLeg.Builder legBuilder = ImmutableEnrichedStructuralVariantLeg.builder().from(leg);

            ObservedRegion fittedRegion = fittedRegions[se];
            PurpleCopyNumber copyNumber = copyNumbers[se];

            // initialise all values
            legBuilder.adjustedAlleleFrequency(0.0);
            legBuilder.adjustedCopyNumber(0.0);
            legBuilder.adjustedCopyNumberChange(0.0);

            if(copyNumber != null)
            {
                double adjustedCN = copyNumber.averageTumorCopyNumber();
                legBuilder.adjustedCopyNumber(copyNumber.averageTumorCopyNumber());

                if(fittedRegion != null)
                {
                    double cnChange = adjustedCN - fittedRegion.refNormalisedCopyNumber();

                    legBuilder.adjustedCopyNumberChange(cnChange);
                }

                double purity = mPurityContext.bestFit().purity();

                if(context.getGenotypes().size() > 1 && purity > 0)
                {
                    // [tumorAF*[2*(1-purity)+adjCNStart*purity] -refAF*2*(1-purity)]/adjCNStart/purity = adjAFStart
                    // rearranged from tumorAF = [refAF*2*(1-purity) + adjAFStart*purity*adjCNStart] / [2*(1-purity) + purity*adjCNStart]
                    Genotype refGenotype = context.getGenotype(0);
                    Genotype tumorGenotype = context.getGenotype(1);

                    double refAF = getOrCalculateAlleleFrequency(refGenotype);
                    double tumorAF = getOrCalculateAlleleFrequency(tumorGenotype);

                    double refPurity = 2 * (1 - purity);
                    double adjustedAF = (tumorAF * (refPurity + adjustedCN * purity) - refAF * refPurity) / (adjustedCN * purity);

                    legBuilder.adjustedAlleleFrequency(adjustedAF);

                    ++legCount;
                    junctionCopyNumberTotal += adjustedAF * adjustedCN;
                }
            }

            if(se == SE_START)
                builder.start(legBuilder.build());
            else
                builder.end(legBuilder.build());
        }

        double junctionCopyNumber = legCount > 0 ? junctionCopyNumberTotal / legCount : 0;
        builder.junctionCopyNumber(junctionCopyNumber);

        addEnrichedVariantContexts(mVariantCollection, builder.build());
    }

    private double getOrCalculateAlleleFrequency(final Genotype genotype)
    {
        if(genotype.hasExtendedAttribute(ALLELE_FRACTION))
            return getGenotypeAttributeAsDouble(genotype, ALLELE_FRACTION, 0);

        int totalReadCoverage = getGenotypeAttributeAsInt(genotype, REFERENCE_BREAKEND_READ_COVERAGE, 0)
                + getGenotypeAttributeAsInt(genotype, REFERENCE_BREAKEND_READPAIR_COVERAGE, 0);

        int variantFrags = getGenotypeAttributeAsInt(genotype, VARIANT_FRAGMENT_BREAKPOINT_COVERAGE, 0) +
                getGenotypeAttributeAsInt(genotype, VARIANT_FRAGMENT_BREAKEND_COVERAGE, 0);

        double total = variantFrags + totalReadCoverage;
        return variantFrags / total;
    }

    private void annotateVariant(final StructuralVariant variant)
    {
        // not expecting too may germline SVs so just search the full collections for each one
        int copyNumberIndex = 0;
        int regionIndex = 0;

        ObservedRegion[] matchedFittedRegions = new ObservedRegion[SE_PAIR];
        PurpleCopyNumber[] matchedCopyNumbers = new PurpleCopyNumber[SE_PAIR];

        for(int se = SE_START; se <= SE_END; ++se)
        {
            StructuralVariantLeg leg = se == SE_START ? variant.start() : variant.end();

            if(leg == null)
                continue;

            if(matchedCopyNumbers[SE_START] != null && breakendMatchesCopyNumber(leg, matchedCopyNumbers[SE_START]))
            {
                // other breakend is within the same CN segment
                matchedCopyNumbers[SE_END] = matchedCopyNumbers[SE_START];
            }
            else
            {
                for(; copyNumberIndex < mCopyNumbers.size(); ++copyNumberIndex)
                {
                    PurpleCopyNumber copyNumber = mCopyNumbers.get(copyNumberIndex);

                    if(breakendMatchesCopyNumber(leg, copyNumber))
                    {
                        matchedCopyNumbers[se] = copyNumber;
                        break;
                    }
                }
            }

            for(; regionIndex < mFittedRegions.size(); ++regionIndex)
            {
                ObservedRegion region = mFittedRegions.get(regionIndex);

                if(!isValidRegion(region))
                    continue;

                ObservedRegion nextRegion = regionIndex < mFittedRegions.size() - 1 ? mFittedRegions.get(regionIndex + 1) : null;
                if(nextRegion != null && !nextRegion.chromosome().equals(region.chromosome()))
                    nextRegion = null;

                ObservedRegion breakendRegion = breakendMatchesRegion(leg, region, nextRegion);

                if(breakendRegion != null)
                {
                    matchedFittedRegions[se] = breakendRegion;
                    break;
                }
            }
        }

        annotateVariant(variant, matchedFittedRegions, matchedCopyNumbers);
    }

    private static boolean breakendMatchesCopyNumber(final StructuralVariantLeg leg, final PurpleCopyNumber copyNumber)
    {
        return copyNumber.chromosome().equals(leg.chromosome()) && positionWithin(leg.position(), copyNumber.start(), copyNumber.end());
    }

    private static boolean isValidRegion(final ObservedRegion region)
    {
        return region.germlineStatus() == HOM_DELETION || region.germlineStatus() == HET_DELETION || region.germlineStatus() == AMPLIFICATION;
    }

    private ObservedRegion breakendMatchesRegion(final StructuralVariantLeg leg, final ObservedRegion region, final ObservedRegion nextRegion)
    {
        // try to match on the start region or the end region, matching on orientation and within max 1 depth window

        for(int se = SE_START; se <= SE_END; ++se)
        {
            int regionStart;
            int regionEnd;
            byte regionOrientation;

            if(se == SE_START)
            {
                regionStart = region.minStart();
                regionEnd = region.maxStart();
                regionOrientation = region.germlineStatus() == AMPLIFICATION ? NEG_ORIENT : POS_ORIENT;
            }
            else
            {
                regionOrientation = region.germlineStatus() == AMPLIFICATION ? POS_ORIENT : NEG_ORIENT;

                if(nextRegion != null)
                {
                    regionStart = nextRegion.minStart() - 1;
                    regionEnd = nextRegion.maxStart() - 1;
                }
                else
                {
                    regionStart = region.end();
                    regionEnd = region.end();
                }
            }

            if(leg.orientation() != regionOrientation)
                continue;

            regionStart -= WINDOW_SIZE;
            regionEnd += WINDOW_SIZE;

            if(positionWithin(leg.position(), regionStart, regionEnd))
            {
                return region;
            }
        }

        return null;
    }

    private static VCFHeader generateOutputHeader(final String purpleVersion, final VCFHeader template)
    {
        return StructuralVariantHeader.generateHeader(purpleVersion, template);
    }
}
