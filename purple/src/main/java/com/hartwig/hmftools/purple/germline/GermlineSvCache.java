package com.hartwig.hmftools.purple.germline;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.purple.GermlineStatus.AMPLIFICATION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HET_DELETION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HOM_DELETION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.UNKNOWN;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ALLELE_FRACTION;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH_PAIR;
import static com.hartwig.hmftools.common.sv.SvVcfTags.TOTAL_FRAGS;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsDouble;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsInt;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.PurpleConstants.WINDOW_SIZE;
import static com.hartwig.hmftools.purple.sv.SomaticSvCache.addEnrichedVariantContexts;

import java.util.Iterator;
import java.util.List;
import java.util.Optional;

import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.sv.ImmutableEnrichedStructuralVariant;
import com.hartwig.hmftools.common.sv.ImmutableEnrichedStructuralVariantLeg;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantHeader;
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;
import com.hartwig.hmftools.common.variant.GenotypeIds;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.purple.PurpleConfig;
import com.hartwig.hmftools.purple.ReferenceData;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.purple.sv.StructuralRefContextEnrichment;
import com.hartwig.hmftools.purple.sv.VariantContextCollection;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;

public class GermlineSvCache
{
    private final Optional<VCFHeader> mVcfHeader;
    private final VariantContextCollection mVariantCollection;
    private final IndexedFastaSequenceFile mRefGenomeFile;
    private final GenotypeIds mGenotypeIds;

    private final PurityContext mPurityContext;
    private final List<ObservedRegion> mFittedRegions;
    private final List<PurpleCopyNumber> mCopyNumbers;

    public GermlineSvCache()
    {
        mVcfHeader = Optional.empty();
        mVariantCollection = new VariantContextCollection(null, false);
        mRefGenomeFile = null;
        mPurityContext = null;
        mFittedRegions = null;
        mCopyNumbers = null;
        mGenotypeIds = null;
    }

    public GermlineSvCache(
            final String version, final String inputVcf, final ReferenceData referenceData, final PurpleConfig config,
            final List<ObservedRegion> fittedRegions, final List<PurpleCopyNumber> copyNumbers, final PurityContext purityContext)
    {
        mPurityContext = purityContext;
        mFittedRegions = fittedRegions;
        mCopyNumbers = copyNumbers;

        VcfFileReader vcfReader = new VcfFileReader(inputVcf);
        mVcfHeader = Optional.of(generateOutputHeader(version, vcfReader.vcfHeader()));

        mGenotypeIds = mVcfHeader.isPresent() ? GenotypeIds.fromVcfHeader(mVcfHeader.get(), config.ReferenceId, config.TumorId) : null;

        mVariantCollection = new VariantContextCollection(mVcfHeader.get(), config.UseGridssSVs);
        mRefGenomeFile = referenceData.RefGenome;

        for(VariantContext context : vcfReader.iterator())
        {
            mVariantCollection.add(context);
        }

        PPL_LOGGER.info("loaded {} germline SVs from {}", mVariantCollection.variants().size(), inputVcf);

        vcfReader.close();
    }

    public List<StructuralVariant> variants() { return mVariantCollection.variants(); }

    public void write(final String outputVcf)
    {
        if(!mVcfHeader.isPresent() || outputVcf.isEmpty())
            return;

        try
        {
            final VariantContextWriter writer = new VariantContextWriterBuilder()
                    .setOutputFile(outputVcf)
                    .setReferenceDictionary(mVcfHeader.get().getSequenceDictionary())
                    .setIndexCreator(new TabixIndexCreator(mVcfHeader.get().getSequenceDictionary(), new TabixFormat()))
                    .setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF)
                    .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                    .build();

            final StructuralRefContextEnrichment refEnricher = new StructuralRefContextEnrichment(mRefGenomeFile, writer::add);

            VCFHeader header = mVcfHeader.get();
            refEnricher.enrichHeader(header);
            writer.writeHeader(header);

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

            writer.close();
        }
        catch(Exception e)
        {
            PPL_LOGGER.error("failed to write germline SV VCF({}): {}", outputVcf, e.toString());
            e.printStackTrace();
            System.exit(1);
        }
    }

    private void annotateVariant(final StructuralVariant variant)
    {
        RegionMatchInfo[] fittedRegions = matchFittedRegions(variant);

        PurpleCopyNumber[] copyNumbers = matchCopyNumbers(variant);

        annotateVariant(variant, fittedRegions, copyNumbers);
    }

    private void annotateVariant(
            final StructuralVariant variant, final RegionMatchInfo[] fittedRegions, final PurpleCopyNumber[] copyNumbers)
    {
        ImmutableEnrichedStructuralVariant.Builder builder = ImmutableEnrichedStructuralVariant.builder().from(variant);

        double junctionCopyNumberTotal = 0;
        int legCount = 0;

        ImmutableEnrichedStructuralVariantLeg.Builder legBuilderStart = null;
        ImmutableEnrichedStructuralVariantLeg.Builder legBuilderEnd = null;
        boolean[] adjustedCnChangeSet = new boolean[] {false, false};

        for(int se = SE_START; se <= SE_END; ++se)
        {
            StructuralVariantLeg leg = se == SE_START ? variant.start() : variant.end();

            if(leg == null)
                continue;

            VariantContext context = se == SE_START ? variant.startContext() : variant.endContext();

            ImmutableEnrichedStructuralVariantLeg.Builder legBuilder = ImmutableEnrichedStructuralVariantLeg.builder().from(leg);

            if(se == SE_START)
                legBuilderStart = legBuilder;
            else
                legBuilderEnd = legBuilder;

            ObservedRegion fittedRegion = fittedRegions[se] != null ? fittedRegions[se].Region : null;
            PurpleCopyNumber copyNumber = copyNumbers[se];

            // initialise all values
            legBuilder.adjustedAlleleFrequency(0.0);
            legBuilder.adjustedCopyNumber(0.0);
            legBuilder.adjustedCopyNumberChange(0.0);

            if(copyNumber != null)
            {
                double adjustedCN = copyNumber.averageTumorCopyNumber();

                if(fittedRegion != null)
                {
                    double cnChange = abs(fittedRegions[se].CopyNumberChange);

                    legBuilder.adjustedCopyNumberChange(cnChange);
                    adjustedCnChangeSet[se] = true;

                    // take the higher of the ref CN regions
                    if(fittedRegion.germlineStatus() == AMPLIFICATION)
                        adjustedCN += cnChange;
                }

                legBuilder.adjustedCopyNumber(adjustedCN);

                double purity = mPurityContext.bestFit().purity();

                if(mGenotypeIds.hasReference() && mGenotypeIds.hasTumor() && purity > 0)
                {
                    // [tumorAF*[2*(1-purity)+adjCNStart*purity] -refAF*2*(1-purity)]/adjCNStart/purity = adjAFStart
                    // rearranged from tumorAF = [refAF*2*(1-purity) + adjAFStart*purity*adjCNStart] / [2*(1-purity) + purity*adjCNStart]

                    // a check has already been made that the germline VCF has the tumor ID in genotype 0 and the ref in genotype 1
                    Genotype refGenotype = context.getGenotype(mGenotypeIds.ReferenceOrdinal);
                    Genotype tumorGenotype = context.getGenotype(mGenotypeIds.TumorOrdinal);

                    double refAF = getOrCalculateAlleleFrequency(refGenotype);
                    double tumorAF = getOrCalculateAlleleFrequency(tumorGenotype);

                    double refPurity = 2 * (1 - purity);
                    double adjustedAF = (tumorAF * (refPurity + adjustedCN * purity) - refAF * refPurity) / (adjustedCN * purity);

                    legBuilder.adjustedAlleleFrequency(adjustedAF);

                    ++legCount;
                    junctionCopyNumberTotal += adjustedAF * adjustedCN;
                }
            }
        }


        double junctionCopyNumber = legCount > 0 ? junctionCopyNumberTotal / legCount : 0;
        builder.junctionCopyNumber(junctionCopyNumber);

        // revert to using simple JCN where the variant didn't match a fitted region (eg from being too short)
        if(!adjustedCnChangeSet[SE_START])
            legBuilderStart.adjustedCopyNumberChange(junctionCopyNumber);

        builder.start(legBuilderStart.build());

        if(legBuilderEnd != null)
        {
            if(!adjustedCnChangeSet[SE_END])
                legBuilderEnd.adjustedCopyNumberChange(junctionCopyNumber);

            builder.end(legBuilderEnd.build());
        }

        addEnrichedVariantContexts(mVariantCollection, builder.build());
    }

    private double getOrCalculateAlleleFrequency(final Genotype genotype)
    {
        if(genotype.hasExtendedAttribute(ALLELE_FRACTION))
            return getGenotypeAttributeAsDouble(genotype, ALLELE_FRACTION, 0);

        int totalReadCoverage = getGenotypeAttributeAsInt(genotype, REF_DEPTH, 0)
                + getGenotypeAttributeAsInt(genotype, REF_DEPTH_PAIR, 0);

        int variantFrags = getGenotypeAttributeAsInt(genotype, TOTAL_FRAGS, 0);

        double total = variantFrags + totalReadCoverage;
        return variantFrags / total;
    }

    private class RegionMatchInfo
    {
        public final ObservedRegion Region;
        public final double CopyNumberChange;

        public RegionMatchInfo(final ObservedRegion region, final double copyNumberChange)
        {
            Region = region;
            CopyNumberChange = copyNumberChange;
        }
    }

    private static boolean isShortLocalVariant(final StructuralVariant variant)
    {
        if(variant.type() == DEL || variant.type() == DUP || variant.type() == INV || variant.type() == INS)
        {
            int length = variant.position(false) - variant.position(true);
            return length <= WINDOW_SIZE * 1.5;
        }

        return false;
    }

    private RegionMatchInfo[] matchFittedRegions(final StructuralVariant variant)
    {
        RegionMatchInfo[] matchedFittedRegions = new RegionMatchInfo[SE_PAIR];

        if(isShortLocalVariant(variant))
            return matchedFittedRegions;

        int regionIndex = 0;
        boolean regionChrMatched = false;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            StructuralVariantLeg leg = se == SE_START ? variant.start() : variant.end();

            if(leg == null)
                continue;

            if(se == SE_END && !leg.chromosome().equals(variant.chromosome(true)))
                regionChrMatched = false;

            for(; regionIndex < mFittedRegions.size(); ++regionIndex)
            {
                ObservedRegion region = mFittedRegions.get(regionIndex);

                if(!regionChrMatched && region.chromosome().equals(leg.chromosome()))
                    regionChrMatched = true;

                if(regionChrMatched && region.start() > leg.position() + WINDOW_SIZE) // stop searching and keep index in place
                    break;

                if(regionIndex == 0)
                    continue;

                if(region.germlineStatus() == UNKNOWN)
                    continue;

                if(!breakendMatchesRegion(leg, region))
                    continue;

                // find the preceding region (again not UNKNOWN)
                int prevIndex = regionIndex - 1;
                ObservedRegion prevRegion = mFittedRegions.get(prevIndex);

                while(prevRegion.germlineStatus() == UNKNOWN)
                {
                    --prevIndex;

                    if(prevIndex < 0)
                    {
                        prevRegion = null;
                        break;
                    }

                    prevRegion = mFittedRegions.get(prevIndex);
                }

                if(prevRegion == null)
                    continue;

                double refCnChange = region.refNormalisedCopyNumber() - prevRegion.refNormalisedCopyNumber();
                byte impliedOrientation = refCnChange > 0 ? NEG_ORIENT : POS_ORIENT;

                if(impliedOrientation != leg.orientation())
                    continue;

                ObservedRegion germlineRegion = isValidRegion(region) ? region : prevRegion;

                matchedFittedRegions[se] = new RegionMatchInfo(germlineRegion, refCnChange);
                break;
            }
        }

        return matchedFittedRegions;
    }

    private static boolean isValidRegion(final ObservedRegion region)
    {
        return region.germlineStatus() == HOM_DELETION || region.germlineStatus() == HET_DELETION || region.germlineStatus() == AMPLIFICATION;
    }

    private boolean breakendMatchesRegion(final StructuralVariantLeg leg, final ObservedRegion region)
    {
        if(!leg.chromosome().equals(region.chromosome()))
            return false;

        // try to match on the start region within max 1 depth window
        int regionStart = region.minStart() - WINDOW_SIZE;
        int regionEnd = region.maxStart() + WINDOW_SIZE;

        return positionWithin(leg.position(), regionStart, regionEnd);
    }

    private PurpleCopyNumber[] matchCopyNumbers(final StructuralVariant variant)
    {
        // not expecting too may germline SVs so just search the full collections for each one
        int copyNumberIndex = 0;
        boolean copyNumberChrMatched = false;

        PurpleCopyNumber[] matchedCopyNumbers = new PurpleCopyNumber[SE_PAIR];

        for(int se = SE_START; se <= SE_END; ++se)
        {
            StructuralVariantLeg leg = se == SE_START ? variant.start() : variant.end();

            if(leg == null)
                continue;

            if(se == SE_END && !leg.chromosome().equals(variant.chromosome(true)))
                copyNumberChrMatched = false;

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

                    if(!copyNumberChrMatched && copyNumber.chromosome().equals(leg.chromosome()))
                        copyNumberChrMatched = true;

                    if(breakendMatchesCopyNumber(leg, copyNumber))
                    {
                        matchedCopyNumbers[se] = copyNumber;
                        break;
                    }

                    if(copyNumberChrMatched && copyNumber.start() > leg.position())
                        break;
                }
            }
        }

        return matchedCopyNumbers;
    }

    private static boolean breakendMatchesCopyNumber(final StructuralVariantLeg leg, final PurpleCopyNumber copyNumber)
    {
        return copyNumber.chromosome().equals(leg.chromosome()) && positionWithin(leg.position(), copyNumber.start(), copyNumber.end());
    }

    private static VCFHeader generateOutputHeader(final String purpleVersion, final VCFHeader template)
    {
        return StructuralVariantHeader.generateHeader(purpleVersion, template);
    }
}
