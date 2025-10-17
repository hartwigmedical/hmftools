package com.hartwig.hmftools.purple.germline;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.purple.GermlineStatus.AMPLIFICATION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HET_DELETION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HOM_DELETION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.UNKNOWN;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.sv.SvUtils.SV_GERMLINE_AD_THRESHOLD;
import static com.hartwig.hmftools.common.sv.SvUtils.SV_GERMLINE_AF_THRESHOLD;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ALLELE_FRACTION;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH_PAIR;
import static com.hartwig.hmftools.common.sv.SvVcfTags.TOTAL_FRAGS;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_REV;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsDouble;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsInt;
import static com.hartwig.hmftools.purple.PurpleConstants.GERMLINE_SV_TINC_FACTOR;
import static com.hartwig.hmftools.purple.PurpleConstants.GERMLINE_SV_TINC_HOTSPOT_MULTIPLIER;
import static com.hartwig.hmftools.purple.PurpleConstants.GERMLINE_SV_TINC_MARGIN;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.PurpleConstants.WINDOW_SIZE;
import static com.hartwig.hmftools.purple.sv.SomaticSvCache.addEnrichedVariantContexts;

import java.util.Iterator;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.sv.ImmutableEnrichedStructuralVariant;
import com.hartwig.hmftools.common.sv.ImmutableEnrichedStructuralVariantLeg;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;
import com.hartwig.hmftools.common.variant.GenotypeIds;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.purple.PurpleConfig;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.purple.sv.SomaticSvCache;
import com.hartwig.hmftools.purple.sv.VariantContextCollection;

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
    private final VCFHeader mVcfHeader;
    private final VariantContextCollection mVariantCollection;
    private final GenotypeIds mGenotypeIds;

    private PurityContext mPurityContext;
    private List<ObservedRegion> mFittedRegions;
    private List<PurpleCopyNumber> mCopyNumbers;

    public GermlineSvCache(final String version, final String inputVcf, final PurpleConfig config)
    {
        if(!inputVcf.isEmpty())
        {
            VcfFileReader vcfReader = new VcfFileReader(inputVcf);
            mVcfHeader = SomaticSvCache.generateOutputHeader(version, vcfReader.vcfHeader());

            mGenotypeIds = GenotypeIds.fromVcfHeader(mVcfHeader, config.ReferenceId, config.TumorId);
            mVariantCollection = new VariantContextCollection(mVcfHeader);

            for(VariantContext context : vcfReader.iterator())
            {
                mVariantCollection.addVariant(context);
            }

            PPL_LOGGER.info("loaded {} germline SVs from {}", mVariantCollection.variants().size(), inputVcf);

            vcfReader.close();
        }
        else
        {
            mVcfHeader = null;
            mVariantCollection = new VariantContextCollection(null);
            mGenotypeIds = null;
        }

        mPurityContext = null;
        mFittedRegions = null;
        mCopyNumbers = null;
    }

    public GermlineSvCache()
    {
        this(null, "", null);
    }

    public List<StructuralVariant> germlineVariants() { return mVariantCollection.variants(); }

    public void checkTincVariants(final SomaticSvCache svCache, double tincLevel)
    {
        if(tincLevel == 0)
            return;

        // switch germline variants with tumor evidence across to the somatic cache

        // convert to SVs so paired breakends can be accessed together
        List<StructuralVariant> variants = Lists.newArrayList(mVariantCollection.variants());

        boolean hasTransfers = false;

        double tincAdjustment = GERMLINE_SV_TINC_FACTOR * tincLevel + GERMLINE_SV_TINC_MARGIN;

        List<StructuralVariant> transferredVariants = Lists.newArrayList();

        for(StructuralVariant variant : variants)
        {
            VariantContext variantContext = variant.startContext();

            if(variantContext.isFiltered())
                continue;

            Genotype refGenotype = variantContext.getGenotype(mGenotypeIds.ReferenceOrdinal);

            int refNonAltFrags = getGenotypeAttributeAsInt(refGenotype, REF_DEPTH, 0)
                    + getGenotypeAttributeAsInt(refGenotype, REF_DEPTH_PAIR, 0);

            int refAltFrags = getGenotypeAttributeAsInt(refGenotype, TOTAL_FRAGS, 0);

            double refDepth = refAltFrags + refNonAltFrags;
            double refAF = refAltFrags / refDepth;

            Genotype tumorGenotype = variantContext.getGenotype(mGenotypeIds.TumorOrdinal);
            int tumorAltFrags = getGenotypeAttributeAsInt(tumorGenotype, TOTAL_FRAGS, 0);
            double tumorAF = getOrCalculateAlleleFrequency(tumorGenotype);

            double tincRecovery = tumorAF * tincAdjustment;

            if(variant.hotspot())
                tincRecovery *= GERMLINE_SV_TINC_HOTSPOT_MULTIPLIER;

            double adjustedRefAF = refAF - tincRecovery;
            double adjustedRefAltFrags = (int)round(refAltFrags - refDepth * tincRecovery);

            // re-test the reduced ref values

            // a breakend is germline filtered if germlineAF/tumorAF >= 0.1 and germlineAD / tumorAD >= 0.01

            double adRatio = adjustedRefAltFrags / (double)max(tumorAltFrags, 1);

            boolean germlineFiltered = adjustedRefAF >= SV_GERMLINE_AF_THRESHOLD * tumorAF && adRatio >= SV_GERMLINE_AD_THRESHOLD;

            if(germlineFiltered)
                continue;

            hasTransfers = true;

            svCache.addTincVariant(variant);
            transferredVariants.add(variant);
        }

        if(hasTransfers)
        {
            mVariantCollection.clear();

            PPL_LOGGER.info("tincLevel({}) recovered {} germline SVs", tincLevel, transferredVariants.size());

            for(StructuralVariant variant : variants)
            {
                if(transferredVariants.contains(variant))
                    continue;

                mVariantCollection.addVariant(variant.startContext());

                if(variant.endContext() != null)
                    mVariantCollection.addVariant(variant.endContext());
            }
        }
    }

    public void annotateCopyNumberInfo(
            final List<ObservedRegion> fittedRegions, final List<PurpleCopyNumber> copyNumbers, final PurityContext purityContext)
    {
        mPurityContext = purityContext;
        mFittedRegions = fittedRegions;
        mCopyNumbers = copyNumbers;

        // note this routine clears and rebuilds the variant context collection, so a copy must be taken of the initial variants
        List<StructuralVariant> variants = Lists.newArrayList(mVariantCollection.variants());

        mVariantCollection.clear();

        for(StructuralVariant variant : variants)
        {
            RegionMatchInfo[] matchedFittedRegions = matchFittedRegions(variant);

            PurpleCopyNumber[] matchedCopyNumbers = matchCopyNumbers(variant);

            annotateVariant(variant, matchedFittedRegions, matchedCopyNumbers);
        }
    }

    public void write(final String outputVcf)
    {
        if(mVcfHeader == null || outputVcf.isEmpty())
            return;

        try
        {
            VariantContextWriter writer = new VariantContextWriterBuilder()
                    .setOutputFile(outputVcf)
                    .setReferenceDictionary(mVcfHeader.getSequenceDictionary())
                    .setIndexCreator(new TabixIndexCreator(mVcfHeader.getSequenceDictionary(), new TabixFormat()))
                    .setOutputFileType(VariantContextWriterBuilder.OutputType.BLOCK_COMPRESSED_VCF)
                    .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                    .build();

            writer.writeHeader(mVcfHeader);

            // now write the variants
            Iterator<VariantContext> variantIter = mVariantCollection.iterator();

            while(variantIter.hasNext())
            {
                VariantContext variant = variantIter.next();

                if(variant.isFiltered())
                    continue;

                writer.add(variant);
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
                byte impliedOrientation = refCnChange > 0 ? ORIENT_REV : ORIENT_FWD;

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
}
