package com.hartwig.hmftools.purple.somatic;

import static java.lang.Math.abs;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.round;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.variant.CodingEffect.NONE;
import static com.hartwig.hmftools.common.variant.CodingEffect.UNDEFINED;
import static com.hartwig.hmftools.common.variant.VariantType.SNP;
import static com.hartwig.hmftools.common.variant.VariantVcfTags.GNOMAD_FREQ;
import static com.hartwig.hmftools.common.variant.VariantVcfTags.getGenotypeAttributeAsDouble;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.config.PurpleConstants.MB_PER_GENOME;
import static com.hartwig.hmftools.purple.config.PurpleConstants.TARGET_REGIONS_CN_DIFF;
import static com.hartwig.hmftools.purple.config.PurpleConstants.TARGET_REGIONS_CN_PERC_DIFF;
import static com.hartwig.hmftools.purple.config.PurpleConstants.TARGET_REGIONS_MAX_AF;
import static com.hartwig.hmftools.purple.config.PurpleConstants.TARGET_REGIONS_MAX_AF_DIFF;
import static com.hartwig.hmftools.purple.config.TargetRegionsData.TMB_GENE_EXCLUSIONS;

import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.purple.config.TargetRegionsData;

import htsjdk.variant.vcf.VCFConstants;

public class TumorMutationalLoad
{
    private final TargetRegionsData mTargetRegions;
    private double mLoad;
    private double mBurden;
    private int mUnclearBurdenVariants;
    private int mUnclearLoadVariants;

    public TumorMutationalLoad(final TargetRegionsData targetRegions)
    {
        mTargetRegions = targetRegions;
        mLoad = 0;
        mBurden = 0;
        mUnclearBurdenVariants = 0;
        mUnclearLoadVariants = 0;
    }

    public double load() { return mLoad; }
    public double burden() { return mBurden; }

    public int tml()
    {
        if(!mTargetRegions.hasTargetRegions())
            return (int)round(mLoad);
        else
            return mTargetRegions.calcTml(mBurden);
    }

    public double burdenPerMb() { return mBurden / MB_PER_GENOME; }

    public void processVariant(final SomaticVariant variant, double purity)
    {
        final VariantImpact variantImpact = variant.variantImpact();

        boolean isUnclearGermline = false;

        if(mTargetRegions.hasTargetRegions())
        {
            if(!mTargetRegions.inTargetRegions(variant.chromosome(), variant.position()))
                return;

            if(variant.isHotspot())
                return;

            if(variant.type() != SNP)
                return;

            if(variantImpact.WorstCodingEffect == NONE || variantImpact.WorstCodingEffect == UNDEFINED)
                return;

            if(TMB_GENE_EXCLUSIONS.contains(variantImpact.CanonicalGeneName))
                return;

            double gnomadFreq = variant.context().getAttributeAsDouble(GNOMAD_FREQ, 0);
            if(gnomadFreq > 0)
                return;

            double rawAf = getGenotypeAttributeAsDouble(variant.context().getGenotype(0), VCFConstants.ALLELE_FREQUENCY_KEY, 0);

            if(rawAf > TARGET_REGIONS_MAX_AF)
                return;

            // - VCN <= Major Allele CN + min(20%,0.5)
            double variantCn = variant.copyNumber(); // of the segment it's on
            double segmentCn = variant.decorator().adjustedCopyNumber(); // of the segment it's on
            double minorAlleleCn = variant.decorator().minorAlleleCopyNumber();
            double majorAlleleCn = segmentCn - minorAlleleCn;
            double diffThreshold = min(TARGET_REGIONS_CN_DIFF, majorAlleleCn * TARGET_REGIONS_CN_PERC_DIFF);
            if(variantCn > majorAlleleCn + diffThreshold)
                return;

            double refPurity = 1 - purity;

            double denom = 2 * refPurity + segmentCn * purity;
            double minorVAF = (refPurity + minorAlleleCn * purity) / denom;
            double majorVAF = (refPurity + majorAlleleCn * purity) / denom;

            isUnclearGermline = abs(majorVAF - rawAf) < TARGET_REGIONS_MAX_AF_DIFF || abs(minorVAF - rawAf) < TARGET_REGIONS_MAX_AF_DIFF;

            PPL_LOGGER.debug(format("var(%s) af(%.2f) copyNumber(vcn=%.2f segCn=%.2f majorCn=%.2f minorVaf=%.2f majorVaf=%.2f) status(%s) for target-regions TMB",
                    variant.toString(), rawAf, variantCn, segmentCn, majorAlleleCn, minorVAF, majorVAF,
                    isUnclearGermline ? "unclear" : "somatic"));
        }

        if(isUnclearGermline)
            ++mUnclearBurdenVariants;
        else
            mBurden++;

        if(variantImpact.WorstCodingEffect.equals(CodingEffect.MISSENSE))
        {
            if(isUnclearGermline)
                ++mUnclearLoadVariants;
            else
                mLoad++;
        }
    }

    public void calculateUnclearVariants()
    {
        if(!mTargetRegions.hasTargetRegions())
            return;

        mBurden += calcVariantsFromUnclear(mUnclearBurdenVariants);
        mLoad += calcVariantsFromUnclear(mUnclearLoadVariants);
    }

    private static int calcVariantsFromUnclear(int unclearCount)
    {
        return unclearCount > 0 ? (int)round(pow(unclearCount,2) / (double)(unclearCount + 8)) : 0;
    }
}
