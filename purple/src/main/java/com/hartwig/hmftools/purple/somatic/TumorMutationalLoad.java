package com.hartwig.hmftools.purple.somatic;

import static java.lang.Math.pow;

import static com.hartwig.hmftools.common.variant.CodingEffect.NONE;
import static com.hartwig.hmftools.common.variant.CodingEffect.UNDEFINED;
import static com.hartwig.hmftools.common.variant.SomaticLikelihood.HIGH;
import static com.hartwig.hmftools.common.variant.SomaticLikelihood.LOW;
import static com.hartwig.hmftools.common.variant.SomaticLikelihood.MEDIUM;
import static com.hartwig.hmftools.common.variant.PaveVcfTags.GNOMAD_FREQ;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PANEL_SOMATIC_LIKELIHOOD;
import static com.hartwig.hmftools.common.variant.VariantType.SNP;
import static com.hartwig.hmftools.purple.PurpleConstants.CODING_BASES_PER_GENOME;
import static com.hartwig.hmftools.purple.PurpleConstants.TARGETED_TMB_GENE_EXCLUSIONS;
import static com.hartwig.hmftools.purple.PurpleConstants.TUMOR_MUT_LOAD_MIN_VAF;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.SomaticLikelihood;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.purple.targeted.TargetRegionsData;

public class TumorMutationalLoad
{
    private final TargetRegionsData mTargetRegions;
    private final boolean mTumorOnly;
    private int mLoad;
    private int mBurden;
    private int mUnclearVariants;

    public TumorMutationalLoad(final TargetRegionsData targetRegions, boolean tumorOnly)
    {
        mTargetRegions = targetRegions;
        mTumorOnly = tumorOnly;

        mLoad = 0;
        mBurden = 0;
        mUnclearVariants = 0;
    }

    public double burden() { return mBurden; }

    public double calcTml()
    {
        if(!mTargetRegions.hasTargetRegions())
            return mLoad;

        if(mTargetRegions.codingBases() == 0)
            return 0;

        double unclearLoad = 0;

        if(mUnclearVariants > 0)
        {
            double unclearFactor = mTargetRegions.codingBases() / mTargetRegions.codingBaseFactor();
            unclearLoad = pow(mUnclearVariants, 2) / (mUnclearVariants + unclearFactor);
        }

        double combinedLoad = mLoad + unclearLoad;
        double calcTml = combinedLoad * mTargetRegions.tmlRatio() * CODING_BASES_PER_GENOME / mTargetRegions.codingBases();

        PPL_LOGGER.debug(String.format("calculated tml(%.4f) variants(high=%d unclear=%d adjusted=%.1f)",
                calcTml, mLoad, mUnclearVariants, unclearLoad));

        return calcTml;
    }

    public void processVariant(final SomaticVariant variant)
    {
        if(variant.alleleFrequency() < TUMOR_MUT_LOAD_MIN_VAF)
            return;

        if(mTargetRegions.hasTargetRegions())
        {
            processTargetedRegionVariant(variant);
            return;
        }

        ++mBurden;

        if(variant.variantImpact().WorstCodingEffect.equals(CodingEffect.MISSENSE))
            mLoad++;
    }

    private void processTargetedRegionVariant(final SomaticVariant variant)
    {
        // test criteria to count a variant towards TMB
        final VariantImpact variantImpact = variant.variantImpact();

        if(!mTargetRegions.inTargetRegions(variant.chromosome(), variant.position()))
            return;

        if(variant.isHotspot())
            return;

        if(variant.type() != SNP)
            return;

        if(!variantImpact.WorstCodingEffect.equals(CodingEffect.MISSENSE))
            return;

        if(TARGETED_TMB_GENE_EXCLUSIONS.contains(variantImpact.GeneName))
            return;

        double gnomadFreq = variant.context().getAttributeAsDouble(GNOMAD_FREQ, 0);
        if(gnomadFreq > 0)
            return;

        boolean isUnclear = false;

        if(mTumorOnly)
        {
            SomaticLikelihood somaticLikelihood = variant.context().hasAttribute(PANEL_SOMATIC_LIKELIHOOD) ?
                    SomaticLikelihood.valueOf(variant.context().getAttributeAsString(PANEL_SOMATIC_LIKELIHOOD, "")) : LOW;

            if(somaticLikelihood != HIGH)
            {
                if(somaticLikelihood == MEDIUM)
                    isUnclear = true;
                else
                    return;
            }
        }

        if(isUnclear)
            ++mUnclearVariants;
        else
            ++mLoad;

        PPL_LOGGER.trace("variant({}) tml somatic status({})", variant, isUnclear ? "unclear" : "high");
    }
}
