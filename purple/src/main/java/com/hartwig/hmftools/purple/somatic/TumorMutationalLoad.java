package com.hartwig.hmftools.purple.somatic;

import static java.lang.Math.pow;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeConstants.BASES_PER_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeConstants.CODING_BASES_PER_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeConstants.MB_PER_GENOME;
import static com.hartwig.hmftools.common.variant.SomaticLikelihood.HIGH;
import static com.hartwig.hmftools.common.variant.SomaticLikelihood.LOW;
import static com.hartwig.hmftools.common.variant.SomaticLikelihood.MEDIUM;
import static com.hartwig.hmftools.common.variant.PaveVcfTags.GNOMAD_FREQ;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PANEL_SOMATIC_LIKELIHOOD;
import static com.hartwig.hmftools.common.variant.VariantType.SNP;
import static com.hartwig.hmftools.purple.PurpleConstants.MUTATIONAL_LOAD_CODING_BASE_FACTOR;
import static com.hartwig.hmftools.purple.PurpleConstants.TMB_GENE_EXCLUSIONS;
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

    private final long mCodingBases;

    private int mLoad;
    private int mBurden;
    private int mUnclearMissenseVariants;
    private int mUnclearVariants;

    public TumorMutationalLoad(final TargetRegionsData targetRegions, boolean tumorOnly)
    {
        mTargetRegions = targetRegions;
        mTumorOnly = tumorOnly;

        mCodingBases = mTargetRegions.hasTargetRegions() ? mTargetRegions.codingBases() : CODING_BASES_PER_GENOME;

        mLoad = 0;
        mBurden = 0;
        mUnclearVariants = 0;
        mUnclearMissenseVariants = 0;
    }

    public double burden() { return mBurden; }

    public double calcTml()
    {
        if(!mTargetRegions.hasTargetRegions() && !mTumorOnly)
            return mLoad;

        if(mCodingBases == 0)
            return 0;

        double unclearCount = 0;

        if(mUnclearMissenseVariants > 0)
        {
            double unclearFactor = mCodingBases / MUTATIONAL_LOAD_CODING_BASE_FACTOR;
            unclearCount = pow(mUnclearMissenseVariants, 2) / (mUnclearMissenseVariants + unclearFactor);
        }

        double calcTml = mLoad + unclearCount;

        if(mTargetRegions.hasTargetRegions())
        {
            calcTml *= mTargetRegions.tmlRatio() * CODING_BASES_PER_GENOME / mCodingBases;
        }

        PPL_LOGGER.debug(String.format("calculated tml(%.4f) variants(high=%d unclear=%d adjusted=%.1f)",
                calcTml, mLoad, mUnclearMissenseVariants, unclearCount));

        return calcTml;
    }

    public double calcWgsTmb()
    {
        if(mTargetRegions.hasTargetRegions())
        {
            PPL_LOGGER.warn("invalid TMB call in targeted mode");
            return -1;
        }

        if(!mTumorOnly)
        {
            return mBurden / MB_PER_GENOME;
        }

        double unclearCount = 0;

        if(mUnclearVariants > 0)
        {
            double unclearFactor = BASES_PER_GENOME / MUTATIONAL_LOAD_CODING_BASE_FACTOR;
            unclearCount = pow(mUnclearVariants, 2) / (mUnclearVariants + unclearFactor);
        }

        double calcTmb = mBurden + unclearCount;

        PPL_LOGGER.debug(String.format("calculated tmb(%.4f) variants(high=%d unclear=%d adjusted=%.1f)",
                calcTmb, mBurden, mUnclearVariants, unclearCount));

        return calcTmb / MB_PER_GENOME;
    }

    public void processVariant(final SomaticVariant variant)
    {
        if(variant.alleleFrequency() < TUMOR_MUT_LOAD_MIN_VAF)
            return;

        if(mTargetRegions.hasTargetRegions() || mTumorOnly)
        {
            processTargetedRegionOrTumorOnlyVariant(variant);
            return;
        }

        ++mBurden;

        if(isMissense(variant.variantImpact()))
            mLoad++;
    }

    private static boolean isMissense(final VariantImpact variantImpact)
    {
        return variantImpact.WorstCodingEffect.equals(CodingEffect.MISSENSE);
    }

    private void processTargetedRegionOrTumorOnlyVariant(final SomaticVariant variant)
    {
        // test criteria to count a variant towards TMB
        VariantImpact variantImpact = variant.variantImpact();

        if(mTargetRegions.hasTargetRegions() && !mTargetRegions.inTargetRegions(variant.chromosome(), variant.position()))
            return;

        if(variant.isHotspot())
            return;

        if(variant.type() != SNP)
            return;

        double gnomadFreq = variant.context().getAttributeAsDouble(GNOMAD_FREQ, 0);
        if(gnomadFreq > 0)
            return;

        if(TMB_GENE_EXCLUSIONS.contains(variantImpact.GeneName))
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
            ++mBurden;

        if(isMissense(variantImpact))
        {
            if(isUnclear)
                ++mUnclearMissenseVariants;
            else
                ++mLoad;
        }

        // PPL_LOGGER.trace("variant({}) mutational calc status({})", variant, isUnclear ? "unclear" : "high");
    }
}
