package com.hartwig.hmftools.purple.somatic;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.variant.CodingEffect.NONE;
import static com.hartwig.hmftools.common.variant.CodingEffect.UNDEFINED;
import static com.hartwig.hmftools.common.variant.VariantType.SNP;
import static com.hartwig.hmftools.common.variant.VariantVcfTags.GNOMAD_FREQ;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.config.PurpleConstants.MB_PER_GENOME;
import static com.hartwig.hmftools.purple.config.PurpleConstants.TARGET_REGIONS_CN_DIFF;
import static com.hartwig.hmftools.purple.config.PurpleConstants.TARGET_REGIONS_CN_PERC_DIFF;
import static com.hartwig.hmftools.purple.config.TargetRegionsData.TMB_GENE_EXCLUSIONS;

import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.purple.config.TargetRegionsData;

public class TumorMutationalLoad
{
    private final TargetRegionsData mTargetRegions;
    private int mLoad;
    private int mBurden;

    public TumorMutationalLoad(final TargetRegionsData targetRegions)
    {
        mTargetRegions = targetRegions;
        mLoad = 0;
        mBurden = 0;
    }

    public int load()
    {
        return mLoad;
    }
    public int burden() { return mBurden; }

    public int tml() { return mTargetRegions.calcTml(mLoad); }

    public double burdenPerMb() { return mBurden / MB_PER_GENOME; }

    public void processVariant(final SomaticVariant variant)
    {
        final VariantImpact variantImpact = variant.variantImpact();

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

            // - VCN <= Major Allele CN + min(20%,0.5)
            double majorAlleleCn = variant.decorator().adjustedCopyNumber() - variant.decorator().minorAlleleCopyNumber();
            double diffThreshold = min(TARGET_REGIONS_CN_DIFF, majorAlleleCn * TARGET_REGIONS_CN_PERC_DIFF);
            if(variant.copyNumber() > majorAlleleCn + diffThreshold)
                return;

            PPL_LOGGER.debug("var({}) copyNumber({}) for target-regions TMB",
                    variant.toString(), format("vcn=%.2f majorCn=%.2f", variant.copyNumber(), majorAlleleCn));
        }

        mBurden++;

        if(variantImpact.WorstCodingEffect.equals(CodingEffect.MISSENSE))
            mLoad++;
    }
}
