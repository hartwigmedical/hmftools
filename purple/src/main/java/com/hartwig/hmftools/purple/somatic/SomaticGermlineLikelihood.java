package com.hartwig.hmftools.purple.somatic;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsDouble;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PANEL_GERMLINE_VAF_DISTANCE;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PANEL_GERMLINE_VAF_DISTANCE_DESC;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PANEL_SOMATIC_LIKELIHOOD;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PANEL_SOMATIC_LIKELIHOOD_DESC;
import static com.hartwig.hmftools.purple.TargetRegionsData.PANEL_SOMATIC_LIKELIHOOD_DIFF_HIGH;
import static com.hartwig.hmftools.purple.TargetRegionsData.PANEL_SOMATIC_LIKELIHOOD_DIFF_LOW;

import com.hartwig.hmftools.common.variant.GenotypeIds;
import com.hartwig.hmftools.common.variant.SomaticLikelihood;
import com.hartwig.hmftools.purple.PurpleConfig;

import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class SomaticGermlineLikelihood
{
    private final boolean mEnabled;
    private final GenotypeIds mGenotypeIds;

    public SomaticGermlineLikelihood(final PurpleConfig config, final GenotypeIds genotypeIds)
    {
        mEnabled = config.tumorOnlyMode();
        mGenotypeIds = genotypeIds;
    }

    public static VCFHeader enrichHeader(final VCFHeader template)
    {
        template.addMetaDataLine(new VCFInfoHeaderLine(PANEL_GERMLINE_VAF_DISTANCE, 2, VCFHeaderLineType.Float, PANEL_GERMLINE_VAF_DISTANCE_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PANEL_SOMATIC_LIKELIHOOD, 1, VCFHeaderLineType.String, PANEL_SOMATIC_LIKELIHOOD_DESC));
        return template;
    }

    public void processVariant(final SomaticVariant variant, double purity)
    {
        if(!mEnabled)
            return;

        double rawAf = getGenotypeAttributeAsDouble(
                variant.context().getGenotype(mGenotypeIds.TumorOrdinal), VCFConstants.ALLELE_FREQUENCY_KEY, 0);

        double segmentCn = variant.decorator().adjustedCopyNumber();
        double tumorMinorCn = variant.decorator().minorAlleleCopyNumber();
        double tumorMajorCn = segmentCn - tumorMinorCn;

        double refPurity = 1 - purity;
        double referenceMajorCn = 1; // B2
        double referenceMinorCn = 1; // C2
        double germlineCn = referenceMajorCn + referenceMinorCn;

        // B3 = tumorMinorCn
        // C3 = tumorMajorCn
        // B5 = purity
        // B6 = observed VAF = rawAf

        // ExpectedGermlineMajorVaf = (B2*(1-B5)+C3*B5)/((B2+C2)*(1-B5)+(B3+C3)*B5)
        // ExpectedGermlineMinorVaf = =(B2*(1-B5)+B3*B5)/((B2+C2)*(1-B5)+(B3+C3)*B5)

        double denom = germlineCn * refPurity + segmentCn * purity;

        double germlineVafMajor = (referenceMajorCn * refPurity + tumorMajorCn * purity) / denom;
        double germlineVafMinor = (referenceMinorCn * refPurity + tumorMinorCn * purity) / denom;
        double homozygousGermline = 1;

        double germlineMajorDiff = germlineVafMajor - rawAf;
        double germlineMinorDiff = germlineVafMinor - rawAf;
        double germlineHomDiff = homozygousGermline - rawAf;

        double germlineMinDiff = abs(germlineMajorDiff) < abs(germlineMinorDiff) ? germlineMajorDiff : germlineMinorDiff;
        germlineMinDiff = abs(germlineHomDiff) < abs(germlineMinDiff) ? germlineHomDiff : germlineMinDiff;

        // =(B5)/((B2+C2)*(1-B5)+(B3+C3)*B5)
        // somaticMinorVaf = =(B3*B5)/((B2+C2)*(1-B5)+(B3+C3)*B5)
        double somaticVaf1 = purity / denom;
        double somaticVafMajor = (tumorMajorCn * purity) / denom;
        double somaticVafMinor = (tumorMinorCn * purity) / denom;

        double somaticMajorDiff = somaticVafMajor - rawAf;
        double somaticMinorDiff = somaticVafMinor - rawAf;
        double somaticHomDiff = somaticVaf1 - rawAf;

        double somaticMinDiff = abs(somaticMajorDiff) < abs(somaticMinorDiff) ? somaticMajorDiff : somaticMinorDiff;
        somaticMinDiff = abs(somaticHomDiff) < abs(somaticMinDiff) ? somaticHomDiff : somaticMinDiff;

        if(rawAf < somaticVaf1)
            somaticMinDiff = 0; // the subclonal case

        /*
        if(PPL_LOGGER.isTraceEnabled())
        {
            PPL_LOGGER.trace(format("var(%s) af(%.2f) cn(seg=%.2f var=%.2f major=%.2f minor=%.2f) diff(som=%.3f germ=%.3f) germlineStatus(%s)",
                    variant, rawAf, segmentCn, variantCn, tumorMajorCn, tumorMinorCn, somaticMinDiff, germlineMinDiff,
                    isUnclearGermline ? "unclear" : "somatic"));
        }
        */

        double[] diffValues = new double[] { germlineMinDiff, somaticMinDiff };
        variant.context().getCommonInfo().putAttribute(PANEL_GERMLINE_VAF_DISTANCE, diffValues);

        SomaticLikelihood somaticLikelihood;
        double somaticGermlineDiff = abs(somaticMinDiff) - abs(germlineMinDiff);

        if(variant.isHotspot() || somaticGermlineDiff < PANEL_SOMATIC_LIKELIHOOD_DIFF_HIGH)
        {
            somaticLikelihood = SomaticLikelihood.HIGH;
        }
        else if(somaticGermlineDiff > PANEL_SOMATIC_LIKELIHOOD_DIFF_LOW)
        {
            somaticLikelihood = SomaticLikelihood.LOW;
        }
        else
        {
            somaticLikelihood = SomaticLikelihood.MEDIUM;
        }

        variant.context().getCommonInfo().putAttribute(PANEL_SOMATIC_LIKELIHOOD, somaticLikelihood);
    }
}
