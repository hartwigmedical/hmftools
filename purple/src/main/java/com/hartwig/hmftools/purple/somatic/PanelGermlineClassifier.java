package com.hartwig.hmftools.purple.somatic;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsDouble;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PANEL_GERMLINE_VAF_DISTANCE;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PANEL_GERMLINE_VAF_DISTANCE_DESC;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PANEL_SOMATIC_LIKELIHOOD;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PANEL_SOMATIC_LIKELIHOOD_DESC;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.UNCLEAR_GERMLINE_FLAG;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.UNCLEAR_GERMLINE_FLAG_DESCRIPTION;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.config.TargetRegionsData.PANEL_SOMATIC_LIKELIHOOD_DIFF_HIGH;
import static com.hartwig.hmftools.purple.config.TargetRegionsData.PANEL_SOMATIC_LIKELIHOOD_DIFF_LOW;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.variant.PanelSomaticLikelihood;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.purple.config.PurpleConfig;
import com.hartwig.hmftools.purple.config.TargetRegionsData;

import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class PanelGermlineClassifier
{
    private final TargetRegionsData mTargetRegions;
    private final BufferedWriter mWriter;

    public PanelGermlineClassifier(final TargetRegionsData targetRegions, final PurpleConfig config)
    {
        mTargetRegions = targetRegions;

        // mWriter = targetRegions.hasTargetRegions() ? initialiseWriter(config) : null;
        mWriter = null; // disabled, was for temporary analysis
    }

    public static VCFHeader enrichHeader(final VCFHeader template)
    {
        template.addMetaDataLine(new VCFInfoHeaderLine(UNCLEAR_GERMLINE_FLAG, 0, VCFHeaderLineType.Flag, UNCLEAR_GERMLINE_FLAG_DESCRIPTION));
        template.addMetaDataLine(new VCFInfoHeaderLine(PANEL_GERMLINE_VAF_DISTANCE, 2, VCFHeaderLineType.Float, PANEL_GERMLINE_VAF_DISTANCE_DESC));
        template.addMetaDataLine(new VCFInfoHeaderLine(PANEL_SOMATIC_LIKELIHOOD, 1, VCFHeaderLineType.String, PANEL_SOMATIC_LIKELIHOOD_DESC));
        return template;
    }

    public void processVariant(final SomaticVariant variant, double purity)
    {
        if(!mTargetRegions.hasTargetRegions())
            return;

        double rawAf = getGenotypeAttributeAsDouble(variant.context().getGenotype(0), VCFConstants.ALLELE_FREQUENCY_KEY, 0);

        double segmentCn = variant.decorator().adjustedCopyNumber();
        double tumorMinorCn = variant.decorator().minorAlleleCopyNumber();
        double tumorMajorCn = segmentCn - tumorMinorCn;

        double refPurity = 1 - purity;
        double referenceMajorCn = 1; // B2
        double referenceMinorCn = 1; // C2
        double germlineCn = referenceMajorCn + referenceMinorCn;

        // double minorVAF = (refPurity + tumorMinorCn * purity) / denom;
        // double majorVAF = (refPurity + tumorMajorCn * purity) / denom;
        // boolean isUnclearGermline = abs(majorVAF - rawAf) < mTargetRegions.maxAFDiff() || abs(minorVAF - rawAf) < mTargetRegions.maxAFDiff();

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

        boolean isUnclearGermline = abs(germlineMajorDiff) < mTargetRegions.maxAFDiff() || abs(germlineMinorDiff) < mTargetRegions.maxAFDiff();

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

        boolean closestDiffGermline = abs(germlineMinDiff) < abs(somaticMinDiff);

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

        if(isUnclearGermline)
        {
            variant.context().getCommonInfo().putAttribute(UNCLEAR_GERMLINE_FLAG, true);
        }

        PanelSomaticLikelihood somaticLikelihood;
        double somaticGermlineDiff = abs(somaticMinDiff) - abs(germlineMinDiff);

        if(variant.isHotspot() || somaticGermlineDiff < PANEL_SOMATIC_LIKELIHOOD_DIFF_HIGH)
        {
            somaticLikelihood = PanelSomaticLikelihood.HIGH;
        }
        else if(somaticGermlineDiff > PANEL_SOMATIC_LIKELIHOOD_DIFF_LOW)
        {
            somaticLikelihood = PanelSomaticLikelihood.LOW;
        }
        else
        {
            somaticLikelihood = PanelSomaticLikelihood.MEDIUM;
        }

        variant.context().getCommonInfo().putAttribute(PANEL_SOMATIC_LIKELIHOOD, somaticLikelihood);

        writeVariantData(
                variant, purity, rawAf, tumorMajorCn, tumorMinorCn, germlineVafMajor, germlineVafMinor, somaticVafMajor, somaticVafMinor,
                closestDiffGermline, closestDiffGermline ? germlineMinDiff : somaticMinDiff);
    }

    private void writeVariantData(
            final SomaticVariant variant, double purity, double rawAf, double tumorMajorCn, double tumorMinorCn, double germlineVafMajor,
            double germlineVafMinor, double somaticVafMajor, double somaticVafMinor, boolean closestDiffGermline, double minDiff)
    {
        if(mWriter == null)
            return;

        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(format("%s:%s %s>%s", variant.chromosome(), variant.position(), variant.decorator().ref(), variant.decorator().alt()));
            sj.add(variant.type().toString());

            VariantImpact variantImpact = variant.variantImpact();

            if(variantImpact.CanonicalGeneName.isEmpty())
                sj.add("");
            else
                sj.add(format("%s:%s", variant.variantImpact().CanonicalGeneName, variant.variantImpact().CanonicalCodingEffect));

            sj.add(format("%.3f", purity));
            sj.add(format("%.3f", rawAf));
            sj.add(format("%.3f", tumorMajorCn));
            sj.add(format("%.3f", tumorMinorCn));
            sj.add(format("%.3f", germlineVafMajor));
            sj.add(format("%.3f", germlineVafMinor));
            sj.add(format("%.3f", somaticVafMajor));
            sj.add(format("%.3f", somaticVafMinor));
            sj.add(format("%s", closestDiffGermline ? "germline" : "somatic"));
            sj.add(format("%.3f", minDiff));

            mWriter.write(sj.toString());
            mWriter.newLine();
        }
        catch(IOException e)
        {
            PPL_LOGGER.error("failed to write panel germline file: {}", e.toString());
        }
    }

    private BufferedWriter initialiseWriter(final PurpleConfig config)
    {
        try
        {
            String fileName = config.OutputDir + config.TumorId + ".purple.panel_germline_data.tsv";

            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.write("VarInfo\tVarType\tGeneInfo\tPurity\tRawAf\tTumorMajorCn\tTumorMinorCn");
            writer.write("\tGermlineVafMajor\tGermlineVafMinor\tSomaticVafMajor\tSomaticVafMinor\tMinStatus\tMinDiff");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            PPL_LOGGER.error("failed to initialise panel germline file: {}", e.toString());
            return null;
        }
    }

    public void close() { closeBufferedWriter(mWriter); }
}
