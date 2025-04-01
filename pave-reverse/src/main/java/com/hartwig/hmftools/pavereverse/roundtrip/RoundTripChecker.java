package com.hartwig.hmftools.pavereverse.roundtrip;

import static com.hartwig.hmftools.pavereverse.ReversePaveConfig.RPV_LOGGER;

import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.pavereverse.BaseSequenceChange;
import com.hartwig.hmftools.pavereverse.BaseSequenceVariants;
import com.hartwig.hmftools.pavereverse.ReversePave;

import htsjdk.variant.variantcontext.VariantContext;

public class RoundTripChecker
{
    private final ReversePave mReversePave;

    public RoundTripChecker(ReversePave reversePave)
    {
        mReversePave = reversePave;
    }

    private int numberOfCodingVariantsOk = 0;
    private int numberOfCodingVariantsInError = 0;

    public void printResults()
    {
        RPV_LOGGER.info("Coding variants ok: " + numberOfCodingVariantsOk);
        RPV_LOGGER.info("Coding variants in error: " + numberOfCodingVariantsInError);

    }

    public void compareActualChangesWithCalculated(VariantContext context, VariantImpact impact)
    {
        boolean compareProteinImpacts = hasComparableProteinImpact(impact);
        boolean compareCoding = hasCodingImpact(impact);
        if(compareProteinImpacts || compareCoding)
        {
            int start = context.getStart();
            String ref = context.getReference().getBaseString();
            String alt = context.getAltAlleleWithHighestAlleleCount().getBaseString();
            BaseSequenceChange actualChange = new BaseSequenceChange(ref, alt, context.getContig(), start);
            if(compareProteinImpacts)
            {
                //                compareActualWithCalculatedProteinImpact(actualChange, impact);
            }
            if(compareCoding)
            {
                compareActualWithCalculatedChange(actualChange, impact);
            }
        }
    }

    private void compareActualWithCalculatedChange(BaseSequenceChange actualChange, VariantImpact impact)
    {
        String gene = impact.GeneName;
        String transcript = impact.CanonicalTranscript;
        String variant = impact.CanonicalHgvsCoding;
//        if(!variant.contains("_"))
//        {
//            return;
//        }
        //        if(variant.contains("c.-73-4038_-73-4039dupGG") && variant.contains("_"))
        //        {
        //            System.out.println(">>>>DUP!! gene: " + gene + " transcript: " + transcript + " variant: " + variant);
        //        }
        BaseSequenceChange calculated;
        try
        {
            calculated = mReversePave.calculateDnaVariant(gene, transcript, variant);
        }
        catch(Exception e)
        {
            RPV_LOGGER.error("Failed to compute DNA change for " + gene + ", " + variant + " and " + transcript, e);
            return;
        }
        if(!actualChange.equals(calculated))
        {
            String message = "Calculated change does not equal actual change.";
            message += "\nGene: " + gene + ", Transcript: " + transcript + ", Variant: " + variant;
            message += "\n    Actual: " + actualChange;
            message += "\nCalculated: " + calculated;
            RPV_LOGGER.warn(message);
            numberOfCodingVariantsInError++;
        }
        else
        {
//            RPV_LOGGER.info("Coding OK for: " + gene + " " + variant);
            numberOfCodingVariantsOk++;
        }
    }

    private void compareActualWithCalculatedProteinImpact(BaseSequenceChange actualChange, VariantImpact impact)
    {
        String gene = impact.GeneName;
        String transcript = impact.CanonicalTranscript;
        String variant = impact.CanonicalHgvsProtein;
        BaseSequenceVariants calculatedVariants;
        try
        {
            calculatedVariants = mReversePave.calculateProteinVariant(gene, transcript, variant);
        }
        catch(Exception e)
        {
            RPV_LOGGER.error("Failed to compute variants for " + gene + ", " + variant + " and " + transcript, e);
            return;
        }
        boolean calculatedContainsActual = calculatedVariants.changes().contains(actualChange);
        if(!calculatedContainsActual)
        {
            String msg =
                    String.format("Calculated changes do not include actual change. Gene: %s, transcript: %s, variant: %s, actual: %s", gene, transcript, variant, actualChange);
            RPV_LOGGER.warn(msg);
        }
        else
        {
            //            RPV_LOGGER.info("Protein OK: " + actualChange);
        }
    }

    private boolean hasCodingImpact(VariantImpact impact)
    {
        String codingHgvs = impact.CanonicalHgvsCoding;
        return codingHgvs != null && codingHgvs.startsWith("c.");
    }

    private boolean hasComparableProteinImpact(VariantImpact variantImpact)
    {
        final String proteinVariant = variantImpact.CanonicalHgvsProtein;
        if(!proteinVariant.startsWith("p."))
        {
            return false;
        }
        if(proteinVariant.endsWith("="))
        {
            return false;
        }
        if(proteinVariant.equals("p.?"))
        {
            return false;
        }
        return true;
    }
}
