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

    public void compareActualChangesWithCalculated(VariantContext context, VariantImpact impact)
    {
        if(shouldIgnoreVariant(impact))
        {
            return;
        }
        int start = context.getStart();
        String ref = context.getReference().getBaseString();
        String alt = context.getAltAlleleWithHighestAlleleCount().getBaseString();
        BaseSequenceChange actualChange = new BaseSequenceChange(ref, alt, context.getContig(), start);

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
            RPV_LOGGER.info("Ok: " + actualChange);
        }
    }

    private boolean shouldIgnoreVariant(VariantImpact variantImpact)
    {
        final String proteinVariant = variantImpact.CanonicalHgvsProtein;
        if(!proteinVariant.startsWith("p."))
        {
            return true;
        }
        if(proteinVariant.endsWith("="))
        {
            return true;
        }
        if(proteinVariant.equals("p.?"))
        {
            return true;
        }
        return false;
    }
}
