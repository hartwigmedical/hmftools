package com.hartwig.hmftools.pave.reverse;

import static com.hartwig.hmftools.pave.reverse.ReversePaveConfig.RPV_LOGGER;
import static com.hartwig.hmftools.pave.reverse.ReversePaveConstants.APP_NAME;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class ReversePaveApplication
{
    private final ReversePaveConfig mConfig;

    public ReversePaveApplication(final ConfigBuilder configBuilder)
    {
        mConfig = new ReversePaveConfig(configBuilder);
    }

    public void run()
    {
        ReversePave reversePave = new ReversePave(mConfig.mEnsemblCache, mConfig.mRefGenome);
        RoundTripChecker checker = new RoundTripChecker(reversePave);

        VcfFileReader vcfFileReader = new VcfFileReader(mConfig.mVcfFile, true);
        VCFHeader inputHeader = vcfFileReader.vcfHeader();

        for(VariantContext context : vcfFileReader.iterator())
        {
            VariantContext fullContext = context.fullyDecode(inputHeader, false);
            VariantImpact variantImpact = VariantImpactSerialiser.fromVariantContext(fullContext);

            //            if(variantImpact != null && variantImpact.CanonicalEffect.equals("missense_variant"))
            final String proteinVariant = variantImpact.CanonicalHgvsProtein;
            if(proteinVariant.contains("p.") && !proteinVariant.endsWith("="))
            {
                checker.compareActualChangesWithCalculated(fullContext, variantImpact);
            }
        }
        vcfFileReader.close();
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        ReversePaveConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        ReversePaveApplication reversePaveApplication = new ReversePaveApplication(configBuilder);
        reversePaveApplication.run();
    }
}

class RoundTripChecker
{
    @NotNull
    private final ReversePave reversePave;

    RoundTripChecker(@NotNull final ReversePave reversePave)
    {
        this.reversePave = reversePave;
    }

    void compareActualChangesWithCalculated(VariantContext context, VariantImpact impact)
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
            calculatedVariants = reversePave.calculateVariant(gene, transcript, variant);
        }
        catch(Exception e)
        {
            RPV_LOGGER.error("Failed to compute variants for " + gene + ", " + variant + " and " + transcript, e);
            return;
        }
        boolean calculatedContainsActual = calculatedVariants.mChanges.contains(actualChange);
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