package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.common.variant.snpeff.SnpEffUtils.SNPEFF_CANONICAL;
import static com.hartwig.hmftools.common.variant.snpeff.SnpEffUtils.SNPEFF_WORST;

import java.util.List;
import java.util.Set;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.variant.enrich.VariantContextEnrichment;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffVariantImpact;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

public class SnpEffEnrichment implements VariantContextEnrichment
{
    private final Consumer<VariantContext> mConsumer;
    private final SnpEffVariantImpact mSnpEffVariantImapct;

    public SnpEffEnrichment(
            @NotNull final Set<String> driverGenes, @NotNull final List<HmfTranscriptRegion> transcripts,
            @NotNull final Consumer<VariantContext> consumer)
    {
        mConsumer = consumer;
        mSnpEffVariantImapct = new SnpEffVariantImpact(driverGenes, transcripts);
    }

    @Override
    public void accept(@NotNull final VariantContext context)
    {
        final VariantImpact variantImpact = mSnpEffVariantImapct.fromSnpEffAnnotations(context);

        if(!variantImpact.WorstGene.isEmpty())
        {
            context.getCommonInfo().putAttribute(SNPEFF_WORST, VariantImpactSerialiser.worstDetails(variantImpact), true);
        }
        if(!variantImpact.CanonicalGene.isEmpty())
        {
            context.getCommonInfo().putAttribute(SNPEFF_CANONICAL, VariantImpactSerialiser.canonicalDetails(variantImpact), true);
        }
        mConsumer.accept(context);
    }

    @Override
    public void flush()
    {
        // Do nothing
    }

    @NotNull
    @Override
    public VCFHeader enrichHeader(@NotNull final VCFHeader header)
    {
        header.addMetaDataLine(new VCFInfoHeaderLine(
                SNPEFF_WORST, 5, VCFHeaderLineType.String,
                "SnpEff worst transcript summary [Gene, Transcript, Effect, CodingEffect, GenesAffected]"));

        header.addMetaDataLine(new VCFInfoHeaderLine(
                SNPEFF_CANONICAL, 6, VCFHeaderLineType.String,
                "SnpEff canonical transcript summary [Gene, Transcript, Effect, CodingEffect, HgvsCodingImpact, HgvsProteinImpact]"));

        return header;
    }
}
