package com.hartwig.hmftools.purple.somatic;

import java.util.List;
import java.util.Set;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.purple.config.ReferenceData;
import com.hartwig.hmftools.purple.fitting.PeakModel;
import com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment;
import com.hartwig.hmftools.common.variant.enrich.VariantContextEnrichment;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFHeader;

public class SomaticVariantEnrichment implements VariantContextEnrichment
{
    private final SomaticPurityEnrichment mPurityEnrichment;
    private final VariantHotspotEnrichment mHotspotEnrichment;
    private final KataegisEnrichment mKataegisEnrichment;
    private final SomaticRefContextEnrichment mSomaticRefContextEnrichment;
    private final SubclonalLikelihoodEnrichment mSubclonalLikelihoodEnrichment;
    private final SnpEffEnrichment mSnpEffEnrichment;
    private final SomaticGenotypeEnrichment mGenotypeEnrichment;

    private final Consumer<VariantContext> mConsumer;

    public SomaticVariantEnrichment(
            boolean hotspotEnabled, boolean snpEffEnrichmentEnabled, double clonalityBinWidth, final String purpleVersion,
            final String referenceId, final String tumorSample, final ReferenceData refData,
            final PurityAdjuster purityAdjuster, final List<PurpleCopyNumber> copyNumbers, final List<FittedRegion> fittedRegions,
            final Multimap<Chromosome, VariantHotspot> hotspots, final List<PeakModel> peakModel, final Consumer<VariantContext> consumer)
    {
        mConsumer = consumer;

        mGenotypeEnrichment = new SomaticGenotypeEnrichment(referenceId, tumorSample);

        mSubclonalLikelihoodEnrichment = new SubclonalLikelihoodEnrichment(clonalityBinWidth, peakModel);

        mPurityEnrichment = new SomaticPurityEnrichment(purpleVersion, tumorSample, purityAdjuster, copyNumbers, fittedRegions);

        mKataegisEnrichment = new KataegisEnrichment();

        mSomaticRefContextEnrichment = new SomaticRefContextEnrichment(refData.RefGenome, null);

        if(snpEffEnrichmentEnabled)
        {
            final Set<String> somaticGenes = refData.DriverGenes.driverGenes().stream()
                    .filter(DriverGene::reportSomatic).map(DriverGene::gene).collect(Collectors.toSet());

            mSnpEffEnrichment = new SnpEffEnrichment(somaticGenes, refData.GeneTransCache, refData.OtherReportableTranscripts);
        }
        else
        {
            mSnpEffEnrichment = null;
        }

        mHotspotEnrichment = new VariantHotspotEnrichment(hotspots, hotspotEnabled);
    }

    @Override
    public void accept(@NotNull final VariantContext context)
    {
        VariantContext newContext = new VariantContextBuilder(context).make();

        mHotspotEnrichment.processVariant(newContext);

        if(mSnpEffEnrichment != null)
            mSnpEffEnrichment.processVariant(newContext);

        mSomaticRefContextEnrichment.processVariant(newContext);

        mKataegisEnrichment.processVariant(newContext);

        mPurityEnrichment.processVariant(newContext);

        mSubclonalLikelihoodEnrichment.processVariant(newContext);

        mGenotypeEnrichment.processVariant(newContext);

        mConsumer.accept(newContext);
    }

    @Override
    public void flush()
    {
        mKataegisEnrichment.flush();
    }

    @Override
    public VCFHeader enrichHeader(final VCFHeader template)
    {
        VCFHeader header = SomaticRefContextEnrichment.addHeader(template);
        header = KataegisEnrichment.enrichHeader(header);
        header = SubclonalLikelihoodEnrichment.enrichHeader(header);
        header = VariantHotspotEnrichment.enrichHeader(header);

        if(mSnpEffEnrichment != null)
            header = SnpEffEnrichment.enrichHeader(header);

        return mPurityEnrichment.enrichHeader(header);
    }

}
