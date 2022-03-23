package com.hartwig.hmftools.purple.germline;

import java.util.List;
import java.util.Set;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment;
import com.hartwig.hmftools.common.variant.enrich.VariantContextEnrichment;
import com.hartwig.hmftools.purple.somatic.VariantHotspotEnrichment;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.purple.config.ReferenceData;
import com.hartwig.hmftools.purple.somatic.SnpEffEnrichment;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class GermlineVariantEnrichment implements VariantContextEnrichment
{
    private final GermlinePurityEnrichment mPurityEnrichment;
    private final SomaticRefContextEnrichment mRefGenomeEnrichment;
    private final SnpEffEnrichment mSnpEffEnrichment;
    private final GermlineReportedEnrichment mReportableEnrichment;
    private final VariantHotspotEnrichment mHotspotEnrichment;
    private final GermlineGenotypeEnrichment mGenotypeEnrichment;
    private final GermlineRescueLowVAF mLowVafRescueEnrichment;

    private final Consumer<VariantContext> mConsumer;

    public GermlineVariantEnrichment(
            final String purpleVersion, final String referenceSample, final String tumorSample, final ReferenceData refData,
            final PurityAdjuster purityAdjuster, final List<PurpleCopyNumber> copyNumbers,
            final Multimap<Chromosome, VariantHotspot> germlineHotspots, final Set<String> somaticReportedGenes,
            boolean snpEffEnrichmentEnabled, final Consumer<VariantContext> consumer)
    {
        mConsumer = consumer;

        final Set<String> germlineGenes = refData.DriverGenes.driverGenes().stream()
                .filter(DriverGene::reportGermline).map(DriverGene::gene).collect(Collectors.toSet());

        mReportableEnrichment = new GermlineReportedEnrichment(refData.DriverGenes.driverGenes(), somaticReportedGenes);
        mRefGenomeEnrichment = new SomaticRefContextEnrichment(refData.RefGenome, null);

        if(snpEffEnrichmentEnabled)
        {
            mSnpEffEnrichment = new SnpEffEnrichment(germlineGenes, refData.GeneTransCache, refData.OtherReportableTranscripts);
        }
        else
        {
            mSnpEffEnrichment = null;
        }

        mLowVafRescueEnrichment = new GermlineRescueLowVAF(referenceSample);

        mPurityEnrichment = new GermlinePurityEnrichment(purpleVersion, tumorSample, referenceSample, purityAdjuster, copyNumbers);

        mHotspotEnrichment = new VariantHotspotEnrichment(germlineHotspots, true);
        mGenotypeEnrichment = new GermlineGenotypeEnrichment(referenceSample, tumorSample);
    }

    @Override
    public void accept(final VariantContext context)
    {
        // the order matters
        VariantContext newContext = mGenotypeEnrichment.processVariant(context);

        mHotspotEnrichment.processVariant(newContext);

        mPurityEnrichment.processVariant(newContext);
        newContext = mLowVafRescueEnrichment.processVariant(newContext);

        newContext = GermlineLowTumorVCNFilter.processVariant(newContext);

        if(mSnpEffEnrichment != null)
            mSnpEffEnrichment.processVariant(newContext);

        mRefGenomeEnrichment.processVariant(newContext);
        GermlinePathogenicEnrichment.processVariant(newContext);
        mReportableEnrichment.processVariant(newContext);

        mConsumer.accept(newContext);
    }

    @Override
    public void flush()
    {
        mReportableEnrichment.flush();
    }

    @Override
    public VCFHeader enrichHeader(final VCFHeader template)
    {
        VCFHeader header = mPurityEnrichment.enrichHeader(template);
        header = mHotspotEnrichment.enrichHeader(header);
        header = mRefGenomeEnrichment.enrichHeader(header);

        if(mSnpEffEnrichment != null)
            header = mSnpEffEnrichment.enrichHeader(header);

        header = GermlineLowTumorVCNFilter.enrichHeader(header);
        header = mReportableEnrichment.enrichHeader(header);
        header = mGenotypeEnrichment.enrichHeader(header);
        return GermlinePathogenicEnrichment.enrichHeader(header);
    }
}
