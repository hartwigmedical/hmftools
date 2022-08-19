package com.hartwig.hmftools.purple.germline;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.purple.purity.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment;
import com.hartwig.hmftools.purple.somatic.HotspotEnrichment;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.purple.config.ReferenceData;

import org.jetbrains.annotations.Nullable;

import htsjdk.variant.vcf.VCFHeader;

public class GermlineVariantEnrichment
{
    private final GermlinePurityEnrichment mPurityEnrichment;
    private final SomaticRefContextEnrichment mRefGenomeEnrichment;
    private final GermlineReportedEnrichment mReportableEnrichment;
    private final HotspotEnrichment mHotspotEnrichment;
    private final GermlineGenotypeEnrichment mGenotypeEnrichment;
    private final GermlineRescueLowVAF mLowVafRescueEnrichment;

    public GermlineVariantEnrichment(
            final String purpleVersion, final String referenceSample, final String tumorSample, final ReferenceData refData,
            @Nullable final PurityAdjuster purityAdjuster, final List<PurpleCopyNumber> copyNumbers,
            final Multimap<Chromosome, VariantHotspot> germlineHotspots, final Set<String> somaticReportedGenes)
    {
        final Set<String> germlineGenes = refData.DriverGenes.driverGenes().stream()
                .filter(DriverGene::reportGermline).map(DriverGene::gene).collect(Collectors.toSet());

        mReportableEnrichment = new GermlineReportedEnrichment(refData.DriverGenes.driverGenes(), somaticReportedGenes);
        mRefGenomeEnrichment = new SomaticRefContextEnrichment(refData.RefGenome, null);

        mLowVafRescueEnrichment = new GermlineRescueLowVAF(referenceSample);

        mPurityEnrichment = new GermlinePurityEnrichment(purpleVersion, tumorSample, referenceSample, purityAdjuster, copyNumbers);

        mHotspotEnrichment = new HotspotEnrichment(germlineHotspots, true);
        mGenotypeEnrichment = new GermlineGenotypeEnrichment(referenceSample, tumorSample);
    }

    public void enrichVariant(final GermlineVariant variant)
    {
        // enrich the variant's original context so subsequent calls to get copy number and hotspot status are valid
        mGenotypeEnrichment.processVariant(variant);

        mHotspotEnrichment.processVariant(variant.context());

        mPurityEnrichment.processVariant(variant.context());

        mLowVafRescueEnrichment.processVariant(variant);

        GermlineLowTumorVCNFilter.processVariant(variant);

        mRefGenomeEnrichment.processVariant(variant.context());
        GermlinePathogenicEnrichment.processVariant(variant.context());
        mReportableEnrichment.processVariant(variant);
    }

    public void flush()
    {
        mReportableEnrichment.flush();
    }

    public VCFHeader enrichHeader(final VCFHeader template)
    {
        VCFHeader header = mPurityEnrichment.enrichHeader(template);
        header = mHotspotEnrichment.enrichHeader(header);
        header = mRefGenomeEnrichment.enrichHeader(header);

        header = GermlineLowTumorVCNFilter.enrichHeader(header);
        header = mReportableEnrichment.enrichHeader(header);
        header = mGenotypeEnrichment.enrichHeader(header);
        return GermlinePathogenicEnrichment.enrichHeader(header);
    }
}
