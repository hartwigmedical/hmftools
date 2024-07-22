package com.hartwig.hmftools.purple.germline;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.purple.fitting.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.purple.somatic.HotspotEnrichment;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.purple.ReferenceData;

import org.jetbrains.annotations.Nullable;

import htsjdk.variant.vcf.VCFHeader;

public class GermlineVariantEnrichment
{
    private final GermlinePurityEnrichment mPurityEnrichment;
    private final GermlineReportedEnrichment mReportableEnrichment;
    private final HotspotEnrichment mHotspotEnrichment;
    private final GermlineGenotypeEnrichment mGenotypeEnrichment;
    private final GermlineRescueLowVAF mLowVafRescueEnrichment;

    public GermlineVariantEnrichment(
            final String purpleVersion, final String referenceSample, final String tumorSample, final ReferenceData refData,
            @Nullable final PurityAdjuster purityAdjuster, final List<PurpleCopyNumber> copyNumbers,
            final Multimap<Chromosome, VariantHotspot> germlineHotspots, final Set<String> somaticReportedGenes)
    {
        mReportableEnrichment = new GermlineReportedEnrichment(refData.DriverGenes.driverGenes(), somaticReportedGenes);

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

        GermlinePathogenicEnrichment.processVariant(variant.context());
        mReportableEnrichment.processVariant(variant);
    }

    public void flush()
    {
        mReportableEnrichment.flush();
    }

    public void enrichHeader(final VCFHeader header)
    {
        mPurityEnrichment.enrichHeader(header);
        mHotspotEnrichment.enrichHeader(header);
        GermlineLowTumorVCNFilter.enrichHeader(header);
        mReportableEnrichment.enrichHeader(header);
        mGenotypeEnrichment.enrichHeader(header);
        GermlinePathogenicEnrichment.enrichHeader(header);
    }
}
