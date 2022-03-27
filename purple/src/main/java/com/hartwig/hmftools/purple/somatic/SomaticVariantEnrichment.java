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
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class SomaticVariantEnrichment
{
    private final SomaticPurityEnrichment mPurityEnrichment;
    private final HotspotEnrichment mHotspotEnrichment;
    private final KataegisEnrichment mKataegisEnrichment;
    private final SomaticRefContextEnrichment mSomaticRefContextEnrichment;
    private final SubclonalLikelihoodEnrichment mSubclonalLikelihoodEnrichment;
    private final SnpEffEnrichment mSnpEffEnrichment;
    private final SomaticGenotypeEnrichment mGenotypeEnrichment;

    public SomaticVariantEnrichment(
            boolean hotspotEnabled, boolean snpEffEnrichmentEnabled, double clonalityBinWidth, final String purpleVersion,
            final String referenceId, final String tumorSample, final ReferenceData refData,
            final PurityAdjuster purityAdjuster, final List<PurpleCopyNumber> copyNumbers, final List<FittedRegion> fittedRegions,
            final Multimap<Chromosome, VariantHotspot> hotspots, final List<PeakModel> peakModel)
    {
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

        mHotspotEnrichment = new HotspotEnrichment(hotspots, hotspotEnabled);
    }

    public void enrich(final SomaticData variant)
    {
        VariantContext newContext = variant.newContext();

        if(mSnpEffEnrichment != null)
            mSnpEffEnrichment.processVariant(newContext);

        mSomaticRefContextEnrichment.processVariant(newContext);

        mKataegisEnrichment.processVariant(newContext);

        // has occurred earlier now
        // mPurityEnrichment.processVariant(newContext);

        mSubclonalLikelihoodEnrichment.processVariant(newContext);

        mGenotypeEnrichment.processVariant(newContext);
    }

    public void flush()
    {
        mKataegisEnrichment.flush();
    }

    public VCFHeader populateHeader(final VCFHeader template)
    {
        VCFHeader header = SomaticRefContextEnrichment.addHeader(template);
        header = KataegisEnrichment.enrichHeader(header);
        header = SubclonalLikelihoodEnrichment.enrichHeader(header);
        header = HotspotEnrichment.enrichHeader(header);

        if(mSnpEffEnrichment != null)
            header = SnpEffEnrichment.enrichHeader(header);

        return mPurityEnrichment.enrichHeader(header);
    }

}
