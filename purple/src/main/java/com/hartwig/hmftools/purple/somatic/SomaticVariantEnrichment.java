package com.hartwig.hmftools.purple.somatic;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.purple.config.ReferenceData;
import com.hartwig.hmftools.purple.fitting.PeakModel;
import com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment;

import htsjdk.variant.vcf.VCFHeader;

public class SomaticVariantEnrichment
{
    private final SomaticPurityEnrichment mPurityEnrichment;
    private final KataegisEnrichment mKataegisEnrichment;
    private final SomaticRefContextEnrichment mSomaticRefContextEnrichment;
    private final SubclonalLikelihoodEnrichment mSubclonalLikelihoodEnrichment;
    private final SnpEffEnrichment mSnpEffEnrichment;
    private final SomaticGenotypeEnrichment mGenotypeEnrichment;

    public SomaticVariantEnrichment(
            boolean snpEffEnrichmentEnabled, double clonalityBinWidth, final String purpleVersion,
            final String referenceId, final String tumorSample, final ReferenceData refData,
            final PurityAdjuster purityAdjuster, final List<PurpleCopyNumber> copyNumbers, final List<FittedRegion> fittedRegions,
            final List<PeakModel> peakModel)
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
    }

    public void enrich(final SomaticVariant variant)
    {
        if(mSnpEffEnrichment != null)
            mSnpEffEnrichment.processVariant(variant.context());

        mSomaticRefContextEnrichment.processVariant(variant.context());

        mKataegisEnrichment.processVariant(variant.context());

        mSubclonalLikelihoodEnrichment.processVariant(variant);

        mGenotypeEnrichment.processVariant(variant);
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
