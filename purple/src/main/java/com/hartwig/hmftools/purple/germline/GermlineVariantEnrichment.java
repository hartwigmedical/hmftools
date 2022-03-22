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

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class GermlineVariantEnrichment implements VariantContextEnrichment
{
    private final VariantContextEnrichment mPurityEnrichment;
    private final VariantContextEnrichment mRefGenomeEnrichment;
    private final VariantContextEnrichment mPathogenicEnrichment;
    private final VariantContextEnrichment mSnpEffEnrichment;
    private final VariantContextEnrichment mReportableEnrichment;
    private final VariantContextEnrichment mHotspotEnrichment;
    private final VariantContextEnrichment mGenotypeEnrichment;
    private final VariantContextEnrichment mLowTumorVCNEnrichment;
    private final VariantContextEnrichment mLowVafRescueEnrichment;

    public GermlineVariantEnrichment(
            final String purpleVersion, final String referenceSample, final String tumorSample, final ReferenceData refData,
            final PurityAdjuster purityAdjuster, final List<PurpleCopyNumber> copyNumbers,
            final Multimap<Chromosome, VariantHotspot> germlineHotspots, final Set<String> somaticReportedGenes,
            boolean snpEffEnrichmentEnabled, final Consumer<VariantContext> consumer)
    {
        final Set<String> germlineGenes = refData.DriverGenes.driverGenes().stream()
                .filter(DriverGene::reportGermline).map(DriverGene::gene).collect(Collectors.toSet());

        // Hotspot must be before reportable
        mReportableEnrichment = new GermlineReportedEnrichment(refData.DriverGenes.driverGenes(), somaticReportedGenes, consumer);
        mPathogenicEnrichment = new GermlinePathogenicEnrichment(mReportableEnrichment);
        mRefGenomeEnrichment = new SomaticRefContextEnrichment(refData.RefGenome, mPathogenicEnrichment);

        if(snpEffEnrichmentEnabled)
        {
            mSnpEffEnrichment =
                    new SnpEffEnrichment(germlineGenes, refData.GeneTransCache, refData.OtherReportableTranscripts, mRefGenomeEnrichment);
        }
        else
        {
            mSnpEffEnrichment = null;
        }

        VariantContextEnrichment prevConsumer = mSnpEffEnrichment != null ? mSnpEffEnrichment : mRefGenomeEnrichment;

        // Purity must go before lowTumorVCNEnrichment
        // Hotspot must be before lowTumorVCNEnrichment
        mLowTumorVCNEnrichment = new GermlineLowTumorVCNEnrichment(prevConsumer);

        // Purity must go before lowVafRescue
        // Genotype must go before lowVafRescue
        mLowVafRescueEnrichment = new GermlineRescueLowVAFEnrichment(referenceSample, mLowTumorVCNEnrichment);

        // Genotype must go before purity enrichment
        mPurityEnrichment = new GermlinePurityEnrichment(purpleVersion,
                tumorSample,
                referenceSample,
                purityAdjuster,
                copyNumbers,
                mLowVafRescueEnrichment);

        mHotspotEnrichment = new VariantHotspotEnrichment(germlineHotspots, mPurityEnrichment);
        mGenotypeEnrichment = new GermlineGenotypeEnrichment(referenceSample, tumorSample, mHotspotEnrichment);
    }

    @Override
    public void accept(final VariantContext context)
    {
        mGenotypeEnrichment.accept(context);
    }

    @Override
    public void flush()
    {
        mGenotypeEnrichment.flush();
        mHotspotEnrichment.flush();
        mPurityEnrichment.flush();
        mLowVafRescueEnrichment.flush();
        mLowTumorVCNEnrichment.flush();

        if(mSnpEffEnrichment != null)
        {
            mSnpEffEnrichment.flush();
        }

        mRefGenomeEnrichment.flush();
        mPathogenicEnrichment.flush();
        mReportableEnrichment.flush();
    }

    @NotNull
    @Override
    public VCFHeader enrichHeader(final VCFHeader template)
    {
        VCFHeader header = mPurityEnrichment.enrichHeader(template);
        header = mHotspotEnrichment.enrichHeader(header);
        header = mRefGenomeEnrichment.enrichHeader(header);

        if(mSnpEffEnrichment != null)
        {
            header = mSnpEffEnrichment.enrichHeader(header);
        }

        header = mLowVafRescueEnrichment.enrichHeader(header);
        header = mLowTumorVCNEnrichment.enrichHeader(header);
        header = mReportableEnrichment.enrichHeader(header);
        header = mGenotypeEnrichment.enrichHeader(header);
        return mPathogenicEnrichment.enrichHeader(header);
    }
}
