package com.hartwig.hmftools.purple.germline;

import java.util.List;
import java.util.Set;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.variant.enrich.GermlineGenotypeEnrichment;
import com.hartwig.hmftools.common.variant.enrich.GermlineLowTumorVCNEnrichment;
import com.hartwig.hmftools.common.variant.enrich.GermlinePathogenicEnrichment;
import com.hartwig.hmftools.common.variant.enrich.GermlinePurityEnrichment;
import com.hartwig.hmftools.common.variant.enrich.GermlineReportedEnrichment;
import com.hartwig.hmftools.common.variant.enrich.GermlineRescueLowVAFEnrichment;
import com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment;
import com.hartwig.hmftools.common.variant.enrich.VariantContextEnrichment;
import com.hartwig.hmftools.common.variant.enrich.VariantHotspotEnrichment;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffEnrichment;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
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

    public GermlineVariantEnrichment(@NotNull final String purpleVersion, @NotNull final String referenceSample,
            @NotNull final String tumorSample, @NotNull final IndexedFastaSequenceFile reference,
            @NotNull final PurityAdjuster purityAdjuster, @NotNull final List<PurpleCopyNumber> copyNumbers,
            @NotNull final DriverGenePanel genePanel, @NotNull final List<HmfTranscriptRegion> transcripts,
            @NotNull final Multimap<Chromosome, VariantHotspot> germlineHotspots, @NotNull final Set<String> somaticReportedGenes,
            @NotNull final Consumer<VariantContext> consumer)
    {
        final Set<String> germlineGenes =
                genePanel.driverGenes().stream().filter(DriverGene::reportGermline).map(DriverGene::gene).collect(Collectors.toSet());

        // Hotspot must be before reportable
        mReportableEnrichment = new GermlineReportedEnrichment(genePanel.driverGenes(), somaticReportedGenes, consumer);
        mPathogenicEnrichment = new GermlinePathogenicEnrichment(mReportableEnrichment);
        mRefGenomeEnrichment = new SomaticRefContextEnrichment(reference, mPathogenicEnrichment);

        mSnpEffEnrichment = new SnpEffEnrichment(germlineGenes, transcripts, mRefGenomeEnrichment);

        // Purity must go before lowTumorVCNEnrichment
        // Hotspot must be before lowTumorVCNEnrichment
        mLowTumorVCNEnrichment = new GermlineLowTumorVCNEnrichment(mSnpEffEnrichment);

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
    public void accept(@NotNull final VariantContext context)
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
        mSnpEffEnrichment.flush();
        mRefGenomeEnrichment.flush();
        mPathogenicEnrichment.flush();
        mReportableEnrichment.flush();
    }

    @NotNull
    @Override
    public VCFHeader enrichHeader(@NotNull final VCFHeader template)
    {
        VCFHeader header = mPurityEnrichment.enrichHeader(template);
        header = mHotspotEnrichment.enrichHeader(header);
        header = mRefGenomeEnrichment.enrichHeader(header);
        header = mSnpEffEnrichment.enrichHeader(header);
        header = mLowVafRescueEnrichment.enrichHeader(header);
        header = mReportableEnrichment.enrichHeader(header);
        header = mGenotypeEnrichment.enrichHeader(header);
        return mPathogenicEnrichment.enrichHeader(header);
    }
}
