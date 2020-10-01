package com.hartwig.hmftools.common.variant.enrich;

import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.region.CanonicalTranscript;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.variant.clonality.PeakModel;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class VariantContextEnrichmentPurple implements VariantContextEnrichment {

    private final VariantContextEnrichment purityEnrichment;
    private final VariantContextEnrichment hotspotEnrichment;
    private final VariantContextEnrichment kataegisEnrichment;
    private final VariantContextEnrichment somaticRefContextEnrichment;
    private final VariantContextEnrichment subclonalLikelihoodEnrichment;
    private final VariantContextEnrichment snpEffEnrichment;

    public VariantContextEnrichmentPurple(boolean hotspotEnabled, double clonalityMaxPloidy, double clonalityBinWidth,
            @NotNull final String purpleVersion, @NotNull final String tumorSample, @NotNull final IndexedFastaSequenceFile reference,
            @NotNull final PurityAdjuster purityAdjuster, @NotNull final DriverGenePanel genePanel,
            @NotNull final List<PurpleCopyNumber> copyNumbers, @NotNull final List<FittedRegion> fittedRegions,
            @NotNull final List<PeakModel> peakModel, @NotNull final Multimap<Chromosome, VariantHotspot> hotspots, @NotNull final List<CanonicalTranscript> transcripts,
            @NotNull final Consumer<VariantContext> consumer) {
        subclonalLikelihoodEnrichment = new SubclonalLikelihoodEnrichment(clonalityMaxPloidy, clonalityBinWidth, peakModel, consumer);
        purityEnrichment =
                new PurityEnrichment(purpleVersion, tumorSample, purityAdjuster, copyNumbers, fittedRegions, subclonalLikelihoodEnrichment);
        kataegisEnrichment = new KataegisEnrichment(purityEnrichment);
        somaticRefContextEnrichment = new SomaticRefContextEnrichment(reference, kataegisEnrichment);
        snpEffEnrichment = new SnpEffEnrichment(genePanel.driverGenes(), transcripts, somaticRefContextEnrichment);
        if (hotspotEnabled) {
            hotspotEnrichment = new VariantHotspotEnrichment(hotspots, snpEffEnrichment);
        } else {
            hotspotEnrichment = VariantContextEnrichmentFactory.noEnrichment().create(snpEffEnrichment);
        }
    }

    @Override
    public void flush() {
        hotspotEnrichment.flush();
        snpEffEnrichment.flush();
        somaticRefContextEnrichment.flush();
        kataegisEnrichment.flush();
        purityEnrichment.flush();
        subclonalLikelihoodEnrichment.flush();
    }

    @NotNull
    @Override
    public VCFHeader enrichHeader(@NotNull final VCFHeader template) {
        VCFHeader header = somaticRefContextEnrichment.enrichHeader(template);
        header = kataegisEnrichment.enrichHeader(header);
        header = subclonalLikelihoodEnrichment.enrichHeader(header);
        header = hotspotEnrichment.enrichHeader(header);
        header = snpEffEnrichment.enrichHeader(header);
        return purityEnrichment.enrichHeader(header);
    }

    @Override
    public void accept(@NotNull final VariantContext context) {
        hotspotEnrichment.accept(context);
    }
}
