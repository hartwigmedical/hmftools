package com.hartwig.hmftools.purple.somatic;

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
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.purple.fitting.PeakModel;
import com.hartwig.hmftools.common.variant.enrich.KataegisEnrichment;
import com.hartwig.hmftools.common.variant.enrich.SnpEffEnrichment;
import com.hartwig.hmftools.common.variant.enrich.SomaticPurityEnrichment;
import com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment;
import com.hartwig.hmftools.common.variant.enrich.VariantContextEnrichment;
import com.hartwig.hmftools.common.variant.enrich.VariantContextEnrichmentFactory;
import com.hartwig.hmftools.common.variant.enrich.VariantHotspotEnrichment;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class SomaticVariantEnrichment implements VariantContextEnrichment
{

    private final VariantContextEnrichment purityEnrichment;
    private final VariantContextEnrichment hotspotEnrichment;
    private final VariantContextEnrichment kataegisEnrichment;
    private final VariantContextEnrichment somaticRefContextEnrichment;
    private final SubclonalLikelihoodEnrichment subclonalLikelihoodEnrichment;
    private final VariantContextEnrichment snpEffEnrichment;

    public SomaticVariantEnrichment(boolean hotspotEnabled, double clonalityBinWidth, @NotNull final String purpleVersion,
            @NotNull final String tumorSample, @NotNull final IndexedFastaSequenceFile reference,
            @NotNull final PurityAdjuster purityAdjuster, @NotNull final DriverGenePanel genePanel,
            @NotNull final List<PurpleCopyNumber> copyNumbers, @NotNull final List<FittedRegion> fittedRegions,
            @NotNull final Multimap<Chromosome, VariantHotspot> hotspots, @NotNull final List<HmfTranscriptRegion> transcripts,
            @NotNull final List<PeakModel> peakModel, @NotNull final Consumer<VariantContext> consumer) {
        subclonalLikelihoodEnrichment = new SubclonalLikelihoodEnrichment(clonalityBinWidth, peakModel, consumer);
        purityEnrichment = new SomaticPurityEnrichment(purpleVersion,
                tumorSample,
                purityAdjuster,
                copyNumbers,
                fittedRegions,
                subclonalLikelihoodEnrichment);
        kataegisEnrichment = new KataegisEnrichment(purityEnrichment);
        somaticRefContextEnrichment = new SomaticRefContextEnrichment(reference, kataegisEnrichment);
        final Set<String> somaticGenes =
                genePanel.driverGenes().stream().filter(DriverGene::reportSomatic).map(DriverGene::gene).collect(Collectors.toSet());
        snpEffEnrichment = new SnpEffEnrichment(somaticGenes, transcripts, somaticRefContextEnrichment);
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
