package com.hartwig.hmftools.common.variant.enrich;

import java.util.Collections;
import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.genome.region.CanonicalTranscript;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class GermlineVariantEnrichment implements VariantContextEnrichment {

    private final VariantContextEnrichment purityEnrichment;
    private final VariantContextEnrichment refGenomeEnrichment;
    private final VariantContextEnrichment clinvarEnrichment;
    private final VariantContextEnrichment snpEffEnrichment;

    public GermlineVariantEnrichment(@NotNull final String purpleVersion, @NotNull final String tumorSample,
            @NotNull final IndexedFastaSequenceFile reference, @NotNull final PurityAdjuster purityAdjuster,
            @NotNull final List<PurpleCopyNumber> copyNumbers, @NotNull final DriverGenePanel genePanel,
            @NotNull final List<CanonicalTranscript> transcripts, @NotNull final Consumer<VariantContext> consumer) {
        this.clinvarEnrichment = new ClinvarEnrichment(consumer);
        this.refGenomeEnrichment = new SomaticRefContextEnrichment(reference, clinvarEnrichment);
        this.snpEffEnrichment = new SnpEffEnrichment(genePanel.driverGenes(), transcripts, refGenomeEnrichment);


        this.purityEnrichment =
                new PurityEnrichment(purpleVersion, tumorSample, purityAdjuster, copyNumbers, Collections.emptyList(), snpEffEnrichment);
    }

    @Override
    public void accept(@NotNull final VariantContext context) {
        purityEnrichment.accept(context);
    }

    @Override
    public void flush() {
        purityEnrichment.flush();
        snpEffEnrichment.flush();
        refGenomeEnrichment.flush();
        clinvarEnrichment.flush();
    }

    @NotNull
    @Override
    public VCFHeader enrichHeader(@NotNull final VCFHeader template) {
        VCFHeader header = purityEnrichment.enrichHeader(template);
        header = refGenomeEnrichment.enrichHeader(header);
        header = snpEffEnrichment.enrichHeader(header);
        return clinvarEnrichment.enrichHeader(header);
    }
}
