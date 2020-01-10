package com.hartwig.hmftools.sage.snpeff;

import static java.util.Collections.singletonList;

import static com.hartwig.hmftools.common.variant.VariantConsequence.INFRAME_DELETION;
import static com.hartwig.hmftools.common.variant.VariantConsequence.INFRAME_INSERTION;
import static com.hartwig.hmftools.common.variant.VariantConsequence.SPLICE_DONOR_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.sufferConsequences;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.region.HmfExonRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.sage.vcf.SageVCF;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class SagePostProcessTest {

    private final Map<String, HmfTranscriptRegion> geneMap = HmfGenePanelSupplier.allGenesMap37();
    private final int phase = 1;

    @Test
    public void testInframeIndelStringsProduceExpectedCodingEffect() {
        assertEquals(INFRAME_INSERTION, sufferConsequences(singletonList(SagePostProcess.INFRAME_INSERTION)).get(0));
        assertEquals(INFRAME_DELETION, sufferConsequences(singletonList(SagePostProcess.INFRAME_DELETION)).get(0));

        assertEquals(CodingEffect.MISSENSE, CodingEffect.effect("Gene", singletonList(INFRAME_INSERTION)));
        assertEquals(CodingEffect.MISSENSE, CodingEffect.effect("Gene", singletonList(INFRAME_DELETION)));
    }

    @Test
    public void testSpliceDonorStringsProduceExpectedCodingEffect() {
        assertEquals(SPLICE_DONOR_VARIANT, sufferConsequences(singletonList(SagePostProcess.SPLICE_DONOR_VARIANT)).get(0));
        assertEquals(CodingEffect.SPLICE, CodingEffect.effect("Gene", singletonList(SPLICE_DONOR_VARIANT)));
    }

    @Test
    public void testPhasedInframeDeletes() {

        final HmfTranscriptRegion transcript = geneMap.get("TERT");
        final HmfExonRegion tertExon = transcript.exonByIndex(2);
        final String annotation = snpeffAnnotation("G", "TERT", transcript.transcriptID());
        final VariantContext first = create(tertExon.chromosome(), tertExon.start(), "GATAC", "G", phase, annotation);
        final VariantContext second = create(tertExon.chromosome(), tertExon.start() + 10, "GAT", "G", phase, annotation);

        processAndCheckPhased(true, first, second);
    }

    @Test
    public void testPhasedInframeInserts() {

        final HmfTranscriptRegion transcript = geneMap.get("TERT");
        final HmfExonRegion tertExon = transcript.exonByIndex(2);
        final String annotation = snpeffAnnotation("G", "TERT", transcript.transcriptID());
        final VariantContext first = create(tertExon.chromosome(), tertExon.start(), "G", "GATAC", phase, annotation);
        final VariantContext second = create(tertExon.chromosome(), tertExon.start() + 10, "G", "GAT", phase, annotation);

        processAndCheckPhased(true, first, second);
    }

    @Test
    public void testPhasedInframeInsertAndDelete() {

        final HmfTranscriptRegion transcript = geneMap.get("TERT");
        final HmfExonRegion tertExon = transcript.exonByIndex(2);
        final String annotation = snpeffAnnotation("G", "TERT", transcript.transcriptID());
        final VariantContext first = create(tertExon.chromosome(), tertExon.start(), "G", "GATAC", phase, annotation);
        final VariantContext second = create(tertExon.chromosome(), tertExon.start() + 10, "GA", "G", phase, annotation);

        processAndCheckPhased(true, first, second);
    }

    @Test
    public void testInframeInsertAndDeleteAreNotPhased() {
        final HmfTranscriptRegion transcript = geneMap.get("TERT");
        final HmfExonRegion tertExon = transcript.exonByIndex(2);
        final String annotation = snpeffAnnotation("G", "TERT", transcript.transcriptID());
        final VariantContext first = create(tertExon.chromosome(), tertExon.start(), "G", "GATAC", phase, annotation);
        final VariantContext second = create(tertExon.chromosome(), tertExon.start() + 10, "GA", "G", phase + 1, annotation);

        processAndCheckPhased(false, first, second);
    }

    @Test
    public void testInsertAndDeleteDoNotProduceInframe() {

        final HmfTranscriptRegion transcript = geneMap.get("TERT");
        final HmfExonRegion tertExon = transcript.exonByIndex(2);
        final String annotation = snpeffAnnotation("G", "TERT", transcript.transcriptID());
        final VariantContext first = create(tertExon.chromosome(), tertExon.start(), "G", "GATAC", phase, annotation);
        final VariantContext second = create(tertExon.chromosome(), tertExon.start() + 10, "GAT", "G", phase, annotation);

        processAndCheckPhased(false, first, second);
    }

    private void processAndCheckPhased(boolean expectedFlag, VariantContext... contexts) {
        final List<VariantContext> catcher = Lists.newArrayList();
        final SagePostProcess postProcess = new SagePostProcess(geneMap, catcher::add);

        for (VariantContext context : contexts) {
            postProcess.accept(context);
        }

        postProcess.close();

        for (int i = 0; i < contexts.length; i++) {
            assertEquals(contexts[i], catcher.get(i));
            assertEquals(expectedFlag, catcher.get(i).hasAttribute(SagePostProcessVCF.PHASED_INFRAME_INDEL));
        }
    }

    @NotNull
    private VariantContext create(final String chr, long position, final String ref, @NotNull final String alt, int lps,
            String... annotations) {

        Allele refAllele = Allele.create(ref, true);
        Allele altAllele = Allele.create(alt, false);
        List<Allele> alleles = Lists.newArrayList(refAllele, altAllele);

        Genotype refGenotype = new GenotypeBuilder("NORMAL").DP(30).AD(new int[] { 1, 29 }).alleles(alleles).make();
        Genotype altGenotype = new GenotypeBuilder("TUMOR").DP(90).AD(new int[] { 20, 70 }).alleles(alleles).make();

        VariantContextBuilder builder = new VariantContextBuilder().chr(chr)
                .start(position)
                .genotypes(refGenotype, altGenotype)
                .alleles(alleles)
                .computeEndFromAlleles(alleles, (int) position);

        if (lps > 0) {
            builder.attribute(SageVCF.PHASE, lps);
        }

        builder.attribute("ANN", Lists.newArrayList(annotations));

        return builder.make();
    }

    @NotNull
    private String snpeffAnnotation(String alt, String gene, String transcript) {
        return alt + "|upstream_gene_variant|MODIFIER|" + gene + "|ENSG00000164362|transcript|" + transcript
                + "|protein_coding||c.-139_-138delCCinsTT|||||80|";
    }

}
