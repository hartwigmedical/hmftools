package com.hartwig.hmftools.sage.snpeff;

import static java.util.Collections.singletonList;

import static com.hartwig.hmftools.common.sage.SagePostProcessVCF.HMF_CANONICAL_EFFECT;
import static com.hartwig.hmftools.common.sage.SagePostProcessVCF.HMF_CANONICAL_IMPACT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.INFRAME_DELETION;
import static com.hartwig.hmftools.common.variant.VariantConsequence.INFRAME_INSERTION;
import static com.hartwig.hmftools.common.variant.VariantConsequence.SPLICE_DONOR_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.sufferConsequences;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.region.CanonicalTranscript;
import com.hartwig.hmftools.common.genome.region.CanonicalTranscriptFactory;
import com.hartwig.hmftools.common.genome.region.HmfExonRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.enrich.SomaticRefContextEnrichment;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationFactory;
import com.hartwig.hmftools.sage.vcf.SageVCF;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class SagePostProcessTest {

    private final List<CanonicalTranscript> transcripts = CanonicalTranscriptFactory.create37();
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

        postProcess(first, second);
        assertPhasedIndel(first);
        assertPhasedIndel(second);
    }

    @Test
    public void testPhasedInframeInserts() {

        final HmfTranscriptRegion transcript = geneMap.get("TERT");
        final HmfExonRegion tertExon = transcript.exonByIndex(2);
        final String annotation = snpeffAnnotation("G", "TERT", transcript.transcriptID());
        final VariantContext first = create(tertExon.chromosome(), tertExon.start(), "G", "GATAC", phase, annotation);
        final VariantContext second = create(tertExon.chromosome(), tertExon.start() + 10, "G", "GAT", phase, annotation);

        postProcess(first, second);
        assertPhasedIndel(first);
        assertPhasedIndel(second);

    }

    @Test
    public void testPhasedInframeInsertAndDelete() {

        final HmfTranscriptRegion transcript = geneMap.get("TERT");
        final HmfExonRegion tertExon = transcript.exonByIndex(2);
        final String annotation = snpeffAnnotation("G", "TERT", transcript.transcriptID());
        final VariantContext first = create(tertExon.chromosome(), tertExon.start(), "G", "GATAC", phase, annotation);
        final VariantContext second = create(tertExon.chromosome(), tertExon.start() + 10, "GA", "G", phase, annotation);

        postProcess(first, second);
        assertPhasedIndel(first);
        assertPhasedIndel(second);
    }

    @Test
    public void testInframeInsertAndDeleteAreNotPhased() {
        final HmfTranscriptRegion transcript = geneMap.get("TERT");
        final HmfExonRegion tertExon = transcript.exonByIndex(2);
        final String annotation = snpeffAnnotation("G", "TERT", transcript.transcriptID());
        final VariantContext first = create(tertExon.chromosome(), tertExon.start(), "G", "GATAC", phase, annotation);
        final VariantContext second = create(tertExon.chromosome(), tertExon.start() + 10, "GA", "G", phase + 1, annotation);

        postProcess(first, second);
        assertStuff(CodingEffect.NONE, "upstream_gene_variant", first);
        assertStuff(CodingEffect.NONE, "upstream_gene_variant", second);
    }

    @Test
    public void testInsertAndDeleteDoNotProduceInframe() {

        final HmfTranscriptRegion transcript = geneMap.get("TERT");
        final HmfExonRegion tertExon = transcript.exonByIndex(2);
        final String annotation = snpeffAnnotation("G", "TERT", transcript.transcriptID());
        final VariantContext first = create(tertExon.chromosome(), tertExon.start(), "G", "GATAC", phase, annotation);
        final VariantContext second = create(tertExon.chromosome(), tertExon.start() + 10, "GAT", "G", phase, annotation);

        postProcess(first, second);
        assertStuff(CodingEffect.NONE, "upstream_gene_variant", first);
        assertStuff(CodingEffect.NONE, "upstream_gene_variant", second);
    }

    @Test
    public void testInframeDueToMicrohomology() {
        final String annotation =
                "C|frameshift_variant&splice_acceptor_variant&splice_region_variant&intron_variant|HIGH|KIT|ENSG00000157404|transcript|ENST00000288135|protein_coding|11/21|c.1648-2_1672delAGAAACCCATGTATGAAGTACAGTGGA|p.Lys550fs||1648/2931|550/976||";
        final VariantContext original = create("4", 55593579, "CAGAAACCCATGTATGAAGTACAGTGGA", "C", "AG", 0, annotation);

        final SnpEffAnnotation originalAnnotation = SnpEffAnnotationFactory.fromContext(original).get(0);
        final CodingEffect originalEffect = CodingEffect.effect("KIT", originalAnnotation.consequences());
        assertEquals(CodingEffect.NONSENSE_OR_FRAMESHIFT, originalEffect);

        postProcess(original);

        assertStuff(CodingEffect.MISSENSE, "inframe_mh", original);
    }

    private void postProcess(VariantContext... contexts) {
        final List<VariantContext> catcher = Lists.newArrayList();
        final SagePostProcess postProcess = new SagePostProcess(geneMap, transcripts, catcher::add);

        for (VariantContext context : contexts) {
            postProcess.accept(context);
        }

        postProcess.close();
    }

    private void assertPhasedIndel(VariantContext context) {
        assertStuff(CodingEffect.MISSENSE, "inframe_phased", context);
    }

    private void assertStuff(CodingEffect expectedEffect, String expectedImpact, VariantContext context) {
        assertEquals(expectedImpact, context.getAttributeAsStringList(HMF_CANONICAL_EFFECT, "").get(0));
        assertEquals(expectedEffect.toString(), context.getAttribute(HMF_CANONICAL_IMPACT));
    }

    @NotNull
    private VariantContext create(final String chr, long position, final String ref, @NotNull final String alt, int lps,
            String... annotations) {
        return create(chr, position, ref, alt, Strings.EMPTY, lps, annotations);
    }

    @NotNull
    private VariantContext create(final String chr, long position, final String ref, @NotNull final String alt, @NotNull final String mh,
            int lps, String... annotations) {

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

        if (!mh.isEmpty()) {
            builder.attribute(SomaticRefContextEnrichment.MICROHOMOLOGY_FLAG, mh);
        }

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
