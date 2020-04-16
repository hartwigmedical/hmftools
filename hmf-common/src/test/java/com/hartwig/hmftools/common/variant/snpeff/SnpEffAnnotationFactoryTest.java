package com.hartwig.hmftools.common.variant.snpeff;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.CodingEffect;

import org.junit.Test;

public class SnpEffAnnotationFactoryTest {

    @Test
    public void testSpliceDonorPlus5() {
        final String annotation =
                "A|splice_region_variant&intron_variant|LOW|TP53|ENSG00000141510|transcript|ENST00000269305|protein_coding|4/10|c.375+5G>T||||||,"
                        + "A|splice_region_variant&intron_variant|LOW|TP53|ENSG00000141510|transcript|ENST00000413465|protein_coding|3/6|c.375+5G>T||||||,"
                        + "A|splice_region_variant&intron_variant|LOW|TP53|ENSG00000141510|transcript|ENST00000359597|protein_coding|3/8|c.375+5G>T||||||,"
                        + "A|splice_region_variant&intron_variant|LOW|TP53|ENSG00000141510|transcript|ENST00000420246|protein_coding|4/11|c.375+5G>T||||||,"
                        + "A|splice_region_variant&intron_variant|LOW|TP53|ENSG00000141510|transcript|ENST00000455263|protein_coding|4/11|c.375+5G>T||||||,"
                        + "A|splice_region_variant&intron_variant|LOW|TP53|ENSG00000141510|transcript|ENST00000445888|protein_coding|4/10|c.375+5G>T||||||,"
                        + "A|splice_region_variant&intron_variant|LOW|TP53|ENSG00000141510|transcript|ENST00000505014|retained_intron|3/4|n.631+5G>T||||||,"
                        + "A|splice_region_variant&intron_variant|LOW|TP53|ENSG00000141510|transcript|ENST00000508793|protein_coding|4/4|c.375+5G>T||||||WARNING_TRANSCRIPT_INCOMPLETE,"
                        + "A|splice_region_variant&intron_variant|LOW|TP53|ENSG00000141510|transcript|ENST00000604348|protein_coding|4/4|c.375+5G>T||||||WARNING_TRANSCRIPT_NO_STOP_CODON,"
                        + "A|splice_region_variant&intron_variant|LOW|TP53|ENSG00000141510|transcript|ENST00000503591|protein_coding|5/5|c.375+5G>T||||||WARNING_TRANSCRIPT_INCOMPLETE,"
                        + "A|upstream_gene_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000504290|retained_intron||n.-496G>T|||||496|,"
                        + "A|upstream_gene_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000504937|retained_intron||n.-496G>T|||||496|,"
                        + "A|upstream_gene_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000510385|retained_intron||n.-496G>T|||||496|,"
                        + "A|upstream_gene_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000574684|processed_transcript||n.-870G>T|||||870|,"
                        + "A|intron_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000509690|protein_coding|1/5|c.-21-753G>T||||||WARNING_TRANSCRIPT_NO_STOP_CODON,"
                        + "A|intron_variant|MODIFIER|TP53|ENSG00000141510|transcript|ENST00000514944|protein_coding|3/5|c.96+393G>T||||||WARNING_TRANSCRIPT_INCOMPLETE";

        final List<String> annotationList = Lists.newArrayList(annotation.split(","));
        final List<SnpEffAnnotation> annotations = SnpEffAnnotationFactory.fromAnnotationList(annotationList);
        final SnpEffAnnotation first = annotations.get(0);
        final CodingEffect codingEffect = CodingEffect.effect("RANDOM", first.consequences());
        assertEquals(CodingEffect.SPLICE, codingEffect);
    }

}
