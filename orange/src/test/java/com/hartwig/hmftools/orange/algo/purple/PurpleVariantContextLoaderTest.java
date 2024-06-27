package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.ImmutableList;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.genotype.GenotypeStatus;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.impact.VariantTranscriptImpact;

import org.junit.Test;

public class PurpleVariantContextLoaderTest
{
    private static final String PURPLE_VARIANT_FILE = Resources.getResource("purple/variants.vcf").getPath();

    @Test
    public void testCanReadVcfFile() throws IOException
    {
        PurpleVariantContextLoader contextLoader = PurpleVariantContextLoader.withPassingOnlyFilter();
        List<PurpleVariantContext> variants = contextLoader.fromVCFFile("tumor_sample", "ref_sample", null, PURPLE_VARIANT_FILE);

        assertEquals(1, variants.size());

        PurpleVariantContext variant = variants.get(0);

        assertEquals(String.valueOf(4), variant.chromosome());
        assertEquals(57181855, variant.position());
        assertEquals(153, variant.allelicDepth().TotalReadCount);
        assertEquals(80, variant.allelicDepth().AlleleReadCount);
        assertEquals(VariantType.SNP, variant.type());
        assertEquals("CRACD", variant.gene());
        assertEquals("C", variant.ref());
        assertEquals("T", variant.alt());
        assertEquals("ENST00000504228", variant.canonicalTranscript());
        assertEquals("synonymous_variant", variant.canonicalEffect());
        assertEquals(CodingEffect.SYNONYMOUS, variant.canonicalCodingEffect());
        assertEquals("c.2187C>T", variant.canonicalHgvsCodingImpact());
        assertEquals("p.Ser729=", variant.canonicalHgvsProteinImpact());
        assertFalse(variant.spliceRegion());
        assertEquals(CodingEffect.SYNONYMOUS, variant.canonicalCodingEffect());

        List<VariantTranscriptImpact> otherImpacts = variant.otherImpacts();
        assertEquals(2, otherImpacts.size());

        // Test if the extra transcripts are loaded properly,
        List<VariantTranscriptImpact> expectedTranscriptImpacts = List.of(
                new VariantTranscriptImpact(
                        "ENSG00000109265",
                        "CRACD",
                        "ENST00000264229",
                        "synonymous_variant",
                        false,
                        "c.2187C>T",
                        "p.Ser729="),
                new VariantTranscriptImpact(
                        "ENSG00000109265",
                        "CRACD",
                        "ENST00000514330",
                        "upstream_gene_variant",
                        true,
                        "",
                        ""
                ));

        for(int i = 0; i < otherImpacts.size(); i++)
        {
            assertVariantTranscriptImpactEquals(expectedTranscriptImpacts.get(i), otherImpacts.get(i));
        }
        // Assertions above already imply it, but explicitly assert the canonical impact isn't included in the otherImpacts field
        assertFalse(otherImpacts.stream().anyMatch(impact -> impact.Transcript.equals(variant.canonicalTranscript())));

        assertEquals(Hotspot.NON_HOTSPOT, variant.hotspot());
        assertFalse(variant.reported());
        assertNull(variant.rnaDepth());
        assertEquals(3.84, variant.adjustedCopyNumber(), 0);
        assertEquals(0.5256, variant.adjustedVAF(), 0);
        assertEquals(1.84, variant.minorAlleleCopyNumber(), 0);
        assertEquals(2.02, variant.variantCopyNumber(), 0);
        assertFalse(variant.biallelic());
        assertEquals(GenotypeStatus.HOM_REF, variant.genotypeStatus());
        assertEquals(2, variant.repeatCount());
        assertEquals(0, variant.subclonalLikelihood(), 0);
        assertEquals(ImmutableList.builder().add(3666).build(), variant.localPhaseSets());
        assertEquals(ImmutableList.builder().add("ENST00000504228").add("ENST00000264229").build(), variant.reportableTranscripts());
    }

    private static void assertVariantTranscriptImpactEquals(VariantTranscriptImpact first, VariantTranscriptImpact second)
    {
        assertEquals(first.GeneId, second.GeneId);
        assertEquals(first.GeneName, second.GeneName);
        assertEquals(first.Transcript, second.Transcript);
        assertEquals(first.Effects, second.Effects);
        assertEquals(first.SpliceRegion, second.SpliceRegion);
        assertEquals(first.HgvsCoding, second.HgvsCoding);
        assertEquals(first.HgvsProtein, second.HgvsProtein);
    }
}