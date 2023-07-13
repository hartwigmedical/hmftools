package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleVariantTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.PurpleCodingEffect;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.datamodel.purple.PurpleVariantEffect;
import com.hartwig.hmftools.datamodel.purple.PurpleVariantTranscriptImpact;
import com.hartwig.hmftools.orange.algo.pave.PaveAlgo;
import com.hartwig.hmftools.orange.algo.pave.TestEnsemblDataCacheFactory;
import com.hartwig.hmftools.orange.conversion.PurpleConversion;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PurpleVariantFactoryTest
{

    private static final String SOMATIC_VARIANT_FILE = Resources.getResource("purple/variants.vcf").getPath();

    @NotNull
    private static PurpleVariantFactory createTestFactory()
    {
        return new PurpleVariantFactory(new PaveAlgo(TestEnsemblDataCacheFactory.createDummyCache()));
    }

    @Test
    public void testCanReadVcfFile() throws IOException
    {
        var testFactory = createTestFactory();
        List<PurpleVariant> variants = testFactory.fromVCFFile("COLO829v003T", "COLO829v003R", null, SOMATIC_VARIANT_FILE);

        assertEquals(1, variants.size());

        var variant = variants.get(0);

        assertEquals(String.valueOf(4), variant.chromosome());
        assertEquals(57181855, variant.position());
        assertEquals("C", variant.ref());
        assertEquals("T", variant.alt());
        var canonicalImpact = variant.canonicalImpact();

        assertNotNull(canonicalImpact);
        assertEquals("ENST00000504228", canonicalImpact.transcript());
        assertEquals("c.2187C>T", canonicalImpact.hgvsCodingImpact());
        assertEquals("p.Ser729=", canonicalImpact.hgvsProteinImpact());
        assertEquals(PurpleCodingEffect.SYNONYMOUS, canonicalImpact.codingEffect());

        List<PurpleVariantTranscriptImpact> purpleVariantTranscriptImpacts = variant.variantTranscriptImpacts();
        assertEquals(4, purpleVariantTranscriptImpacts.size());

        // Test if the extra transcripts are loaded properly
        var expectedTranscriptImpacts = List.of(
                ImmutablePurpleVariantTranscriptImpact.builder().
                        transcript("ENST00000264229")
                        .effects(List.of(PurpleVariantEffect.SYNONYMOUS))
                        .hgvsCodingImpact("c.2187C>T")
                        .hgvsProteinImpact("p.Ser729=")
                        .build(),
                ImmutablePurpleVariantTranscriptImpact.builder().
                        transcript("ENST00000504228")
                        .effects(List.of(PurpleVariantEffect.SYNONYMOUS))
                        .hgvsCodingImpact("c.2187C>T")
                        .hgvsProteinImpact("p.Ser729=")
                        .build(),
                ImmutablePurpleVariantTranscriptImpact.builder().
                        transcript("ENST00000514330")
                        .effects(List.of(PurpleVariantEffect.UPSTREAM_GENE))
                        .hgvsCodingImpact("")
                        .hgvsProteinImpact("")
                        .build(),
                ImmutablePurpleVariantTranscriptImpact.builder().
                        transcript("ENST00000541073")
                        .effects(List.of(PurpleVariantEffect.SYNONYMOUS))
                        .hgvsCodingImpact("c.2166C>T")
                        .hgvsProteinImpact("p.Ser722=]")
                        .build());

        assertEquals(expectedTranscriptImpacts, purpleVariantTranscriptImpacts);
    }
}