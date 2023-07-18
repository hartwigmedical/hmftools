package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.genotype.GenotypeStatus;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.ImmutableAllelicDepthImpl;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.impact.VariantTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.PurpleCodingEffect;
import com.hartwig.hmftools.datamodel.purple.PurpleGenotypeStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.datamodel.purple.PurpleVariantEffect;
import com.hartwig.hmftools.datamodel.purple.PurpleVariantType;
import com.hartwig.hmftools.orange.algo.pave.PaveAlgo;
import com.hartwig.hmftools.orange.algo.pave.TestEnsemblDataCacheFactory;

import org.junit.Test;

public class PurpleVariantFactoryTest
{

    @Test
    public void testCanCreateVariantFromContext()
    {
        var otherImpacts = List.of(
                new VariantTranscriptImpact(
                        "[ENSG00000109265",
                        "CRACD",
                        "ENST00000264229",
                        "synonymous_variant",
                        false,
                        "c.2187C>T",
                        "p.Ser729=]"
                ));

        PurpleVariantContext context = ImmutablePurpleVariantContext.builder()
                .chromosome("chromosome")
                .position(1)
                .totalReadCount(1)
                .alleleReadCount(1)
                .type(VariantType.INDEL)
                .gene("gene")
                .ref("ref")
                .alt("alt")
                .canonicalTranscript("canonicalTranscript")
                .canonicalEffect("frameshift_variant&missense_variant")
                .canonicalCodingEffect(CodingEffect.NONE)
                .canonicalHgvsCodingImpact("coding")
                .canonicalHgvsProteinImpact("protein")
                .spliceRegion(true)
                .worstCodingEffect(CodingEffect.UNDEFINED)
                .otherImpacts(otherImpacts)
                .hotspot(Hotspot.HOTSPOT)
                .reported(true)
                .tumorDepth(ImmutableAllelicDepthImpl.builder().alleleReadCount(1).totalReadCount(1).build())
                .rnaDepth(null)
                .adjustedCopyNumber(1)
                .adjustedVAF(1)
                .minorAlleleCopyNumber(1)
                .variantCopyNumber(1)
                .biallelic(false)
                .genotypeStatus(GenotypeStatus.UNKNOWN)
                .repeatCount(1)
                .subclonalLikelihood(1)
                .localPhaseSets(List.of(1, 2, 3))
                .build();

        PaveAlgo paveAlgo = new PaveAlgo(TestEnsemblDataCacheFactory.createDummyCache());
        PurpleVariant purpleVariant = new PurpleVariantFactory(paveAlgo).fromPurpleVariantContext(context);

        assertEquals(PurpleVariantType.INDEL, purpleVariant.type());
        assertEquals("gene", purpleVariant.gene());
        assertEquals("chromosome", purpleVariant.chromosome());
        assertEquals(1, purpleVariant.position());
        assertEquals("ref", purpleVariant.ref());
        assertEquals("alt", purpleVariant.alt());
        assertEquals(PurpleCodingEffect.UNDEFINED, purpleVariant.worstCodingEffect());

        var canonicalImpact = purpleVariant.canonicalImpact();
        assertEquals("canonicalTranscript", canonicalImpact.transcript());
        assertEquals("coding", canonicalImpact.hgvsCodingImpact());
        assertEquals("protein", canonicalImpact.hgvsProteinImpact());
        assertNull(canonicalImpact.affectedCodon());
        assertNull(canonicalImpact.affectedExon());
        assertTrue(canonicalImpact.spliceRegion());
        assertEquals(Set.of(PurpleVariantEffect.FRAMESHIFT, PurpleVariantEffect.MISSENSE), canonicalImpact.effects());

        List<PurpleTranscriptImpact> purpleOtherImpacts = purpleVariant.otherImpacts();
        assertEquals(1, purpleOtherImpacts.size());
        var purpleOtherImpact = purpleOtherImpacts.get(0);
        assertEquals("ENST00000264229", purpleOtherImpact.transcript());
        assertEquals("c.2187C>T", purpleOtherImpact.hgvsCodingImpact());
        assertEquals("p.Ser729=", purpleOtherImpact.hgvsProteinImpact());
        assertEquals(Set.of(PurpleVariantEffect.SYNONYMOUS), purpleOtherImpact.effects());
        assertEquals(PurpleCodingEffect.SYNONYMOUS, purpleOtherImpact.codingEffect());
        assertFalse(purpleOtherImpact.spliceRegion());
        assertNull(purpleOtherImpact.affectedCodon());
        assertNull(purpleOtherImpact.affectedExon());

        assertEquals(com.hartwig.hmftools.datamodel.purple.Hotspot.HOTSPOT, purpleVariant.hotspot());
        assertTrue(purpleVariant.reported());
        assertEquals(1, purpleVariant.tumorDepth().alleleReadCount());
        assertEquals(1, purpleVariant.tumorDepth().totalReadCount());
        assertNull(purpleVariant.rnaDepth());
        assertEquals(1, purpleVariant.adjustedCopyNumber(), 0);
        assertEquals(1, purpleVariant.adjustedVAF(), 0);
        assertEquals(1, purpleVariant.minorAlleleCopyNumber(), 0);
        assertEquals(1, purpleVariant.variantCopyNumber(), 0);
        assertFalse(purpleVariant.biallelic());
        assertEquals(PurpleGenotypeStatus.UNKNOWN, purpleVariant.genotypeStatus());
        assertEquals(1, purpleVariant.repeatCount());
        assertEquals(1, purpleVariant.subclonalLikelihood(), 0);
        assertEquals(List.of(1, 2, 3), purpleVariant.localPhaseSets());
    }
}
