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
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.impact.VariantTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.HotspotType;
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
        List<VariantTranscriptImpact> otherImpacts = List.of(
                new VariantTranscriptImpact(
                        "ENSG00000109265",
                        "CRACD",
                        "ENST00000264229",
                        "synonymous_variant",
                        false,
                        "c.2187C>T",
                        "p.Ser729="
                ));

        PurpleVariantContext context = ImmutablePurpleVariantContext.builder()
                .chromosome("chromosome")
                .position(57181855)
                .totalReadCount(153)
                .alleleReadCount(80)
                .type(VariantType.INDEL)
                .gene("CRACD")
                .ref("C")
                .alt("T")
                .canonicalTranscript("ENST00000504228")
                .canonicalEffect("frameshift_variant&missense_variant")
                .canonicalCodingEffect(CodingEffect.NONE)
                .canonicalHgvsCodingImpact("c.2187C>T")
                .canonicalHgvsProteinImpact("p.Ser729=")
                .spliceRegion(true)
                .worstCodingEffect(CodingEffect.UNDEFINED)
                .otherImpacts(otherImpacts)
                .hotspot(Hotspot.HOTSPOT)
                .reported(true)
                .rnaDepth(null)
                .adjustedCopyNumber(3.84)
                .adjustedVAF(0.5256)
                .minorAlleleCopyNumber(1.84)
                .variantCopyNumber(2.02)
                .biallelic(false)
                .genotypeStatus(GenotypeStatus.UNKNOWN)
                .repeatCount(2)
                .subclonalLikelihood(1)
                .localPhaseSets(List.of(1, 2, 3))
                .build();

        PaveAlgo paveAlgo = new PaveAlgo(TestEnsemblDataCacheFactory.createDummyCache(), false);
        PurpleVariant purpleVariant = new PurpleVariantFactory(paveAlgo).fromPurpleVariantContext(context);

        assertEquals(PurpleVariantType.INDEL, purpleVariant.type());
        assertEquals("CRACD", purpleVariant.gene());
        assertEquals("chromosome", purpleVariant.chromosome());
        assertEquals(57181855, purpleVariant.position());
        assertEquals("C", purpleVariant.ref());
        assertEquals("T", purpleVariant.alt());
        assertEquals(PurpleCodingEffect.UNDEFINED, purpleVariant.worstCodingEffect());

        PurpleTranscriptImpact canonicalImpact = purpleVariant.canonicalImpact();
        assertEquals("ENST00000504228", canonicalImpact.transcript());
        assertEquals("c.2187C>T", canonicalImpact.hgvsCodingImpact());
        assertEquals("p.Ser729=", canonicalImpact.hgvsProteinImpact());
        assertNull(canonicalImpact.affectedCodon());
        assertNull(canonicalImpact.affectedExon());
        assertTrue(canonicalImpact.inSpliceRegion());
        assertEquals(Set.of(PurpleVariantEffect.FRAMESHIFT, PurpleVariantEffect.MISSENSE), canonicalImpact.effects());

        List<PurpleTranscriptImpact> purpleOtherImpacts = purpleVariant.otherImpacts();
        assertEquals(1, purpleOtherImpacts.size());
        PurpleTranscriptImpact purpleOtherImpact = purpleOtherImpacts.get(0);
        assertEquals("ENST00000264229", purpleOtherImpact.transcript());
        assertEquals("c.2187C>T", purpleOtherImpact.hgvsCodingImpact());
        assertEquals("p.Ser729=", purpleOtherImpact.hgvsProteinImpact());
        assertEquals(Set.of(PurpleVariantEffect.SYNONYMOUS), purpleOtherImpact.effects());
        assertEquals(PurpleCodingEffect.SYNONYMOUS, purpleOtherImpact.codingEffect());
        assertFalse(purpleOtherImpact.inSpliceRegion());
        assertNull(purpleOtherImpact.affectedCodon());
        assertNull(purpleOtherImpact.affectedExon());

        assertEquals(HotspotType.HOTSPOT, purpleVariant.hotspot());
        assertTrue(purpleVariant.reported());
        assertEquals(80, purpleVariant.tumorDepth().alleleReadCount());
        assertEquals(153, purpleVariant.tumorDepth().totalReadCount());
        assertNull(purpleVariant.rnaDepth());
        assertEquals(3.84, purpleVariant.adjustedCopyNumber(), 0);
        assertEquals(0.5256, purpleVariant.adjustedVAF(), 0);
        assertEquals(1.84, purpleVariant.minorAlleleCopyNumber(), 0);
        assertEquals(2.02, purpleVariant.variantCopyNumber(), 0);
        assertFalse(purpleVariant.biallelic());
        assertEquals(PurpleGenotypeStatus.UNKNOWN, purpleVariant.genotypeStatus());
        assertEquals(2, purpleVariant.repeatCount());
        assertEquals(1, purpleVariant.subclonalLikelihood(), 0);
        assertEquals(List.of(1, 2, 3), purpleVariant.localPhaseSets());
    }
}
