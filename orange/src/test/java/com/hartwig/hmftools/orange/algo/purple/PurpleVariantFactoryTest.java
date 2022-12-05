package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genotype.GenotypeStatus;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.test.SomaticVariantTestFactory;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.ImmutableAllelicDepthImpl;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.impact.AltTranscriptReportableInfo;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PurpleVariantFactoryTest {

    private static final double EPSILON = 1.0E-10;

    @Test
    public void canHandleNullVariants() {
        assertNull(PurpleVariantFactory.create(null));
    }

    @Test
    public void canCreateSimplePurpleVariant() {
        SomaticVariant somaticVariant = SomaticVariantTestFactory.builder()
                .totalReadCount(20)
                .alleleReadCount(10)
                .chromosome("1")
                .position(15)
                .type(VariantType.INDEL)
                .gene("gene")
                .ref("C")
                .alt("CA")
                .canonicalTranscript("canonical transcript")
                .canonicalEffect("missense_variant")
                .canonicalCodingEffect(CodingEffect.MISSENSE)
                .canonicalHgvsCodingImpact("canonical hgvs coding")
                .canonicalHgvsProteinImpact("canonical hgvs protein")
                .qual(1.5)
                .filter("PASS")
                .genesAffected(2)
                .spliceRegion(false)
                .otherReportedEffects(Strings.EMPTY)
                .worstCodingEffect(CodingEffect.NONSENSE_OR_FRAMESHIFT)
                .hotspot(Hotspot.NEAR_HOTSPOT)
                .recovered(true)
                .mappability(0.5)
                .adjustedCopyNumber(0.6)
                .adjustedVAF(0.7)
                .minorAlleleCopyNumber(0.2)
                .variantCopyNumber(2.5)
                .biallelic(false)
                .reported(true)
                .genotypeStatus(GenotypeStatus.HOM_REF)
                .germlineStatus(GermlineStatus.NOISE)
                .trinucleotideContext("TTT")
                .microhomology("GC")
                .repeatSequence("AT")
                .repeatCount(3)
                .kataegis("kataegis")
                .tier(VariantTier.HOTSPOT)
                .subclonalLikelihood(0.35)
                .rnaDepth(createDepth(18, 14))
                .referenceDepth(createDepth(28, 24))
                .localPhaseSets(Lists.newArrayList(1, 2))
                .build();

        List<PurpleVariant> converted = PurpleVariantFactory.create(Lists.newArrayList(somaticVariant));
        assertEquals(1, converted.size());

        PurpleVariant variant = converted.get(0);
        assertEquals(VariantType.INDEL, variant.type());
        assertEquals("gene", variant.gene());
        assertEquals(2, variant.genesAffected());
        assertEquals("1", variant.chromosome());
        assertEquals(15, variant.position());
        assertEquals("C", variant.ref());
        assertEquals("CA", variant.alt());
        assertEquals(CodingEffect.NONSENSE_OR_FRAMESHIFT, variant.worstCodingEffect());
        assertEquals("canonical transcript", variant.canonicalImpact().transcript());
        assertEquals("canonical hgvs coding", variant.canonicalImpact().hgvsCodingImpact());
        assertEquals("canonical hgvs protein", variant.canonicalImpact().hgvsProteinImpact());
        assertNull(variant.canonicalImpact().affectedCodon());
        assertNull(variant.canonicalImpact().affectedExon());
        assertFalse(variant.canonicalImpact().spliceRegion());
        assertEquals(1, variant.canonicalImpact().effects().size());
        assertTrue(variant.canonicalImpact().effects().contains(VariantEffect.MISSENSE));
        assertEquals(CodingEffect.MISSENSE, variant.canonicalImpact().codingEffect());
        assertTrue(variant.otherImpacts().isEmpty());
        assertEquals(Hotspot.NEAR_HOTSPOT, variant.hotspot());
        assertTrue(variant.reported());
        assertFalse(variant.filtered());
        assertEquals("PASS", variant.filter());
        assertTrue(variant.recovered());
        assertEquals(20, variant.tumorDepth().totalReadCount());
        assertEquals(10, variant.tumorDepth().alleleReadCount());
        assertEquals(18, variant.rnaDepth().totalReadCount());
        assertEquals(14, variant.rnaDepth().alleleReadCount());
        assertEquals(28, variant.referenceDepth().totalReadCount());
        assertEquals(24, variant.referenceDepth().alleleReadCount());
        assertEquals(0.6, variant.adjustedCopyNumber(), EPSILON);
        assertEquals(0.7, variant.adjustedVAF(), EPSILON);
        assertEquals(0.2, variant.minorAlleleCopyNumber(), EPSILON);
        assertEquals(2.5, variant.variantCopyNumber(), EPSILON);
        assertFalse(variant.biallelic());
        assertEquals(GenotypeStatus.HOM_REF, variant.genotypeStatus());
        assertEquals(GermlineStatus.NOISE, variant.germlineStatus());
        assertEquals("TTT", variant.trinucleotideContext());
        assertEquals("GC", variant.microhomology());
        assertEquals("AT", variant.repeatSequence());
        assertEquals(3, variant.repeatCount());
        assertEquals("kataegis", variant.kataegis());
        assertEquals(VariantTier.HOTSPOT, variant.tier());
        assertEquals(0.35, variant.subclonalLikelihood(), EPSILON);
        assertEquals(2, variant.localPhaseSets().size());
        assertTrue(variant.localPhaseSets().contains(1));
        assertTrue(variant.localPhaseSets().contains(2));
    }

    @Test
    public void canPopulateOtherTranscriptImpacts() {
        AltTranscriptReportableInfo info =
                new AltTranscriptReportableInfo("transcript", "hgvs coding", "hgvs protein", "missense_variant", CodingEffect.MISSENSE);
        SomaticVariant somaticVariant = SomaticVariantTestFactory.builder().otherReportedEffects(info.serialise()).build();

        PurpleVariant variant = PurpleVariantFactory.create(Lists.newArrayList(somaticVariant)).get(0);

        assertEquals(1, variant.otherImpacts().size());
        PurpleTranscriptImpact impact = variant.otherImpacts().get(0);
        assertEquals("transcript", impact.transcript());
        assertEquals("hgvs coding", impact.hgvsCodingImpact());
        assertEquals("hgvs protein", impact.hgvsProteinImpact());
        assertNull(impact.affectedCodon());
        assertNull(impact.affectedExon());
        assertFalse(impact.spliceRegion());
        assertEquals(1, impact.effects().size());
        assertTrue(impact.effects().contains(VariantEffect.MISSENSE));
        assertEquals(CodingEffect.MISSENSE, impact.codingEffect());
    }

    @NotNull
    private static AllelicDepth createDepth(int totalReadCount, int alleleReadCount) {
        return ImmutableAllelicDepthImpl.builder().totalReadCount(totalReadCount).alleleReadCount(alleleReadCount).build();
    }
}