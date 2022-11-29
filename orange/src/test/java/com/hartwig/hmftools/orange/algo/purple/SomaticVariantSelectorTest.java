package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneTestFactory;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.VariantType;

import org.junit.Test;

public class SomaticVariantSelectorTest {

    @Test
    public void canSelectUnreportedNearHotspots() {
        PurpleVariant hotspot = PurpleVariantTestFactory.builder().hotspot(Hotspot.HOTSPOT).reported(false).build();
        PurpleVariant nearHotspot = PurpleVariantTestFactory.builder().hotspot(Hotspot.NEAR_HOTSPOT).reported(false).build();
        PurpleVariant nonHotspot = PurpleVariantTestFactory.builder().hotspot(Hotspot.NON_HOTSPOT).reported(false).build();

        List<PurpleVariant> variants =
                SomaticVariantSelector.selectInterestingUnreportedVariants(Lists.newArrayList(hotspot, nearHotspot, nonHotspot),
                        Lists.newArrayList(),
                        Lists.newArrayList());

        assertEquals(2, variants.size());
        assertTrue(variants.contains(hotspot));
        assertTrue(variants.contains(nearHotspot));
    }

    @Test
    public void canSelectVariantsWithReportedPhaseSet() {
        PurpleVariant withMatch = PurpleVariantTestFactory.builder().gene("gene").reported(false).addLocalPhaseSets(1).build();
        PurpleVariant withoutMatch = PurpleVariantTestFactory.builder().gene("gene").reported(false).addLocalPhaseSets(2).build();
        PurpleVariant withoutPhase = PurpleVariantTestFactory.builder().gene("gene").reported(false).build();

        PurpleVariant withPhase = PurpleVariantTestFactory.builder().addLocalPhaseSets(1).build();
        PurpleVariant noPhase = PurpleVariantTestFactory.builder().localPhaseSets(null).build();

        List<PurpleVariant> variants =
                SomaticVariantSelector.selectInterestingUnreportedVariants(Lists.newArrayList(withMatch, withoutMatch, withoutPhase),
                        Lists.newArrayList(withPhase, noPhase),
                        Lists.newArrayList());

        assertEquals(1, variants.size());
        assertTrue(variants.contains(withMatch));
    }

    @Test
    public void canSelectVariantsRelevantForCuppa() {
        String cuppaGene = SomaticVariantSelector.CUPPA_GENES.iterator().next();
        PurpleVariant cuppaRelevant =
                PurpleVariantTestFactory.builder().gene(cuppaGene).type(VariantType.INDEL).repeatCount(4).reported(false).build();

        PurpleVariant cuppaTooManyRepeats =
                PurpleVariantTestFactory.builder().gene(cuppaGene).type(VariantType.INDEL).repeatCount(10).reported(false).build();

        PurpleVariant wrongGene =
                PurpleVariantTestFactory.builder().gene("wrong gene").type(VariantType.INDEL).repeatCount(4).reported(false).build();

        List<PurpleVariant> variants = SomaticVariantSelector.selectInterestingUnreportedVariants(Lists.newArrayList(cuppaRelevant,
                cuppaTooManyRepeats,
                wrongGene), Lists.newArrayList(), Lists.newArrayList());

        assertEquals(1, variants.size());
        assertTrue(variants.contains(cuppaRelevant));
    }

    @Test
    public void canSelectSynonymousVariantsForDriverGenes() {
        String gene = "driver";
        DriverGene driverGene = DriverGeneTestFactory.builder()
                .gene(gene)
                .reportNonsenseAndFrameshift(true)
                .reportSplice(false)
                .reportMissenseAndInframe(true)
                .build();

        PurpleVariant nonsense = PurpleVariantTestFactory.builder()
                .gene(gene)
                .canonicalImpact(PurpleVariantTestFactory.impactBuilder().codingEffect(CodingEffect.SYNONYMOUS).build())
                .worstCodingEffect(CodingEffect.NONSENSE_OR_FRAMESHIFT)
                .build();

        PurpleVariant splice = PurpleVariantTestFactory.builder()
                .gene(gene)
                .canonicalImpact(PurpleVariantTestFactory.impactBuilder().codingEffect(CodingEffect.SYNONYMOUS).build())
                .worstCodingEffect(CodingEffect.SPLICE)
                .build();

        PurpleVariant missense = PurpleVariantTestFactory.builder()
                .gene(gene)
                .canonicalImpact(PurpleVariantTestFactory.impactBuilder().codingEffect(CodingEffect.SYNONYMOUS).build())
                .worstCodingEffect(CodingEffect.MISSENSE)
                .build();

        PurpleVariant wrongGene = PurpleVariantTestFactory.builder()
                .gene("wrong gene")
                .canonicalImpact(PurpleVariantTestFactory.impactBuilder().codingEffect(CodingEffect.SYNONYMOUS).build())
                .worstCodingEffect(CodingEffect.MISSENSE)
                .build();

        List<PurpleVariant> variants =
                SomaticVariantSelector.selectInterestingUnreportedVariants(Lists.newArrayList(nonsense, splice, missense, wrongGene),
                        Lists.newArrayList(),
                        Lists.newArrayList(driverGene));

        assertEquals(2, variants.size());
        assertTrue(variants.contains(nonsense));
        assertTrue(variants.contains(missense));
    }

    @Test
    public void canSelectSpliceRegionVariants() {
        DriverGene driverGeneYes = DriverGeneTestFactory.builder().gene("driver 1").reportSplice(true).build();
        DriverGene driverGeneNo = DriverGeneTestFactory.builder().gene("driver 2").reportSplice(false).build();

        PurpleVariant reportedGene = PurpleVariantTestFactory.builder()
                .gene(driverGeneYes.gene())
                .canonicalImpact(PurpleVariantTestFactory.impactBuilder().spliceRegion(true).build())
                .build();

        PurpleVariant nonReportedGene = PurpleVariantTestFactory.builder()
                .gene(driverGeneNo.gene())
                .canonicalImpact(PurpleVariantTestFactory.impactBuilder().spliceRegion(true).build())
                .build();

        PurpleVariant otherGene = PurpleVariantTestFactory.builder()
                .gene("other gene")
                .canonicalImpact(PurpleVariantTestFactory.impactBuilder().spliceRegion(true).build())
                .build();

        List<PurpleVariant> variants =
                SomaticVariantSelector.selectInterestingUnreportedVariants(Lists.newArrayList(reportedGene, nonReportedGene, otherGene),
                        Lists.newArrayList(),
                        Lists.newArrayList(driverGeneYes, driverGeneNo));

        assertEquals(1, variants.size());
        assertTrue(variants.contains(reportedGene));
    }
}