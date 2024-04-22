package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneTestFactory;
import com.hartwig.hmftools.datamodel.purple.HotspotType;
import com.hartwig.hmftools.datamodel.purple.PurpleCodingEffect;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.datamodel.purple.PurpleVariantType;

import org.junit.Test;

public class SomaticVariantSelectorTest
{
    @Test
    public void canSelectUnreportedNearHotspots()
    {
        PurpleVariant hotspot = TestPurpleVariantFactory.builder().hotspot(HotspotType.HOTSPOT).build();
        PurpleVariant nearHotspot = TestPurpleVariantFactory.builder().hotspot(HotspotType.NEAR_HOTSPOT).build();
        PurpleVariant nonHotspot = TestPurpleVariantFactory.builder().hotspot(HotspotType.NON_HOTSPOT).build();

        List<PurpleVariant> variants =
                SomaticVariantSelector.selectInterestingUnreportedVariants(Lists.newArrayList(hotspot, nearHotspot, nonHotspot),
                        Lists.newArrayList(),
                        Lists.newArrayList());

        assertEquals(2, variants.size());
        assertTrue(variants.contains(hotspot));
        assertTrue(variants.contains(nearHotspot));
    }

    @Test
    public void canSelectVariantsWithReportedPhaseSet()
    {
        PurpleVariant withMatch = TestPurpleVariantFactory.builder().gene("gene").addLocalPhaseSets(1).build();
        PurpleVariant withoutMatch = TestPurpleVariantFactory.builder().gene("gene").addLocalPhaseSets(2).build();
        PurpleVariant withoutPhase = TestPurpleVariantFactory.builder().gene("gene").build();

        PurpleVariant withPhase = TestPurpleVariantFactory.builder().addLocalPhaseSets(1).build();
        PurpleVariant noPhase = TestPurpleVariantFactory.builder().localPhaseSets(null).build();

        List<PurpleVariant> variants =
                SomaticVariantSelector.selectInterestingUnreportedVariants(Lists.newArrayList(withMatch, withoutMatch, withoutPhase),
                        Lists.newArrayList(withPhase, noPhase),
                        Lists.newArrayList());

        assertEquals(1, variants.size());
        assertTrue(variants.contains(withMatch));
    }

    @Test
    public void canSelectVariantsRelevantForCuppa()
    {
        String cuppaGene = SomaticVariantSelector.CUPPA_GENES.iterator().next();
        PurpleVariant cuppaRelevant =
                TestPurpleVariantFactory.builder().gene(cuppaGene).type(PurpleVariantType.INDEL).repeatCount(4).build();

        PurpleVariant cuppaTooManyRepeats =
                TestPurpleVariantFactory.builder().gene(cuppaGene).type(PurpleVariantType.INDEL).repeatCount(10).build();

        PurpleVariant wrongGene =
                TestPurpleVariantFactory.builder().gene("wrong gene").type(PurpleVariantType.INDEL).repeatCount(4).build();

        List<PurpleVariant> variants = SomaticVariantSelector.selectInterestingUnreportedVariants(Lists.newArrayList(cuppaRelevant,
                cuppaTooManyRepeats,
                wrongGene), Lists.newArrayList(), Lists.newArrayList());

        assertEquals(1, variants.size());
        assertTrue(variants.contains(cuppaRelevant));
    }

    @Test
    public void canSelectSynonymousVariantsForDriverGenes()
    {
        String gene = "driver";
        DriverGene driverGene = DriverGeneTestFactory.builder()
                .gene(gene)
                .reportNonsenseAndFrameshift(true)
                .reportSplice(false)
                .reportMissenseAndInframe(true)
                .build();

        PurpleVariant nonsense = TestPurpleVariantFactory.builder()
                .gene(gene)
                .canonicalImpact(TestPurpleVariantFactory.impactBuilder().codingEffect(PurpleCodingEffect.SYNONYMOUS).build())
                .worstCodingEffect(PurpleCodingEffect.NONSENSE_OR_FRAMESHIFT)
                .build();

        PurpleVariant splice = TestPurpleVariantFactory.builder()
                .gene(gene)
                .canonicalImpact(TestPurpleVariantFactory.impactBuilder().codingEffect(PurpleCodingEffect.SYNONYMOUS).build())
                .worstCodingEffect(PurpleCodingEffect.SPLICE)
                .build();

        PurpleVariant missense = TestPurpleVariantFactory.builder()
                .gene(gene)
                .canonicalImpact(TestPurpleVariantFactory.impactBuilder().codingEffect(PurpleCodingEffect.SYNONYMOUS).build())
                .worstCodingEffect(PurpleCodingEffect.MISSENSE)
                .build();

        PurpleVariant wrongGene = TestPurpleVariantFactory.builder()
                .gene("wrong gene")
                .canonicalImpact(TestPurpleVariantFactory.impactBuilder().codingEffect(PurpleCodingEffect.SYNONYMOUS).build())
                .worstCodingEffect(PurpleCodingEffect.MISSENSE)
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
    public void canSelectSpliceRegionVariants()
    {
        DriverGene driverGeneYes = DriverGeneTestFactory.builder().gene("driver 1").reportSplice(true).build();
        DriverGene driverGeneNo = DriverGeneTestFactory.builder().gene("driver 2").reportSplice(false).build();

        PurpleVariant reportedGene = TestPurpleVariantFactory.builder()
                .gene(driverGeneYes.gene())
                .canonicalImpact(TestPurpleVariantFactory.impactBuilder().inSpliceRegion(true).build())
                .build();

        PurpleVariant nonReportedGene = TestPurpleVariantFactory.builder()
                .gene(driverGeneNo.gene())
                .canonicalImpact(TestPurpleVariantFactory.impactBuilder().inSpliceRegion(true).build())
                .build();

        PurpleVariant otherGene = TestPurpleVariantFactory.builder()
                .gene("other gene")
                .canonicalImpact(TestPurpleVariantFactory.impactBuilder().inSpliceRegion(true).build())
                .build();

        List<PurpleVariant> variants =
                SomaticVariantSelector.selectInterestingUnreportedVariants(Lists.newArrayList(reportedGene, nonReportedGene, otherGene),
                        Lists.newArrayList(),
                        Lists.newArrayList(driverGeneYes, driverGeneNo));

        assertEquals(1, variants.size());
        assertTrue(variants.contains(reportedGene));
    }
}