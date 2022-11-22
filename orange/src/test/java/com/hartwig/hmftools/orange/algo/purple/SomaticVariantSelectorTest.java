package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneTestFactory;
import com.hartwig.hmftools.common.test.SomaticVariantTestFactory;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.common.variant.ReportableVariantTestFactory;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class SomaticVariantSelectorTest {

    @Test
    public void canSelectUnreportedNearHotspots() {
        SomaticVariant hotspot =
                SomaticVariantTestFactory.builder().canonicalEffect("hotspot").hotspot(Hotspot.HOTSPOT).reported(false).build();
        SomaticVariant nearHotspot =
                SomaticVariantTestFactory.builder().canonicalEffect("near hotspot").hotspot(Hotspot.NEAR_HOTSPOT).reported(false).build();
        SomaticVariant nonHotspot =
                SomaticVariantTestFactory.builder().canonicalEffect("non hotspot").hotspot(Hotspot.NON_HOTSPOT).reported(false).build();

        List<ReportableVariant> variants =
                SomaticVariantSelector.selectInterestingUnreportedVariants(Lists.newArrayList(hotspot, nearHotspot, nonHotspot),
                        Lists.newArrayList(),
                        Lists.newArrayList());

        assertEquals(2, variants.size());
        assertNotNull(findByCanonicalEffect(variants, "hotspot"));
        assertNotNull(findByCanonicalEffect(variants, "near hotspot"));
    }

    @Test
    public void canSelectVariantsWithReportedPhaseSet() {
        SomaticVariant withMatch =
                SomaticVariantTestFactory.builder().gene("gene").canonicalEffect("match").reported(false).addLocalPhaseSets(1).build();
        SomaticVariant withoutMatch =
                SomaticVariantTestFactory.builder().gene("gene").canonicalEffect("no match").reported(false).addLocalPhaseSets(2).build();
        SomaticVariant withoutPhase = SomaticVariantTestFactory.builder().gene("gene").canonicalEffect("no phase").reported(false).build();

        ReportableVariant withPhase = ReportableVariantTestFactory.builder().localPhaseSet(1).build();
        ReportableVariant noPhase = ReportableVariantTestFactory.builder().localPhaseSet(null).build();

        List<ReportableVariant> variants =
                SomaticVariantSelector.selectInterestingUnreportedVariants(Lists.newArrayList(withMatch, withoutMatch, withoutPhase),
                        Lists.newArrayList(withPhase, noPhase),
                        Lists.newArrayList());

        assertEquals(1, variants.size());
        assertNotNull(findByCanonicalEffect(variants, "match"));
    }

    @Test
    public void canSelectVariantsRelevantForCuppa() {
        String cuppaGene = SomaticVariantSelector.CUPPA_GENES.iterator().next();
        SomaticVariant cuppaRelevant = SomaticVariantTestFactory.builder()
                .canonicalEffect("relevant")
                .gene(cuppaGene)
                .type(VariantType.INDEL)
                .repeatCount(4)
                .reported(false)
                .build();

        SomaticVariant cuppaTooManyRepeats = SomaticVariantTestFactory.builder()
                .canonicalEffect("too many repeats")
                .gene(cuppaGene)
                .type(VariantType.INDEL)
                .repeatCount(10)
                .reported(false)
                .build();

        SomaticVariant wrongGene = SomaticVariantTestFactory.builder()
                .canonicalEffect("wrong")
                .gene("wrong gene")
                .type(VariantType.INDEL)
                .repeatCount(4)
                .reported(false)
                .build();

        List<ReportableVariant> variants = SomaticVariantSelector.selectInterestingUnreportedVariants(Lists.newArrayList(cuppaRelevant,
                cuppaTooManyRepeats,
                wrongGene), Lists.newArrayList(), Lists.newArrayList());

        assertEquals(1, variants.size());
        assertNotNull(findByCanonicalEffect(variants, "relevant"));
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

        SomaticVariant nonsense = SomaticVariantTestFactory.builder()
                .gene(gene)
                .canonicalEffect("nonsense")
                .canonicalCodingEffect(CodingEffect.SYNONYMOUS)
                .worstCodingEffect(CodingEffect.NONSENSE_OR_FRAMESHIFT)
                .build();

        SomaticVariant splice = SomaticVariantTestFactory.builder()
                .gene(gene)
                .canonicalEffect("splice")
                .canonicalCodingEffect(CodingEffect.SYNONYMOUS)
                .worstCodingEffect(CodingEffect.SPLICE)
                .build();

        SomaticVariant missense = SomaticVariantTestFactory.builder()
                .gene(gene)
                .canonicalEffect("missense")
                .canonicalCodingEffect(CodingEffect.SYNONYMOUS)
                .worstCodingEffect(CodingEffect.MISSENSE)
                .build();

        SomaticVariant wrongGene = SomaticVariantTestFactory.builder()
                .gene("wrong gene")
                .canonicalEffect("wrong gene")
                .canonicalCodingEffect(CodingEffect.SYNONYMOUS)
                .worstCodingEffect(CodingEffect.MISSENSE)
                .build();

        List<ReportableVariant> variants =
                SomaticVariantSelector.selectInterestingUnreportedVariants(Lists.newArrayList(nonsense, splice, missense, wrongGene),
                        Lists.newArrayList(),
                        Lists.newArrayList(driverGene));

        assertEquals(2, variants.size());
        assertNotNull(findByCanonicalEffect(variants, "nonsense"));
        assertNotNull(findByCanonicalEffect(variants, "missense"));
    }

    @Test
    public void canSelectSpliceRegionVariants() {
        DriverGene driverGeneYes = DriverGeneTestFactory.builder().gene("driver 1").reportSplice(true).build();
        DriverGene driverGeneNo = DriverGeneTestFactory.builder().gene("driver 2").reportSplice(false).build();

        SomaticVariant reportedGene =
                SomaticVariantTestFactory.builder().spliceRegion(true).canonicalEffect("all correct").gene(driverGeneYes.gene()).build();

        SomaticVariant nonReportedGene = SomaticVariantTestFactory.builder()
                .spliceRegion(true)
                .canonicalEffect("non reported gene")
                .gene(driverGeneNo.gene())
                .build();

        SomaticVariant otherGene =
                SomaticVariantTestFactory.builder().spliceRegion(true).canonicalEffect("other gene").gene(driverGeneNo.gene()).build();

        List<ReportableVariant> variants =
                SomaticVariantSelector.selectInterestingUnreportedVariants(Lists.newArrayList(reportedGene, nonReportedGene, otherGene),
                        Lists.newArrayList(),
                        Lists.newArrayList(driverGeneYes, driverGeneNo));

        assertEquals(1, variants.size());
        assertNotNull(findByCanonicalEffect(variants, "all correct"));
    }

    @Nullable
    private static ReportableVariant findByCanonicalEffect(@NotNull List<ReportableVariant> variants,
            @NotNull String canonicalEffectToFind) {
        for (ReportableVariant variant : variants) {
            if (variant.canonicalEffect().equals(canonicalEffectToFind)) {
                return variant;
            }
        }
        return null;
    }
}