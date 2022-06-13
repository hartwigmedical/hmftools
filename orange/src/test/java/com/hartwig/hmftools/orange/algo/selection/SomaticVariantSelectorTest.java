package com.hartwig.hmftools.orange.algo.selection;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneTestFactory;
import com.hartwig.hmftools.common.protect.ProtectEventGenerator;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.protect.ProtectTestFactory;
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
    public void canSelectUnreportedHotspots() {
        SomaticVariant hotspot =
                SomaticVariantTestFactory.builder().canonicalHgvsCodingImpact("1").hotspot(Hotspot.HOTSPOT).reported(false).build();
        SomaticVariant nonHotspot =
                SomaticVariantTestFactory.builder().canonicalHgvsCodingImpact("2").hotspot(Hotspot.NON_HOTSPOT).reported(false).build();

        List<ReportableVariant> variants = SomaticVariantSelector.selectNonDrivers(Lists.newArrayList(hotspot, nonHotspot),
                Lists.newArrayList(),
                Lists.newArrayList(),
                Lists.newArrayList());

        assertEquals(1, variants.size());
        ReportableVariant variant = variants.get(0);
        assertEquals("1", variant.canonicalHgvsCodingImpact());
        assertEquals(Hotspot.HOTSPOT, variant.hotspot());
    }

    @Test
    public void canSelectVariantsWithEvidence() {
        SomaticVariant withEvidence = SomaticVariantTestFactory.builder().canonicalHgvsCodingImpact("1").reported(false).build();
        SomaticVariant withoutEvidence = SomaticVariantTestFactory.builder().canonicalHgvsCodingImpact("2").reported(false).build();

        ProtectEvidence evidence =
                ProtectTestFactory.builder().gene(withEvidence.gene()).event(ProtectEventGenerator.variantEvent(withEvidence)).build();

        List<ReportableVariant> variants = SomaticVariantSelector.selectNonDrivers(Lists.newArrayList(withEvidence, withoutEvidence),
                Lists.newArrayList(),
                Lists.newArrayList(evidence),
                Lists.newArrayList());

        assertEquals(1, variants.size());
        ReportableVariant variant = variants.get(0);
        assertEquals("1", variant.canonicalHgvsCodingImpact());
    }

    @Test
    public void canSelectVariantsWithReportedPhaseSet() {
        SomaticVariant withMatch = SomaticVariantTestFactory.builder().gene("gene 1").reported(false).addLocalPhaseSets(1).build();
        SomaticVariant withoutMatch = SomaticVariantTestFactory.builder().gene("gene 2").reported(false).addLocalPhaseSets(2).build();
        SomaticVariant withoutPhase = SomaticVariantTestFactory.builder().gene("gene 3").reported(false).build();

        ReportableVariant withPhase = ReportableVariantTestFactory.builder().localPhaseSet(1).build();
        ReportableVariant noPhase = ReportableVariantTestFactory.builder().localPhaseSet(null).build();

        List<ReportableVariant> variants =
                SomaticVariantSelector.selectNonDrivers(Lists.newArrayList(withMatch, withoutMatch, withoutPhase),
                        Lists.newArrayList(withPhase, noPhase),
                        Lists.newArrayList(),
                        Lists.newArrayList());

        assertEquals(1, variants.size());
        ReportableVariant variant = variants.get(0);
        assertEquals("gene 1", variant.gene());
    }

    @Test
    public void canSelectVariantsRelevantForCuppa() {
        String cuppaGene = SomaticVariantSelector.CUPPA_GENES.iterator().next();
        SomaticVariant cuppaRelevant =
                SomaticVariantTestFactory.builder().gene(cuppaGene).type(VariantType.INDEL).repeatCount(4).reported(false).build();

        SomaticVariant cuppaTooManyRepeats =
                SomaticVariantTestFactory.builder().gene(cuppaGene).type(VariantType.INDEL).repeatCount(10).reported(false).build();

        SomaticVariant wrongGene =
                SomaticVariantTestFactory.builder().gene("wrong gene").type(VariantType.INDEL).repeatCount(4).reported(false).build();

        List<ReportableVariant> variants =
                SomaticVariantSelector.selectNonDrivers(Lists.newArrayList(cuppaRelevant, cuppaTooManyRepeats, wrongGene),
                        Lists.newArrayList(),
                        Lists.newArrayList(),
                        Lists.newArrayList());

        assertEquals(1, variants.size());
        ReportableVariant variant = variants.get(0);
        assertEquals(cuppaGene, variant.gene());
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
                SomaticVariantSelector.selectNonDrivers(Lists.newArrayList(nonsense, splice, missense, wrongGene),
                        Lists.newArrayList(),
                        Lists.newArrayList(),
                        Lists.newArrayList(driverGene));

        assertEquals(2, variants.size());

        assertNotNull(findByCanonicalEffect(variants, "nonsense"));
        assertNotNull(findByCanonicalEffect(variants, "missense"));
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