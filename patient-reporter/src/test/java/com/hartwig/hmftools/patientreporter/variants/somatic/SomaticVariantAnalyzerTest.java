package com.hartwig.hmftools.patientreporter.variants.somatic;

import static org.junit.Assert.assertEquals;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientreporter.variants.driver.DriverGeneView;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.ImmutableSomaticVariantImpl;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantTestBuilderFactory;
import com.hartwig.hmftools.patientreporter.PatientReporterTestFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class SomaticVariantAnalyzerTest {

    private static final String PASS_FILTER = "PASS";

    private static final CodingEffect SPLICE = CodingEffect.SPLICE;
    private static final CodingEffect MISSENSE = CodingEffect.MISSENSE;
    private static final CodingEffect SYNONYMOUS = CodingEffect.SYNONYMOUS;

    private static final String RIGHT_GENE = "RIGHT";
    private static final String WRONG_GENE = "WRONG";

    @Test
    public void onlyReportsAndCountsRelevantVariants() {
        DriverGeneView driverGeneView = PatientReporterTestFactory.createTestDriverGeneView(RIGHT_GENE, "AnyGene");

        List<SomaticVariant> variants =
                Lists.newArrayList(builder().gene(RIGHT_GENE).canonicalCodingEffect(MISSENSE).worstCodingEffect(MISSENSE).build(),
                        builder().gene(RIGHT_GENE).canonicalCodingEffect(SYNONYMOUS).worstCodingEffect(SYNONYMOUS).build(),
                        builder().gene(RIGHT_GENE).canonicalCodingEffect(SPLICE).worstCodingEffect(SPLICE).build(),
                        builder().gene(RIGHT_GENE)
                                .canonicalCodingEffect(SYNONYMOUS)
                                .worstCodingEffect(SYNONYMOUS)
                                .hotspot(Hotspot.HOTSPOT)
                                .build(),
                        builder().gene(WRONG_GENE).canonicalCodingEffect(MISSENSE).worstCodingEffect(MISSENSE).build(),
                        builder().gene(WRONG_GENE).canonicalCodingEffect(SYNONYMOUS).worstCodingEffect(SYNONYMOUS).build());

        SomaticVariantAnalysis analysis = SomaticVariantAnalyzer.run(variants, driverGeneView, Collections.emptyList());

        // Report the missense variant on RIGHT_GENE plus the synonymous hotspot.
        assertEquals(2, analysis.variantsToReport().size());
    }

    @NotNull
    private static ImmutableSomaticVariantImpl.Builder builder() {
        return SomaticVariantTestBuilderFactory.create().filter(PASS_FILTER);
    }
}
