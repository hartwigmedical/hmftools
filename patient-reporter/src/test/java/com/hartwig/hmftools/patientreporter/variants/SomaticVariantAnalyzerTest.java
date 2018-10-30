package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.ImmutableEnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantTestBuilderFactory;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneFusion;

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
        List<EnrichedSomaticVariant> variants =
                Lists.newArrayList(builder().gene(RIGHT_GENE).canonicalCodingEffect(MISSENSE).worstCodingEffect(MISSENSE).build(),
                        builder().gene(RIGHT_GENE).canonicalCodingEffect(SYNONYMOUS).worstCodingEffect(SYNONYMOUS).build(),
                        builder().gene(RIGHT_GENE).canonicalCodingEffect(SPLICE).worstCodingEffect(SPLICE).build(),
                        builder().gene(WRONG_GENE).canonicalCodingEffect(MISSENSE).worstCodingEffect(MISSENSE).build(),
                        builder().gene(WRONG_GENE).canonicalCodingEffect(SYNONYMOUS).worstCodingEffect(SYNONYMOUS).build());

        //TODO: create genecopyNumber list
        List<GeneCopyNumber> geneCopyNumbers = Lists.newArrayList();
        List<GeneFusion> geneFusion = Lists.newArrayList();

//        SomaticVariantAnalysis analysis = SomaticVariantAnalyzer.run(variants,
//                Sets.newHashSet(RIGHT_GENE),
//                Maps.newHashMap(),
//                null, exomeGeneCopyNumbers, geneFusion);
//
//        assertEquals(2, analysis.tumorMutationalLoad());
//        assertEquals(2, analysis.reportableSomaticVariants().size());
//
//        Map<String, DriverCategory> driverCategoryMap = Maps.newHashMap();
//        driverCategoryMap.put(RIGHT_GENE, DriverCategory.ONCO);
//        SomaticVariantAnalysis analysisOnco = SomaticVariantAnalyzer.run(variants,
//                Sets.newHashSet(RIGHT_GENE),
//                driverCategoryMap,
//                null, exomeGeneCopyNumbers, geneFusion);
//
//        assertEquals(2, analysisOnco.tumorMutationalLoad());
//        assertEquals(1, analysisOnco.reportableSomaticVariants().size());
    }

    @NotNull
    private static ImmutableEnrichedSomaticVariant.Builder builder() {
        return SomaticVariantTestBuilderFactory.createEnriched().filter(PASS_FILTER);
    }
}
