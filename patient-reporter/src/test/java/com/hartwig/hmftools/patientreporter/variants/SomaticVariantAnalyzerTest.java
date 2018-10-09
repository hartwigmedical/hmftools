package com.hartwig.hmftools.patientreporter.variants;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.ImmutableEnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantTestBuilderFactory;

import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Ignore;
import org.junit.Test;

public class SomaticVariantAnalyzerTest {

    private static final String PASS_FILTER = "PASS";

    private static final CodingEffect SPLICE = CodingEffect.SPLICE;
    private static final CodingEffect MISSENSE = CodingEffect.MISSENSE;
    private static final CodingEffect SYNONYMOUS = CodingEffect.SYNONYMOUS;
    private static final String RIGHT_GENE = "RIGHT";
    private static final String WRONG_GENE = "WRONG";

    @Test
    public void onlyReportsAndCountsRelevantVariants() throws IOException, ParseException {
        List<EnrichedSomaticVariant> variants =
                Lists.newArrayList(builder().gene(RIGHT_GENE).canonicalCodingEffect(MISSENSE).worstCodingEffect(MISSENSE).build(),
                        builder().gene(RIGHT_GENE).canonicalCodingEffect(SYNONYMOUS).worstCodingEffect(SYNONYMOUS).build(),
                        builder().gene(RIGHT_GENE).canonicalCodingEffect(SPLICE).worstCodingEffect(SPLICE).build(),
                        builder().gene(WRONG_GENE).canonicalCodingEffect(MISSENSE).worstCodingEffect(MISSENSE).build(),
                        builder().gene(WRONG_GENE).canonicalCodingEffect(SYNONYMOUS).worstCodingEffect(SYNONYMOUS).build());

//        SomaticVariantAnalysis analysis = SomaticVariantAnalyzer.run(variants,
//                Sets.newHashSet(RIGHT_GENE),
//                Maps.newHashMap(),
//                Strings.EMPTY,
//                Lists.newArrayList(),
//                Strings.EMPTY);
//
//        assertEquals(2, analysis.tumorMutationalLoad());
//        assertEquals(2, analysis.variantsToReport().size());
//
//        Map<String, DriverCategory> driverCategoryMap = Maps.newHashMap();
//        driverCategoryMap.put(RIGHT_GENE, DriverCategory.ONCO);
//        SomaticVariantAnalysis analysisOnco = SomaticVariantAnalyzer.run(variants,
//                Sets.newHashSet(RIGHT_GENE),
//                driverCategoryMap,
//                Strings.EMPTY,
//                Lists.newArrayList(),
//                Strings.EMPTY);
//
//        assertEquals(2, analysisOnco.tumorMutationalLoad());
//        assertEquals(1, analysisOnco.variantsToReport().size());
    }

    @NotNull
    private static ImmutableEnrichedSomaticVariant.Builder builder() {
        return SomaticVariantTestBuilderFactory.createEnriched().filter(PASS_FILTER);
    }
}
