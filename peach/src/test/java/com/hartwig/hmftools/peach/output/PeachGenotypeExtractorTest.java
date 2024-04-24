package com.hartwig.hmftools.peach.output;

import static com.hartwig.hmftools.peach.TestUtils.getTestResourcePath;

import static junit.framework.TestCase.assertEquals;

import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.peach.ImmutablePeachGenotype;
import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.peach.HaplotypeAnalysis;
import com.hartwig.hmftools.peach.PeachQCStatus;
import com.hartwig.hmftools.peach.data_loader.DrugInfoLoader;
import com.hartwig.hmftools.peach.data_loader.HaplotypeFunctionLoader;
import com.hartwig.hmftools.peach.effect.DrugInfo;
import com.hartwig.hmftools.peach.effect.DrugInfoStore;
import com.hartwig.hmftools.peach.effect.HaplotypeFunction;
import com.hartwig.hmftools.peach.effect.HaplotypeFunctionStore;
import com.hartwig.hmftools.peach.haplotype.HaplotypeCombination;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PeachGenotypeExtractorTest
{
    @Test
    public void testEmpty()
    {
        DrugInfoStore drugInfoStore = createTestDrugInfoStore();
        HaplotypeFunctionStore haplotypeFunctionStore = createHaplotypeFunctionStore();

        assertEquals(Collections.emptyList(), PeachGenotypeExtractor.extract(new HashMap<>(), null, null));
        assertEquals(Collections.emptyList(), PeachGenotypeExtractor.extract(new HashMap<>(), null, haplotypeFunctionStore));
        assertEquals(Collections.emptyList(), PeachGenotypeExtractor.extract(new HashMap<>(), drugInfoStore, null));
        assertEquals(Collections.emptyList(), PeachGenotypeExtractor.extract(new HashMap<>(), drugInfoStore, haplotypeFunctionStore));
    }

    @Test
    public void testNonEmptyNoDrugOrFunctionInfo()
    {
        Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis = createTestFakeGeneToHaplotypeAnalysis();
        List<PeachGenotype> genotypes = PeachGenotypeExtractor.extract(geneToHaplotypeAnalysis, null, null);
        List<PeachGenotype> expectedGenotypes = List.of(
                ImmutablePeachGenotype.builder()
                        .gene("FAKE1")
                        .allele("*1")
                        .alleleCount(2)
                        .function("")
                        .linkedDrugs("")
                        .urlPrescriptionInfo("")
                        .build(),
                ImmutablePeachGenotype.builder()
                        .gene("FAKE2")
                        .allele("Unresolved Haplotype")
                        .alleleCount(2)
                        .function("")
                        .linkedDrugs("")
                        .urlPrescriptionInfo("")
                        .build(),
                ImmutablePeachGenotype.builder()
                        .gene("FAKE3")
                        .allele("*3")
                        .alleleCount(1)
                        .function("")
                        .linkedDrugs("")
                        .urlPrescriptionInfo("")
                        .build(),
                ImmutablePeachGenotype.builder()
                        .gene("FAKE3")
                        .allele("*4")
                        .alleleCount(1)
                        .function("")
                        .linkedDrugs("")
                        .urlPrescriptionInfo("")
                        .build()
        );
        assertEquals(expectedGenotypes, genotypes);
    }

    @Test
    public void testNonEmptyNoMatchingDrugOrFunction()
    {
        DrugInfoStore drugInfoStore = createTestDrugInfoStore();
        HaplotypeFunctionStore haplotypeFunctionStore = createHaplotypeFunctionStore();

        Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis = createTestFakeGeneToHaplotypeAnalysis();
        List<PeachGenotype> genotypes = PeachGenotypeExtractor.extract(geneToHaplotypeAnalysis, drugInfoStore, haplotypeFunctionStore);
        List<PeachGenotype> expectedGenotypes = List.of(
                ImmutablePeachGenotype.builder()
                        .gene("FAKE1")
                        .allele("*1")
                        .alleleCount(2)
                        .function("")
                        .linkedDrugs("")
                        .urlPrescriptionInfo("")
                        .build(),
                ImmutablePeachGenotype.builder()
                        .gene("FAKE2")
                        .allele("Unresolved Haplotype")
                        .alleleCount(2)
                        .function("")
                        .linkedDrugs("")
                        .urlPrescriptionInfo("")
                        .build(),
                ImmutablePeachGenotype.builder()
                        .gene("FAKE3")
                        .allele("*3")
                        .alleleCount(1)
                        .function("")
                        .linkedDrugs("")
                        .urlPrescriptionInfo("")
                        .build(),
                ImmutablePeachGenotype.builder()
                        .gene("FAKE3")
                        .allele("*4")
                        .alleleCount(1)
                        .function("")
                        .linkedDrugs("")
                        .urlPrescriptionInfo("")
                        .build()
        );
        assertEquals(expectedGenotypes, genotypes);
    }

    @Test
    public void testNonEmptyWithMatchingDrugAndFunction()
    {
        DrugInfoStore drugInfoStore = createTestDrugInfoStore();
        HaplotypeFunctionStore haplotypeFunctionStore = createHaplotypeFunctionStore();

        Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis = new HashMap<>(createTestFakeGeneToHaplotypeAnalysis());
        geneToHaplotypeAnalysis.put(
                "DPYD",
                new HaplotypeAnalysis(
                        Map.of("EVENT_DPYD_1", 1),
                        List.of(new HaplotypeCombination(Map.of("*2A", 2))),
                        "*1",
                        "*1",
                        PeachQCStatus.PASS,
                        new HaplotypeCombination(Map.of("*2A", 2))
                )
        );
        geneToHaplotypeAnalysis.put(
                "UGT1A1",
                new HaplotypeAnalysis(
                        Map.of("EVENT_UGT1A1_1", 1),
                        List.of(new HaplotypeCombination(Map.of("*28", 1, "*1", 1))),
                        "*1",
                        "*1",
                        PeachQCStatus.PASS,
                        new HaplotypeCombination(Map.of("*28", 1, "*1", 1))
                )
        );

        List<PeachGenotype> genotypes = PeachGenotypeExtractor.extract(geneToHaplotypeAnalysis, drugInfoStore, haplotypeFunctionStore);
        List<PeachGenotype> expectedGenotypes = List.of(
                ImmutablePeachGenotype.builder()
                        .gene("DPYD")
                        .allele("*2A")
                        .alleleCount(2)
                        .function("No Function")
                        .linkedDrugs("5-Fluorouracil;Capecitabine;Tegafur")
                        .urlPrescriptionInfo("https://www.pharmgkb.org/guidelineAnnotation/PA166104939;https://www.pharmgkb.org/guidelineAnnotation/PA166104963;https://www.pharmgkb.org/guidelineAnnotation/PA166104944")
                        .build(),
                ImmutablePeachGenotype.builder()
                        .gene("FAKE1")
                        .allele("*1")
                        .alleleCount(2)
                        .function("")
                        .linkedDrugs("")
                        .urlPrescriptionInfo("")
                        .build(),
                ImmutablePeachGenotype.builder()
                        .gene("FAKE2")
                        .allele("Unresolved Haplotype")
                        .alleleCount(2)
                        .function("")
                        .linkedDrugs("")
                        .urlPrescriptionInfo("")
                        .build(),
                ImmutablePeachGenotype.builder()
                        .gene("FAKE3")
                        .allele("*3")
                        .alleleCount(1)
                        .function("")
                        .linkedDrugs("")
                        .urlPrescriptionInfo("")
                        .build(),
                ImmutablePeachGenotype.builder()
                        .gene("FAKE3")
                        .allele("*4")
                        .alleleCount(1)
                        .function("")
                        .linkedDrugs("")
                        .urlPrescriptionInfo("")
                        .build(),
                ImmutablePeachGenotype.builder()
                        .gene("UGT1A1")
                        .allele("*1")
                        .alleleCount(1)
                        .function("Normal Function")
                        .linkedDrugs("Irinotecan")
                        .urlPrescriptionInfo("https://www.pharmgkb.org/guidelineAnnotation/PA166104951")
                        .build(),
                ImmutablePeachGenotype.builder()
                        .gene("UGT1A1")
                        .allele("*28")
                        .alleleCount(1)
                        .function("Reduced Function")
                        .linkedDrugs("Irinotecan")
                        .urlPrescriptionInfo("https://www.pharmgkb.org/guidelineAnnotation/PA166104951")
                        .build()
        );
        assertEquals(expectedGenotypes, genotypes);
    }

    @NotNull
    private static Map<String, HaplotypeAnalysis> createTestFakeGeneToHaplotypeAnalysis()
    {
        HaplotypeAnalysis fake1HaplotypeAnalysis = new HaplotypeAnalysis(
                new HashMap<>(),
                List.of(new HaplotypeCombination(Map.of("*1", 2))),
                "*1",
                "*1",
                PeachQCStatus.PASS,
                new HaplotypeCombination(Map.of("*1", 2))
        );
        HaplotypeAnalysis fake2HaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of("EVENT_1", 1, "EVENT_3", 2),
                List.of(
                        new HaplotypeCombination(Map.of("2373C>T", 1, "*9", 1)),
                        new HaplotypeCombination(Map.of("*2", 1, "*5", 1))
                ),
                "*9",
                "*1",
                PeachQCStatus.FAIL_NO_UNIQUE_BEST_COMBINATION_FOUND,
                null
        );
        HaplotypeAnalysis fake3HaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of("EVENT_1", 2),
                List.of(new HaplotypeCombination(Map.of("*4", 1, "*3", 1))),
                "*9",
                "*1",
                PeachQCStatus.PASS,
                new HaplotypeCombination(Map.of("*4", 1, "*3", 1))
        );
        return Map.of("FAKE3", fake3HaplotypeAnalysis, "FAKE2", fake2HaplotypeAnalysis, "FAKE1", fake1HaplotypeAnalysis);
    }

    @NotNull
    private static HaplotypeFunctionStore createHaplotypeFunctionStore()
    {
        String functionResourcePath = getTestResourcePath("function.tsv");
        List<HaplotypeFunction> functions = HaplotypeFunctionLoader.loadFunctions(functionResourcePath);
        return new HaplotypeFunctionStore(functions);
    }

    @NotNull
    private static DrugInfoStore createTestDrugInfoStore()
    {
        String drugResourcePath = getTestResourcePath("drugs.tsv");
        List<DrugInfo> drugInfos = DrugInfoLoader.loadDrugInfos(drugResourcePath);
        return new DrugInfoStore(drugInfos);
    }
}
