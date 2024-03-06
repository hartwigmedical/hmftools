package com.hartwig.hmftools.peach.output;

import static com.hartwig.hmftools.peach.TestUtils.getTestResourcePath;

import static junit.framework.TestCase.assertEquals;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.peach.HaplotypeAnalysis;
import com.hartwig.hmftools.peach.data_loader.DrugInfoLoader;
import com.hartwig.hmftools.peach.data_loader.HaplotypeFunctionLoader;
import com.hartwig.hmftools.peach.effect.DrugInfo;
import com.hartwig.hmftools.peach.effect.DrugInfoStore;
import com.hartwig.hmftools.peach.effect.HaplotypeFunction;
import com.hartwig.hmftools.peach.effect.HaplotypeFunctionStore;
import com.hartwig.hmftools.peach.haplotype.HaplotypeCombination;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class BestHaplotypeCombinationsFileTest
{
    private static final String EXPECTED_HEADER = "gene\tallele\tcount\tfunction\tlinkedDrugs\tprescriptionUrls";

    @Test
    public void testEmpty()
    {
        DrugInfoStore drugInfoStore = createTestDrugInfoStore();
        HaplotypeFunctionStore haplotypeFunctionStore = createHaplotypeFunctionStore();

        assertEquals(List.of(EXPECTED_HEADER), BestHaplotypeCombinationsFile.toLines(new HashMap<>(), null, null));
        assertEquals(List.of(EXPECTED_HEADER), BestHaplotypeCombinationsFile.toLines(new HashMap<>(), null, haplotypeFunctionStore));
        assertEquals(List.of(EXPECTED_HEADER), BestHaplotypeCombinationsFile.toLines(new HashMap<>(), drugInfoStore, null));
        assertEquals(List.of(EXPECTED_HEADER), BestHaplotypeCombinationsFile.toLines(new HashMap<>(), drugInfoStore, haplotypeFunctionStore));
    }

    @Test
    public void testNonEmptyNoDrugOrFunctionInfo()
    {
        Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis = createTestFakeGeneToHaplotypeAnalysis();
        List<String> outputLines = BestHaplotypeCombinationsFile.toLines(geneToHaplotypeAnalysis, null, null);
        List<String> expectedLines =
                List.of(EXPECTED_HEADER, "FAKE1\t*1\t2\t\t\t", "FAKE2\tUnresolved Haplotype\t2\t\t\t", "FAKE3\t*3\t1\t\t\t", "FAKE3\t*4\t1\t\t\t");
        assertEquals(expectedLines, outputLines);
    }

    @Test
    public void testNonEmptyNoMatchingDrugOrFunction()
    {
        DrugInfoStore drugInfoStore = createTestDrugInfoStore();
        HaplotypeFunctionStore haplotypeFunctionStore = createHaplotypeFunctionStore();

        Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis = createTestFakeGeneToHaplotypeAnalysis();
        List<String> outputLines = BestHaplotypeCombinationsFile.toLines(geneToHaplotypeAnalysis, drugInfoStore, haplotypeFunctionStore);
        List<String> expectedLines =
                List.of(EXPECTED_HEADER, "FAKE1\t*1\t2\t\t\t", "FAKE2\tUnresolved Haplotype\t2\t\t\t", "FAKE3\t*3\t1\t\t\t", "FAKE3\t*4\t1\t\t\t");
        assertEquals(expectedLines, outputLines);
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
                        "*1"
                )
        );
        geneToHaplotypeAnalysis.put(
                "UGT1A1",
                new HaplotypeAnalysis(
                        Map.of("EVENT_UGT1A1_1", 1),
                        List.of(new HaplotypeCombination(Map.of("*28", 1, "*1", 1))),
                        "*1",
                        "*1"
                )
        );

        List<String> outputLines = BestHaplotypeCombinationsFile.toLines(geneToHaplotypeAnalysis, drugInfoStore, haplotypeFunctionStore);
        List<String> expectedLines = List.of(
                EXPECTED_HEADER,
                "DPYD\t*2A\t2\tNo Function\t5-Fluorouracil;Capecitabine;Tegafur\thttps://www.pharmgkb.org/guidelineAnnotation/PA166104939;https://www.pharmgkb.org/guidelineAnnotation/PA166104963;https://www.pharmgkb.org/guidelineAnnotation/PA166104944",
                "FAKE1\t*1\t2\t\t\t",
                "FAKE2\tUnresolved Haplotype\t2\t\t\t",
                "FAKE3\t*3\t1\t\t\t",
                "FAKE3\t*4\t1\t\t\t",
                "UGT1A1\t*1\t1\tNormal Function\tIrinotecan\thttps://www.pharmgkb.org/guidelineAnnotation/PA166104951",
                "UGT1A1\t*28\t1\tReduced Function\tIrinotecan\thttps://www.pharmgkb.org/guidelineAnnotation/PA166104951"
        );
        assertEquals(expectedLines, outputLines);
    }

    @NotNull
    private static Map<String, HaplotypeAnalysis> createTestFakeGeneToHaplotypeAnalysis()
    {
        HaplotypeAnalysis fake1HaplotypeAnalysis = new HaplotypeAnalysis(
                new HashMap<>(),
                List.of(new HaplotypeCombination(Map.of("*1", 2))),
                "*1",
                "*1"
        );
        HaplotypeAnalysis fake2HaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of("EVENT_1", 1, "EVENT_3", 2),
                List.of(
                        new HaplotypeCombination(Map.of("2373C>T", 1, "*9", 1)),
                        new HaplotypeCombination(Map.of("*2", 1, "*5", 1))
                ),
                "*9",
                "*1"
        );
        HaplotypeAnalysis fake3HaplotypeAnalysis = new HaplotypeAnalysis(
                Map.of("EVENT_1", 2),
                List.of(new HaplotypeCombination(Map.of("*4", 1, "*3", 1))),
                "*9",
                "*1"
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
