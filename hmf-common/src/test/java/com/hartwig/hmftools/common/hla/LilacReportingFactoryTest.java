package com.hartwig.hmftools.common.hla;

import java.util.List;
import java.util.Map;

import static junit.framework.TestCase.assertEquals;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.lilac.LilacTestFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class LilacReportingFactoryTest {

    @NotNull
    private static LilacSummaryData createTestLilacData() {
        List<LilacAllele> alleles = Lists.newArrayList();
        alleles.add(LilacTestFactory.builder().allele("A*03:01").somaticMissense(2D).tumorCopyNumber(4.7).build());
        alleles.add(LilacTestFactory.builder().allele("A*03:01").somaticMissense(0D).tumorCopyNumber(1.5).build());
        alleles.add(LilacTestFactory.builder().allele("B*18:02").somaticSplice(1D).tumorCopyNumber(1.2).somaticNonsenseOrFrameshift(1D).build());
        alleles.add(LilacTestFactory.builder().allele("B*35:02").tumorCopyNumber(1.1).build());
        alleles.add(LilacTestFactory.builder().allele("C*10:12").somaticSynonymous(1D).tumorCopyNumber(0).build());
        alleles.add(LilacTestFactory.builder().allele("C*16:02").tumorCopyNumber(0).build());

        return ImmutableLilacSummaryData.builder().qc("PASS").alleles(alleles).build();
    }

    @Test
    public void testGene() {
        Map<String, List<LilacAllele>> mapLilac = LilacReportingFactory.generateLilacMap(createTestLilacData());
        assertEquals(LilacReportingFactory.extractHLAGene(mapLilac.get("A*03:01").get(0)), "HLA-A");
        assertEquals(LilacReportingFactory.extractHLAGene(mapLilac.get("B*18:02").get(0)), "HLA-B");
        assertEquals(LilacReportingFactory.extractHLAGene(mapLilac.get("B*35:02").get(0)), "HLA-B");
        assertEquals(LilacReportingFactory.extractHLAGene(mapLilac.get("C*10:12").get(0)), "HLA-C");
        assertEquals(LilacReportingFactory.extractHLAGene(mapLilac.get("C*16:02").get(0)), "HLA-C");
    }

    @Test
    public void testMutationString() {
        Map<String, List<LilacAllele>> mapLilac = LilacReportingFactory.generateLilacMap(createTestLilacData());
        assertEquals(LilacReportingFactory.mutationString(mapLilac.get("A*03:01").get(0)), "2 missense");
        assertEquals(LilacReportingFactory.mutationString(mapLilac.get("B*18:02").get(0)), "1 nonsense or frameshift, 1 splice");
        assertEquals(LilacReportingFactory.mutationString(mapLilac.get("B*35:02").get(0)), "No");
        assertEquals(LilacReportingFactory.mutationString(mapLilac.get("C*10:12").get(0)), "1 synonymous");
        assertEquals(LilacReportingFactory.mutationString(mapLilac.get("C*16:02").get(0)), "No");
    }

    @Test
    public void testReliableInterpretation() {
        Map<String, List<LilacAllele>> mapLilac = LilacReportingFactory.generateLilacMap(createTestLilacData());
        assertEquals(LilacReportingFactory.HLApresenceInTumor(mapLilac.get("A*03:01").get(0),
                LilacReportingFactory.mutationString(mapLilac.get("A*03:01").get(0)),
                true), "Yes, but mutation(s) detected");
        assertEquals(LilacReportingFactory.HLApresenceInTumor(mapLilac.get("B*18:02").get(0),
                LilacReportingFactory.mutationString(mapLilac.get("B*18:02").get(0)),
                true), "Yes, but mutation(s) detected");
        assertEquals(LilacReportingFactory.HLApresenceInTumor(mapLilac.get("B*35:02").get(0),
                LilacReportingFactory.mutationString(mapLilac.get("B*35:02").get(0)),
                true), "Yes");
        assertEquals(LilacReportingFactory.HLApresenceInTumor(mapLilac.get("C*10:12").get(0),
                LilacReportingFactory.mutationString(mapLilac.get("C*10:12").get(0)),
                true), "Unknown");
        assertEquals(LilacReportingFactory.HLApresenceInTumor(mapLilac.get("C*16:02").get(0),
                LilacReportingFactory.mutationString(mapLilac.get("C*16:02").get(0)),
                true), "No");
    }

    @Test
    public void testUnreliableInterpretation() {
        Map<String, List<LilacAllele>> mapLilac = LilacReportingFactory.generateLilacMap(createTestLilacData());
        assertEquals(LilacReportingFactory.HLApresenceInTumor(mapLilac.get("A*03:01").get(0),
                LilacReportingFactory.mutationString(mapLilac.get("A*03:01").get(0)),
                false), "Unknown");
        assertEquals(LilacReportingFactory.HLApresenceInTumor(mapLilac.get("B*18:02").get(0),
                LilacReportingFactory.mutationString(mapLilac.get("B*18:02").get(0)),
                false), "Unknown");
        assertEquals(LilacReportingFactory.HLApresenceInTumor(mapLilac.get("B*35:02").get(0),
                LilacReportingFactory.mutationString(mapLilac.get("B*35:02").get(0)),
                false), "Unknown");
        assertEquals(LilacReportingFactory.HLApresenceInTumor(mapLilac.get("C*10:12").get(0),
                LilacReportingFactory.mutationString(mapLilac.get("C*10:12").get(0)),
                false), "Unknown");
        assertEquals(LilacReportingFactory.HLApresenceInTumor(mapLilac.get("C*16:02").get(0),
                LilacReportingFactory.mutationString(mapLilac.get("C*16:02").get(0)),
                false), "Unknown");
    }

    @Test
    public void testConvertDataReliable() {
        LilacReportingData lilacReportingData = LilacReportingFactory.convertToReportData(createTestLilacData(), true);
        Map<String, List<LilacReporting>> lilacReporting = lilacReportingData.lilacReporting();

        LilacReporting lilacReporting1 = extractLilacReporting("A*03:01", lilacReporting.get("HLA-A"));
        assertEquals(lilacReporting1.lilacGermlineAllele().gene(), "HLA-A");
        assertEquals(lilacReporting1.lilacGermlineAllele().germlineAllele(), "A*03:01");
        assertEquals(lilacReporting1.somaticMutations(), "2 missense");
        assertEquals(lilacReporting1.interpretation(), "Yes, but mutation(s) detected");
        assertEquals(lilacReporting1.tumorCopies(), 6,2);
        assertEquals(lilacReporting1.germlineCopies(), 2D);

        LilacReporting lilacReporting2 = extractLilacReporting("B*18:02", lilacReporting.get("HLA-B"));
        assertEquals(lilacReporting2.lilacGermlineAllele().gene(), "HLA-B");
        assertEquals(lilacReporting2.lilacGermlineAllele().germlineAllele(), "B*18:02");
        assertEquals(lilacReporting2.somaticMutations(), "1 nonsense or frameshift, 1 splice");
        assertEquals(lilacReporting2.interpretation(), "Yes, but mutation(s) detected");
        assertEquals(lilacReporting2.tumorCopies(), 1.2);
        assertEquals(lilacReporting2.germlineCopies(), 1D);

        LilacReporting lilacReporting3 = extractLilacReporting("B*35:02", lilacReporting.get("HLA-B"));
        assertEquals(lilacReporting3.lilacGermlineAllele().gene(), "HLA-B");
        assertEquals(lilacReporting3.lilacGermlineAllele().germlineAllele(), "B*35:02");
        assertEquals(lilacReporting3.somaticMutations(), "No");
        assertEquals(lilacReporting3.interpretation(), "Yes");
        assertEquals(lilacReporting3.tumorCopies(), 1,1);
        assertEquals(lilacReporting3.germlineCopies(), 1D);

        LilacReporting lilacReporting4 = extractLilacReporting("C*10:12", lilacReporting.get("HLA-C"));
        assertEquals(lilacReporting4.lilacGermlineAllele().gene(), "HLA-C");
        assertEquals(lilacReporting4.lilacGermlineAllele().germlineAllele(), "C*10:12");
        assertEquals(lilacReporting4.somaticMutations(), "1 synonymous");
        assertEquals(lilacReporting4.interpretation(), "Unknown");
        assertEquals(lilacReporting4.tumorCopies(), 0D);
        assertEquals(lilacReporting4.germlineCopies(), 1D);

        LilacReporting lilacReporting5 = extractLilacReporting("C*16:02", lilacReporting.get("HLA-C"));
        assertEquals(lilacReporting5.lilacGermlineAllele().gene(), "HLA-C");
        assertEquals(lilacReporting5.lilacGermlineAllele().germlineAllele(), "C*16:02");
        assertEquals(lilacReporting5.somaticMutations(), "No");
        assertEquals(lilacReporting5.interpretation(), "No");
        assertEquals(lilacReporting5.tumorCopies(), 0D);
        assertEquals(lilacReporting5.germlineCopies(), 1D);
    }

    @NotNull
    public LilacReporting extractLilacReporting(@NotNull String germlineAllele, @NotNull List<LilacReporting> lilacReportingData) {
        for (LilacReporting lilacReporting: lilacReportingData) {
            if (lilacReporting.lilacGermlineAllele().germlineAllele().equals(germlineAllele)){
                return lilacReporting;
            }
        }
        throw new IllegalStateException("Could not find lilac reporting: " + germlineAllele);
    }
}