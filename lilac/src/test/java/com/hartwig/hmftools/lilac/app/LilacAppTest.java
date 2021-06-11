package com.hartwig.hmftools.lilac.app;

import static com.hartwig.hmftools.lilac.LilacConstants.EXPECTED_ALLELE_COUNT;
import static com.hartwig.hmftools.lilac.LilacConstants.STOP_LOSS_ON_C_ALLELE;
import static com.hartwig.hmftools.lilac.LilacConstants.longGeneName;
import static com.hartwig.hmftools.lilac.ReferenceData.STOP_LOSS_ON_C_INDEL;
import static com.hartwig.hmftools.lilac.misc.LilacTestUtils.createFragment;
import static com.hartwig.hmftools.lilac.misc.LilacTestUtils.disableLogging;
import static com.hartwig.hmftools.lilac.qc.LilacQCStatus.PASS;
import static com.hartwig.hmftools.lilac.qc.LilacQCStatus.WARN_UNMATCHED_HAPLOTYPE;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.lilac.LilacApplication;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.ReferenceData;
import com.hartwig.hmftools.lilac.coverage.HlaComplexCoverage;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.qc.LilacQC;
import com.hartwig.hmftools.lilac.qc.SolutionSummary;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

import org.junit.Test;

public class LilacAppTest
{
    private static final String SAMPLE_TEST = "SAMPLE_ID";

    private static List<String> COMMON_ALLELES = Lists.newArrayList(
            "A*01:01", "A*02:01", "A*03:01", "B*07:02", "B*08:01", "C*01:02", "C*02:02");

    private int mReadIdCounter = 0;

    @Test
    public void basicApplicationTest()
    {
        disableLogging();

        LilacApplication lilac = new LilacApplication(new LilacConfig(SAMPLE_TEST));

        MockBamReader refBamReader = new MockBamReader();
        MockBamReader tumorBamReader = new MockBamReader();

        lilac.setBamReaders(refBamReader, tumorBamReader);

        ReferenceData refData = lilac.getReferenceData();
        loadTestReferenceData(refData);

        // build support for the 3 common alleles
        int fragCount = 50;
        int startLoci = 20;
        int endLoci = 1050;
        int length = 150;
        int gap = 20;
        HlaAllele a1 = refData.findAllele("A*01:01:01", false);
        refBamReader.Fragments.addAll(createFragments(refData, a1, fragCount, startLoci, endLoci, length, gap));

        HlaAllele a2 = refData.findAllele("A*02:01:01", false);
        refBamReader.Fragments.addAll(createFragments(refData, a2, fragCount, startLoci, endLoci, length, gap));

        HlaAllele b1 = refData.findAllele("B*07:02:01", false);
        refBamReader.Fragments.addAll(createFragments(refData, b1, fragCount, startLoci, endLoci, length, gap));

        HlaAllele b2 = refData.findAllele("B*08:01:01", false);
        refBamReader.Fragments.addAll(createFragments(refData, b2, fragCount, startLoci, endLoci, length, gap));

        HlaAllele c1 = refData.findAllele("C*01:02:01", false);
        refBamReader.Fragments.addAll(createFragments(refData, c1, fragCount, startLoci, endLoci, length, gap));

        HlaAllele c2 = refData.findAllele("C*02:02:01", false);
        refBamReader.Fragments.addAll(createFragments(refData, c2, fragCount, startLoci, endLoci, length, gap));

        lilac.run();

        // check various outputs

        assertEquals(1, lilac.getRankedComplexes().size());
        HlaComplexCoverage winningComplex = lilac.getRankedComplexes().get(0);
        assertEquals(0, winningComplex.homozygousCount());
        assertEquals(1, winningComplex.recoveredCount());
        assertEquals(-6.0, winningComplex.cohortFrequencyTotal(), 0.01);

        SolutionSummary solutionSummary = lilac.getSolutionSummary();

        assertEquals(EXPECTED_ALLELE_COUNT, solutionSummary.ReferenceCoverage.getAlleles().size());
        assertTrue(solutionSummary.ReferenceCoverage.getAlleles().contains(a1.asFourDigit()));
        assertTrue(solutionSummary.ReferenceCoverage.getAlleles().contains(a2.asFourDigit()));
        assertTrue(solutionSummary.ReferenceCoverage.getAlleles().contains(b1.asFourDigit()));
        assertTrue(solutionSummary.ReferenceCoverage.getAlleles().contains(b2.asFourDigit()));
        assertTrue(solutionSummary.ReferenceCoverage.getAlleles().contains(c1.asFourDigit()));
        assertTrue(solutionSummary.ReferenceCoverage.getAlleles().contains(c2.asFourDigit()));

        LilacQC qcMetrics = lilac.getSummaryMetrics();

        assertTrue(qcMetrics.Status.contains(PASS));
        assertFalse(qcMetrics.HasHlaY);
        assertEquals(0, qcMetrics.AminoAcidQC.UnusedAminoAcids);
        assertTrue(qcMetrics.HaplotypeQC.UnmatchedHaplotypes.isEmpty());
        assertEquals(fragCount * EXPECTED_ALLELE_COUNT, qcMetrics.CoverageQC.TotalFragments);
        assertEquals(fragCount * EXPECTED_ALLELE_COUNT, qcMetrics.CoverageQC.FittedFragments);
    }

    @Test
    public void knownStopLossTest()
    {
        LilacApplication lilac = new LilacApplication(new LilacConfig(SAMPLE_TEST));

        MockBamReader refBamReader = new MockBamReader();
        MockBamReader tumorBamReader = new MockBamReader();

        lilac.setBamReaders(refBamReader, tumorBamReader);

        ReferenceData refData = lilac.getReferenceData();
        loadTestReferenceData(refData);

        // build support for the 3 common alleles
        int fragCount = 50;
        int startLoci = 20;
        int endLoci = 1050;
        int length = 150;
        int gap = 20;
        HlaAllele a1 = refData.findAllele("A*01:01:01", false);
        refBamReader.Fragments.addAll(createFragments(refData, a1, fragCount, startLoci, endLoci, length, gap));

        HlaAllele a2 = refData.findAllele("A*02:01:01", false);
        refBamReader.Fragments.addAll(createFragments(refData, a2, fragCount, startLoci, endLoci, length, gap));

        HlaAllele b1 = refData.findAllele("B*07:02:01", false);
        refBamReader.Fragments.addAll(createFragments(refData, b1, fragCount, startLoci, endLoci, length, gap));

        HlaAllele b2 = refData.findAllele("B*08:01:01", false);
        refBamReader.Fragments.addAll(createFragments(refData, b2, fragCount, startLoci, endLoci, length, gap));

        HlaAllele c1 = refData.findAllele("C*01:02:01", false);
        refBamReader.Fragments.addAll(createFragments(refData, c1, fragCount, startLoci, endLoci, length, gap));

        HlaAllele c2 = refData.findAllele("C*02:02:01", false);
        refBamReader.Fragments.addAll(createFragments(refData, c2, fragCount, startLoci, endLoci, length, gap));

        HlaAllele c3 = refData.findAllele(STOP_LOSS_ON_C_ALLELE, false);
        List<Fragment> stopLossFrags = createFragments(refData, c3, 20, 900, 1090, length, gap);
        refBamReader.Fragments.addAll(stopLossFrags);
        refBamReader.Fragments.addAll(createFragments(refData, c3, fragCount, startLoci, endLoci, length, gap));
        refBamReader.StopLossFragments.put(STOP_LOSS_ON_C_INDEL, stopLossFrags);

        lilac.run();

        // check various outputs

        assertEquals(1, lilac.getRankedComplexes().size());
        HlaComplexCoverage winningComplex = lilac.getRankedComplexes().get(0);
        assertEquals(0, winningComplex.homozygousCount());
        assertEquals(1, winningComplex.recoveredCount());
        assertEquals(-9.0, winningComplex.cohortFrequencyTotal(), 0.01);

        SolutionSummary solutionSummary = lilac.getSolutionSummary();

        assertEquals(EXPECTED_ALLELE_COUNT, solutionSummary.ReferenceCoverage.getAlleles().size());
        assertTrue(solutionSummary.ReferenceCoverage.getAlleles().contains(a1.asFourDigit()));
        assertTrue(solutionSummary.ReferenceCoverage.getAlleles().contains(a2.asFourDigit()));
        assertTrue(solutionSummary.ReferenceCoverage.getAlleles().contains(b1.asFourDigit()));
        assertTrue(solutionSummary.ReferenceCoverage.getAlleles().contains(b2.asFourDigit()));
        assertTrue(solutionSummary.ReferenceCoverage.getAlleles().contains(c2.asFourDigit()));
        assertTrue(solutionSummary.ReferenceCoverage.getAlleles().contains(c3.asFourDigit()));

        LilacQC qcMetrics = lilac.getSummaryMetrics();

        // C*01:01 will cause unmatched haplotypes
        assertTrue(qcMetrics.Status.contains(WARN_UNMATCHED_HAPLOTYPE));
        assertFalse(qcMetrics.HasHlaY);
        assertEquals(0, qcMetrics.AminoAcidQC.UnusedAminoAcids);
        assertEquals(5, qcMetrics.HaplotypeQC.UnmatchedHaplotypes.size());
        assertEquals(370, qcMetrics.CoverageQC.TotalFragments);
        assertEquals(338, qcMetrics.CoverageQC.FittedFragments);
    }

    @Test
    public void wildcardTest()
    {
        disableLogging();

        LilacApplication lilac = new LilacApplication(new LilacConfig(SAMPLE_TEST));

        MockBamReader refBamReader = new MockBamReader();
        MockBamReader tumorBamReader = new MockBamReader();

        lilac.setBamReaders(refBamReader, tumorBamReader);

        ReferenceData refData = lilac.getReferenceData();
        loadTestReferenceData(refData);

        // build support for the 3 common alleles
        int fragCount = 150;
        int startLoci = 20;
        int endLoci = 1050;
        int length = 150;
        int gap = 20;
        HlaAllele a1 = refData.findAllele("A*01:01:01", false);
        refBamReader.Fragments.addAll(createFragments(refData, a1, 10, startLoci, endLoci, length, gap));

        HlaAllele a2 = refData.findAllele("A*02:01:01", false);
        refBamReader.Fragments.addAll(createFragments(refData, a2, fragCount, startLoci, endLoci, length, gap));

        HlaAllele a3 = refData.findAllele("A*03:01:01", false);
        refBamReader.Fragments.addAll(createFragments(refData, a3, fragCount * 2, startLoci, endLoci, length, gap));

        HlaAllele b1 = refData.findAllele("B*07:02:01", false);
        refBamReader.Fragments.addAll(createFragments(refData, b1, 10, startLoci, endLoci, length, gap));

        HlaAllele b2 = refData.findAllele("B*08:01:01", false);
        refBamReader.Fragments.addAll(createFragments(refData, b2, fragCount, startLoci, endLoci, length, gap));

        HlaAllele bwc = refData.findAllele("B*38:58", false);
        refBamReader.Fragments.addAll(createFragments(refData, bwc, fragCount * 2, startLoci, endLoci, length, gap));

        HlaAllele c1 = refData.findAllele("C*01:02:01", false);
        refBamReader.Fragments.addAll(createFragments(refData, c1, fragCount, startLoci, endLoci, length, gap));

        // HlaAllele c2 = refData.findAllele("C*02:02:01", false);
        // refBamReader.Fragments.addAll(createFragments(refData, c2, fragCount, startLoci, endLoci, length, gap));

        lilac.run();

        // check various outputs

        assertEquals(1, lilac.getRankedComplexes().size());
        HlaComplexCoverage winningComplex = lilac.getRankedComplexes().get(0);
        assertEquals(1, winningComplex.homozygousCount());
        assertEquals(0, winningComplex.recoveredCount());
        assertEquals(-9.0, winningComplex.cohortFrequencyTotal(), 0.01);

        SolutionSummary solutionSummary = lilac.getSolutionSummary();

        assertEquals(EXPECTED_ALLELE_COUNT, solutionSummary.ReferenceCoverage.getAlleles().size());
        assertTrue(solutionSummary.ReferenceCoverage.getAlleles().contains(a2.asFourDigit()));
        assertTrue(solutionSummary.ReferenceCoverage.getAlleles().contains(a3.asFourDigit()));
        assertTrue(solutionSummary.ReferenceCoverage.getAlleles().contains(b2.asFourDigit()));
        assertTrue(solutionSummary.ReferenceCoverage.getAlleles().contains(bwc.asFourDigit()));
        assertTrue(solutionSummary.ReferenceCoverage.getAlleles().contains(c1.asFourDigit()));

        LilacQC qcMetrics = lilac.getSummaryMetrics();

        // C*01:01 will cause unmatched haplotypes
        assertTrue(qcMetrics.Status.contains(PASS));
        assertFalse(qcMetrics.HasHlaY);
        assertEquals(0, qcMetrics.AminoAcidQC.UnusedAminoAcids);
        assertTrue(!qcMetrics.HaplotypeQC.UnmatchedHaplotypes.isEmpty());
        assertTrue(qcMetrics.CoverageQC.PercentWildcard > 0);
    }

    private List<Fragment> createFragments(
            ReferenceData refData, final HlaAllele allele, int fragCount, int startLoci, int endLoci, int length, int gap)
    {
        List<Fragment> fragments = Lists.newArrayList();

        HlaSequenceLoci sequenceLoci = refData.NucleotideSequences.stream().filter(x -> x.Allele.matches(allele)).findFirst().orElse(null);

        if(sequenceLoci == null)
            return fragments;

        int startLocus = startLoci;
        for(int i = 0; i < fragCount; ++i)
        {
            if(startLocus + length > endLoci)
            {
                startLocus = startLoci;
            }

            int endLocus = startLocus + length;

            fragments.add(createFragment(
                    String.format("READ_%03d", ++mReadIdCounter),  longGeneName(sequenceLoci.Allele.Gene),
                    sequenceLoci.sequence(), startLocus, endLocus));

            startLocus += gap;
        }

        return fragments;
    }

    private void loadTestReferenceData(final ReferenceData refData)
    {
        final List<String> nucleotides = new BufferedReader(new InputStreamReader(
                ReferenceData.class.getResourceAsStream("/test_allele_nucleotides.csv")))
                .lines().collect(Collectors.toList());

        refData.loadSequenceFile(nucleotides, refData.NucleotideSequences, false);

        final List<String> aminoAcids = new BufferedReader(new InputStreamReader(
                ReferenceData.class.getResourceAsStream("/test_allele_amino_acids.csv")))
                .lines().collect(Collectors.toList());

        refData.loadSequenceFile(aminoAcids, refData.AminoAcidSequences, true);

        for(String alleleStr : COMMON_ALLELES)
        {
            HlaSequenceLoci seqLoci = refData.AminoAcidSequences.stream().filter(x -> x.Allele.matches(alleleStr)).findFirst().orElse(null);

            if(seqLoci == null)
                continue;

            refData.getAlleleFrequencies().getAlleleFrequencies().put(seqLoci.Allele, 0.1);
        }
    }
}
