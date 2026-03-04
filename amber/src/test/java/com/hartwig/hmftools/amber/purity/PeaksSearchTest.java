package com.hartwig.hmftools.amber.purity;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.amber.PositionEvidence;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.junit.Test;

public class PeaksSearchTest extends PurityTestBase
{
    @Test
    public void test()
    {
        List<PositionEvidence> evidence = new ArrayList<>();
        int start = 2_000_000;
        // Simulate 20% contamination.
        for(int i = 0; i < 4000; i++)
        {
            int position = start + i * 1000;
            int mode = i % 16;
            switch(mode)
            {
                case 0:
                    // Tumor ref,ref, Contaminant ref,ref
                    evidence.add(evidenceWithDepthAndAltCount(HumanChromosome._3, position, 1000, 0));
                    break;
                case 1, 3:
                    // Tumor ref, ref, Contaminant ref, alt
                    evidence.add(evidenceWithDepthAndAltCount(HumanChromosome._3, position, 1000, 100));
                    break;
                case 2:
                    // Tumor ref, ref, Contaminant alt, alt
                    evidence.add(evidenceWithDepthAndAltCount(HumanChromosome._3, position, 1000, 200));
                    break;
                case 4, 8:
                    // Tumor ref, alt, Contaminant ref, ref
                    evidence.add(evidenceWithDepthAndAltCount(HumanChromosome._3, position, 1000, 400));
                    break;
                case 5, 7, 9, 11:
                    // Tumor ref, alt, Contaminant ref, alt
                    evidence.add(evidenceWithDepthAndAltCount(HumanChromosome._3, position, 1000, 500));
                    break;
                case 6, 10:
                    // Tumor ref, alt, Contaminant alt, alt
                    evidence.add(evidenceWithDepthAndAltCount(HumanChromosome._3, position, 1000, 600));
                    break;
                case 12:
                    // Tumor alt, alt, Contaminant ref, ref
                    evidence.add(evidenceWithDepthAndAltCount(HumanChromosome._3, position, 1000, 800));
                    break;
                case 13, 14:
                    // Tumor alt, alt, Contaminant ref, alt
                    evidence.add(evidenceWithDepthAndAltCount(HumanChromosome._3, position, 1000, 900));
                    break;
                case 15:
                    // Tumor alt, alt, Contaminant alt, alt
                    evidence.add(evidenceWithDepthAndAltCount(HumanChromosome._3, position, 1000, 1000));
            }
        }
        PeakSearch search = new PeakSearch(evidence, 4);
        List<CandidatePeakEvaluationResult> peaks = search.peaks();
        assertEquals(2, peaks.size());
        assertEquals(0.1, peaks.get(0).candidatePeak().vaf(), 0.01);
        assertEquals(0.2, peaks.get(1).candidatePeak().vaf(), 0.01);

        // Now add a copy number event with a peak at 30%.
        start = 1_000_000;
        for(int i = 0; i < 100; i++)
        {
            int position = start + i * 1000;
            evidence.add(evidenceWithDepthAndAltCount(HumanChromosome._3, position, 1000, 300));
        }
        search = new PeakSearch(evidence, 4);
        peaks = search.peaks();
        assertEquals(3, peaks.size());
        assertEquals(0.1, peaks.get(0).candidatePeak().vaf(), 0.01);
        assertEquals(0.2, peaks.get(1).candidatePeak().vaf(), 0.01);
        assertEquals(0.3, peaks.get(2).candidatePeak().vaf(), 0.01);
    }
}
