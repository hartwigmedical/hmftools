package com.hartwig.hmftools.peach.data_loader;

import static com.hartwig.hmftools.peach.TestUtils.getTestResourcePath;
import static com.hartwig.hmftools.peach.TestUtils.loadTestHaplotypePanel;

import static org.junit.Assert.assertEquals;

import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.peach.panel.HaplotypePanel;

import org.junit.Test;

public class HaplotypeEventLoaderTest
{
    @Test
    public void testLoad()
    {
        String vcfPath = getTestResourcePath("variants.vcf");
        String sampleName = "FAKER";
        Map<Chromosome, Set<Integer>> relevantVariantPositions = Map.of(
                HumanChromosome._1, Set.of(98205966, 98348885),
                HumanChromosome._2, Set.of(234668879, 234668880),
                HumanChromosome._4, Set.of(1733823)
        );

        Map<String, Integer> haplotypeEventsToCount = HaplotypeEventLoader.loadRelevantVariantHaplotypeEvents(
                vcfPath, sampleName, relevantVariantPositions
        );
        Map<String, Integer> expectedHaplotypeEventsToCount = Map.of(
                "VAR_chr1_98205966_GATGA_G", 1,
                "VAR_chr1_98205966_G_C", 1,
                "VAR_chr1_98348885_G_A", 2,
                "VAR_chr2_234668879_C_CAT", 2
        );
        assertEquals(expectedHaplotypeEventsToCount, haplotypeEventsToCount);
    }

    @Test
    public void testLoadZipped()
    {
        String vcfPath = getTestResourcePath("variants.vcf.gz");
        String sampleName = "FAKER";
        Map<Chromosome, Set<Integer>> relevantVariantPositions = Map.of(
                HumanChromosome._1, Set.of(98205966, 98348885),
                HumanChromosome._2, Set.of(234668879, 234668880),
                HumanChromosome._4, Set.of(1733823)
        );

        Map<String, Integer> haplotypeEventsToCount = HaplotypeEventLoader.loadRelevantVariantHaplotypeEvents(
                vcfPath, sampleName, relevantVariantPositions
        );
        Map<String, Integer> expectedHaplotypeEventsToCount = Map.of(
                "VAR_chr1_98205966_GATGA_G", 1,
                "VAR_chr1_98205966_G_C", 1,
                "VAR_chr1_98348885_G_A", 2,
                "VAR_chr2_234668879_C_CAT", 2
        );
        assertEquals(expectedHaplotypeEventsToCount, haplotypeEventsToCount);
    }

    @Test
    public void testLoadWithPanel()
    {
        String vcfPath = getTestResourcePath("variants.vcf");
        String sampleName = "FAKER";
        HaplotypePanel panel = loadTestHaplotypePanel("haplotypes.complicated.37.tsv");
        
        Map<String, Integer> haplotypeEventsToCount = HaplotypeEventLoader.loadRelevantVariantHaplotypeEvents(
                vcfPath, sampleName, panel.getRelevantVariantPositions()
        );
        Map<String, Integer> expectedHaplotypeEventsToCount = Map.of(
                "VAR_chr1_98205966_GATGA_G", 1,
                "VAR_chr1_98205966_G_C", 1,
                "VAR_chr1_98348885_G_A", 2,
                "VAR_chr2_234668879_C_CAT", 2
        );
        assertEquals(expectedHaplotypeEventsToCount, haplotypeEventsToCount);
    }
}
