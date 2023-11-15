package com.hartwig.hmftools.peach;

import com.google.common.io.Resources;
import com.hartwig.hmftools.peach.data_loader.PanelLoader;
import com.hartwig.hmftools.peach.haplotype.HaplotypeCombination;
import com.hartwig.hmftools.peach.panel.HaplotypePanel;
import org.junit.Test;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static junit.framework.TestCase.assertEquals;
import static org.junit.Assert.assertTrue;

public class HaplotypeCallerTest
{
    @Test
    public void testEmptyRun()
    {
        HaplotypePanel haplotypePanel = new HaplotypePanel(new HashMap<>());
        HaplotypeCaller caller = new HaplotypeCaller(haplotypePanel);
        Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis = caller.getGeneToHaplotypeAnalysis(new HashMap<>());
        assertEquals(0, geneToHaplotypeAnalysis.size());
    }

    @Test
    public void testNoEvents()
    {
        HaplotypePanel haplotypePanel = loadTestHaplotypePanel("haplotypes.complicated.37.tsv");
        HaplotypeCaller caller = new HaplotypeCaller(haplotypePanel);
        Map<String, HaplotypeAnalysis> geneToHaplotypeAnalysis = caller.getGeneToHaplotypeAnalysis(new HashMap<>());
        
        HaplotypeAnalysis expectedDpydHaplotypeAnalysis = new HaplotypeAnalysis(
                new HashMap<>(),
                List.of(new HaplotypeCombination(Map.of("*9A", 2))),
                "*9A",
                "*1"
        );
        HaplotypeAnalysis expectedUgt1A1HaplotypeAnalysis = new HaplotypeAnalysis(
                new HashMap<>(),
                List.of(new HaplotypeCombination(Map.of("*1", 2))),
                "*1",
                "*1"
        );
        assertEquals(2, geneToHaplotypeAnalysis.size());
        assertTrue(geneToHaplotypeAnalysis.containsKey("DPYD"));
        assertTrue(geneToHaplotypeAnalysis.containsKey("UGT1A1"));
        
        assertEqualHaplotypeAnalysis(expectedDpydHaplotypeAnalysis, geneToHaplotypeAnalysis.get("DPYD"));
        assertEqualHaplotypeAnalysis(expectedUgt1A1HaplotypeAnalysis, geneToHaplotypeAnalysis.get("UGT1A1"));
    }
    
    private void assertEqualHaplotypeAnalysis(HaplotypeAnalysis expected, HaplotypeAnalysis actual)
    {
        assertEquals(expected.getDefaultHaplotypeName(), actual.getDefaultHaplotypeName());
        assertEquals(expected.getWildTypeHaplotypeName(), actual.getWildTypeHaplotypeName());
        
        assertEquals(expected.getEventIds(), actual.getEventIds());
        for (String eventId : expected.getEventIds())
        {
            assertEquals(
                    String.format("Compare event counts of %s", eventId), 
                    expected.getEventCount(eventId), 
                    actual.getEventCount(eventId)
            );
        }
        assertEquals(expected.getHaplotypeCombinations(), actual.getHaplotypeCombinations());
    }
        

    private HaplotypePanel loadTestHaplotypePanel(String fileName)
    {
        String filePath = Resources.getResource(fileName).getPath();
        return PanelLoader.loadHaplotypePanel(filePath);
    }
}
