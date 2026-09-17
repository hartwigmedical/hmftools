package com.hartwig.hmftools.purple.tools;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.purple.GermlineStatus.AMPLIFICATION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HET_DELETION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HOM_DELETION;

import static org.junit.Assert.assertEquals;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.StringWriter;
import java.util.List;
import java.util.TreeMap;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.junit.Test;

public class EventCountsTest extends ToolsTestBase
{
    @Test
    public void mergeWriteTest() throws IOException
    {
        RegionGeneEvents rge_1_1 = rge(_1, 1000, 2000, AMPLIFICATION);
        RegionGeneEvents rge_1_2 = rge(_1, 1000, 2000, HET_DELETION);
        RegionGeneEvents rge_1_3 = rge(_1, 2000, 3000, AMPLIFICATION);
        RegionGeneEvents rge_1_4 = rge(_1, 2000, 3000, HET_DELETION);
        RegionGeneEvents rge_1_5 = rge(_1, 2000, 3000, HOM_DELETION);
        RegionGeneEvents rge_1_6 = rge(_1, 3000, 4000, AMPLIFICATION);
        RegionGeneEvents rge_1_7 = rge(_1, 3000, 5000, AMPLIFICATION);
        RegionGeneEvents rge_1_8 = rge(_1, 4000, 5000, HET_DELETION);

        RegionGeneEvents rge_2_1 = rge(_2, 1000, 2000, AMPLIFICATION);
        RegionGeneEvents rge_2_2 = rge(_2, 1000, 2000, HET_DELETION);
        RegionGeneEvents rge_2_3 = rge(_2, 2000, 3000, AMPLIFICATION);
        RegionGeneEvents rge_2_4 = rge(_2, 2000, 3000, HET_DELETION);
        RegionGeneEvents rge_2_5 = rge(_2, 2000, 3000, HOM_DELETION);
        RegionGeneEvents rge_2_6 = rge(_2, 3000, 4000, AMPLIFICATION);
        RegionGeneEvents rge_2_7 = rge(_2, 3000, 5000, AMPLIFICATION);
        RegionGeneEvents rge_2_8 = rge(_2, 4000, 5000, HET_DELETION);

        var chr_1_1 = new ChromosomeRegionCounts(_1, List.of(rge_1_1, rge_1_2, rge_1_3, rge_1_4, rge_1_5, rge_1_6, rge_1_7, rge_1_8));
        var chr_2_1 = new ChromosomeRegionCounts(_2, List.of(rge_2_1, rge_2_2, rge_2_3, rge_2_4, rge_2_5, rge_2_6, rge_2_7, rge_2_8));
        TreeMap<HumanChromosome, ChromosomeRegionCounts> map1 = new TreeMap<>();
        map1.put(_1, chr_1_1);
        map1.put(_2, chr_2_1);
        EventCounts counts1 = new EventCounts(map1);

        var chr_1_2 = new ChromosomeRegionCounts(_1, List.of(rge_1_1, rge_1_3, rge_1_5, rge_1_7));
        var chr_2_2 = new ChromosomeRegionCounts(_2, List.of(rge_2_1, rge_2_3, rge_2_5, rge_2_7));
        TreeMap<HumanChromosome, ChromosomeRegionCounts> map2 = new TreeMap<>();
        map2.put(_1, chr_1_2);
        map2.put(_2, chr_2_2);
        EventCounts counts2 = new EventCounts(map2);

        var chr_1_3 = new ChromosomeRegionCounts(_1, List.of(rge_1_1, rge_1_3, rge_1_5, rge_1_7));
        var chr_2_3 = new ChromosomeRegionCounts(_2, List.of(rge_2_1, rge_2_3, rge_2_5, rge_2_7));
        TreeMap<HumanChromosome, ChromosomeRegionCounts> map3 = new TreeMap<>();
        map3.put(_1, chr_1_3);
        map3.put(_2, chr_2_3);
        EventCounts counts3 = new EventCounts(map3);

        var chr_1_4 = new ChromosomeRegionCounts(_1, List.of(rge_1_1, rge_1_3));
        var chr_2_4 = new ChromosomeRegionCounts(_2, List.of(rge_2_1, rge_2_3));
        TreeMap<HumanChromosome, ChromosomeRegionCounts> map4 = new TreeMap<>();
        map4.put(_1, chr_1_4);
        map4.put(_2, chr_2_4);
        EventCounts counts4 = new EventCounts(map4);

        var chr_1_5 = new ChromosomeRegionCounts(_1, List.of(rge_1_1));
        var chr_2_5 = new ChromosomeRegionCounts(_2, List.of(rge_2_1));
        TreeMap<HumanChromosome, ChromosomeRegionCounts> map5 = new TreeMap<>();
        map5.put(_1, chr_1_5);
        map5.put(_2, chr_2_5);
        EventCounts counts5 = new EventCounts(map5);

        EventCounts counts = new EventCounts(new TreeMap<>());
        counts.mergeAdd(counts1);
        counts.mergeAdd(counts2);
        counts.mergeAdd(counts3);
        counts.mergeAdd(counts4);
        counts.mergeAdd(counts5);

        // 1_1 -> 5, 1_3 -> 4, 1_5 -> 4 because 3 added directly and 1 is from 1_4, all others 2, 1.
        StringWriter sw = new StringWriter();
        counts.writeTo(new BufferedWriter(sw), 4, V38);
        String written = sw.toString();
        String[] lines = written.split("\n");
        assertEquals(7, lines.length);
        assertEquals("Chromosome,RegionStart,RegionEnd,Type,Frequency", lines[0]);
        assertEquals("chr1,1000,2000,AMP,5", lines[1]);
        assertEquals("chr1,2000,3000,AMP,4", lines[2]);
        assertEquals("chr1,2000,3000,DEL,4", lines[3]);
    }
}
