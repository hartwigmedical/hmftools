package com.hartwig.hmftools.purple.data;

import static org.junit.Assert.assertEquals;

import java.nio.file.Paths;
import java.util.List;
import java.util.Set;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.utils.pcf.PCFPosition;
import com.hartwig.hmftools.common.utils.pcf.PCFSource;

import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

public class PcfDataTest
{
    private final String CobaltDirectory = Paths.get("src", "test", "resources", "cobalt").toAbsolutePath().toString();
    private final PcfData data = new PcfData(CobaltDirectory);

//    @Test
    public void comparison() throws Exception
    {
        String baselineDir = "/Users/timlavers/work/junk/";
        ListMultimap<Chromosome, PCFPosition> baseline = new PcfData(baselineDir).loadCobaltSegments("COLO829v003T", PCFSource.TUMOR_RATIO);

        String newVersionDir = "/Users/timlavers/work/batches/2025/11/13/1/cobalt/";
        ListMultimap<Chromosome, PCFPosition> newVersion = new PcfData(newVersionDir).loadCobaltSegments("COLO829v003T", PCFSource.TUMOR_RATIO);

        var chr1Baseline = baseline.get(HumanChromosome._1);
        var newChr1Baseline = newVersion.get(HumanChromosome._1);
        System.out.println(chr1Baseline.size());
        System.out.println(newChr1Baseline.size());
        Set<PCFPosition> baselineSet = Sets.newHashSet(baseline.get(HumanChromosome._1));
        Set<PCFPosition> newVersionSet = Sets.newHashSet(newVersion.get(HumanChromosome._1));

//        Sets.difference(baselineSet, newVersionSet).forEach(System.out::println);
        baselineSet.forEach(position -> {
            if (!newVersionSet.contains(position)) {
                System.out.println("In old, not in new: " + position);
            }
        });
        newVersionSet.forEach(position -> {
            if (!baselineSet.contains(position)) {
                System.out.println("In new, not in old: " + position);
            }
        });
//        int limit = Math.min(chr1Baseline.size(), newChr1Baseline.size());
//        for (int i = 72; i < limit; i++)
//        {
//            System.out.println(i);
//            System.out.println("old: " + chr1Baseline.get(i));
//            System.out.println("new: " + newChr1Baseline.get(i));
//            assertEquals(chr1Baseline.get(i).Position, newChr1Baseline.get(i).Position);
//            assertEquals(chr1Baseline.get(i).minPosition(), newChr1Baseline.get(i).minPosition());
//            assertEquals(chr1Baseline.get(i).maxPosition(), newChr1Baseline.get(i).maxPosition());
//        }

    }

    @Test
    public void returnEmptyForNull() throws Exception
    {
        Assert.assertTrue(data.loadCobaltSegments(null, PCFSource.REFERENCE_RATIO).isEmpty());
    }

    @Test
    public void oneSegment() throws Exception
    {
        ListMultimap<Chromosome, PCFPosition> loaded = data.loadCobaltSegments("s1", PCFSource.REFERENCE_RATIO);
        assertEquals(1, loaded.keySet().size());
        List<PCFPosition> positions = loaded.get(HumanChromosome._1);
        assertEquals(2, positions.size());
        assertEquals(1001, positions.get(0).position());
        assertEquals(1, positions.get(0).minPosition());
        assertEquals(1001, positions.get(0).maxPosition());

        assertEquals(4001, positions.get(1).position());
        assertEquals(4001, positions.get(1).minPosition());
        assertEquals(4001, positions.get(1).maxPosition());
        assertEquals(PCFSource.REFERENCE_RATIO, positions.get(0).Source);
    }

    @Test
    public void nonAbuttingSegments() throws Exception
    {
        ListMultimap<Chromosome, PCFPosition> loaded = data.loadCobaltSegments("s2", PCFSource.REFERENCE_RATIO);
        assertEquals(1, loaded.keySet().size());
        List<PCFPosition> positions = loaded.get(HumanChromosome._1);
        assertEquals(6, positions.size());
        assertEquals(1001, positions.get(0).position());
        assertEquals(1, positions.get(0).minPosition());
        assertEquals(1001, positions.get(0).maxPosition());

        assertEquals(2001, positions.get(1).position());
        assertEquals(2001, positions.get(1).minPosition());
        assertEquals(4001, positions.get(1).maxPosition());

        assertEquals(4001, positions.get(2).position());
        assertEquals(2001, positions.get(2).minPosition());
        assertEquals(4001, positions.get(2).maxPosition());

        assertEquals(6001, positions.get(3).position());
        assertEquals(6001, positions.get(3).minPosition());
        assertEquals(9001, positions.get(3).maxPosition());

        assertEquals(9001, positions.get(4).position());
        assertEquals(6001, positions.get(4).minPosition());
        assertEquals(9001, positions.get(4).maxPosition());

        assertEquals(12001, positions.get(5).position());
        assertEquals(12001, positions.get(5).minPosition());
        assertEquals(12001, positions.get(5).maxPosition());
    }

    @Test
    public void abuttingSegments() throws Exception
    {
        ListMultimap<Chromosome, PCFPosition> loaded = data.loadCobaltSegments("s3", PCFSource.REFERENCE_RATIO);
        assertEquals(1, loaded.keySet().size());
        List<PCFPosition> positions = loaded.get(HumanChromosome._1);
        assertEquals(4, positions.size());
        assertEquals(1001, positions.get(0).position());
        assertEquals(1, positions.get(0).minPosition());
        assertEquals(1001, positions.get(0).maxPosition());

        assertEquals(2001, positions.get(1).position());
        assertEquals(2001, positions.get(1).minPosition());
        assertEquals(2001, positions.get(1).maxPosition());

        assertEquals(4001, positions.get(2).position());
        assertEquals(4001, positions.get(2).minPosition());
        assertEquals(4001, positions.get(2).maxPosition());

        assertEquals(6001, positions.get(3).position());
        assertEquals(6001, positions.get(3).minPosition());
        assertEquals(6001, positions.get(3).maxPosition());
    }

    @Test
    public void multipleChromosomes() throws Exception
    {
        ListMultimap<Chromosome, PCFPosition> loaded = data.loadCobaltSegments("s4", PCFSource.REFERENCE_RATIO);
        assertEquals(5, loaded.keySet().size());
        List<PCFPosition> positions1 = loaded.get(HumanChromosome._1);

        assertEquals(5, positions1.size());
        assertEquals(820001, positions1.get(0).position());
        assertEquals(1, positions1.get(0).minPosition());
        assertEquals(820001, positions1.get(0).maxPosition());

        assertEquals(1979001, positions1.get(1).position());
        assertEquals(1979001, positions1.get(1).minPosition());
        assertEquals(1979001, positions1.get(1).maxPosition());

        assertEquals(1982001, positions1.get(2).position());
        assertEquals(1982001, positions1.get(2).minPosition());
        assertEquals(1982001, positions1.get(2).maxPosition());

        assertEquals(3725001, positions1.get(3).position());
        assertEquals(3725001, positions1.get(3).minPosition());
        assertEquals(3725001, positions1.get(3).maxPosition());

        assertEquals(3727001, positions1.get(4).position());
        assertEquals(3727001, positions1.get(4).minPosition());
        assertEquals(3727001, positions1.get(4).maxPosition());

        List<PCFPosition> positionsX = loaded.get(HumanChromosome._X);
        assertEquals(4, positionsX.size());
    }

    @Test
    public void handleOldFormat() throws Exception
    {
        ListMultimap<Chromosome, PCFPosition> loaded = data.loadCobaltSegments("s5",  PCFSource.TUMOR_RATIO);
        assertEquals(2, loaded.keySet().size());
        List<PCFPosition> positions1 = loaded.get(HumanChromosome._1);

        assertEquals(6, positions1.size());
        assertEquals(820001, positions1.get(0).position());
        assertEquals(1, positions1.get(0).minPosition());
        assertEquals(820001, positions1.get(0).maxPosition());
        assertEquals(PCFSource.TUMOR_RATIO, positions1.get(0).Source);
    }

    @Test(expected = org.apache.commons.cli.ParseException.class)
    public void handleMissingFile() throws Exception
    {
        data.loadCobaltSegments("missing", PCFSource.TUMOR_RATIO);
    }
}