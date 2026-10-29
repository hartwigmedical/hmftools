package com.hartwig.hmftools.purple.tools;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._20;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._21;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._22;
import static com.hartwig.hmftools.purple.drivers.AmpDelRegionFrequency.EventType.AMP;
import static com.hartwig.hmftools.purple.drivers.AmpDelRegionFrequency.EventType.DEL;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;
import java.util.SortedMap;

import org.junit.Test;

public class SampleGermlineGeneTaskTest extends ToolsTestBase
{
    @Test
    public void runForSample1()
    {
        EventCounts counts = runForSamples(sample1);
        assertEquals(3, counts.data().size());

        ChromosomeRegionCounts counts20 = counts.data().get(_20);
        SortedMap<RegionAmpDel, Integer> data = counts20.counts();
        assertEquals(4, data.size());

        RegionAmpDel ampDel1 = rad(_20, 1914001, 1915000, DEL);// SIRPA
        int count = data.get(ampDel1);
        assertEquals(1, count);

        ampDel1 = rad(_20, 2692001, 2693000, AMP);// EBF4
        count = data.get(ampDel1);
        assertEquals(1, count);

        ampDel1 = rad(_20, 58692001, 58693000, DEL);// NPEPL1,STX16-NPEPL1
        count = data.get(ampDel1);
        assertEquals(1, count);

        ampDel1 = rad(_20, 58981001, 58982000, DEL);// NELFCD
        count = data.get(ampDel1);
        assertEquals(1, count);

        ChromosomeRegionCounts counts21 = counts.data().get(_21);
        data = counts21.counts();
        assertEquals(1, data.size());

        ampDel1 = rad(_21, 7929001, 8550172, DEL);// MIR6724-1
        count = data.get(ampDel1);
        assertEquals(1, count);

        ChromosomeRegionCounts counts22 = counts.data().get(_22);
        data = counts22.counts();
        assertEquals(1, data.size());

        ampDel1 = rad(_22, 10699413, 12694000, AMP);// ACTR3BP7
        count = data.get(ampDel1);
        assertEquals(1, count);
    }

    @Test
    public void runForSample3()
    {
        EventCounts counts = runForSamples(sample3);
        assertEquals(3, counts.data().size());

        ChromosomeRegionCounts counts20 = counts.data().get(_20);
        SortedMap<RegionAmpDel, Integer> data = counts20.counts();
        assertEquals(2, data.size());
        RegionAmpDel ampDel1 = rad(_20, 58692001, 58693000, DEL);// NPEPL1,STX16-NPEPL1
        int count = data.get(ampDel1);
        assertEquals(1, count);

        ampDel1 = rad(_20, 58981001, 58982000, DEL);// NELFCD
        count = data.get(ampDel1);
        assertEquals(1, count);

        ChromosomeRegionCounts counts21 = counts.data().get(_21);
        data = counts21.counts();
        assertEquals(1, data.size());

        ampDel1 = rad(_21, 7929001, 8550172, DEL);// MIR6724-1,MIR3648-1
        count = data.get(ampDel1);
        assertEquals(1, count);

        ChromosomeRegionCounts counts22 = counts.data().get(_22);
        data = counts22.counts();
        assertEquals(1, data.size());

        ampDel1 = rad(_22, 10699413, 12694000, AMP);// ACTR3BP7
        count = data.get(ampDel1);
        assertEquals(1, count);
    }

    @Test
    public void runFor2Samples()
    {
        EventCounts counts = runForSamples(sample1, sample3);
        assertEquals(3, counts.data().size());

        ChromosomeRegionCounts counts20 = counts.data().get(_20);
        SortedMap<RegionAmpDel, Integer> data = counts20.counts();
        assertEquals(4, data.size());

        RegionAmpDel ampDel1 = rad(_20, 1914001, 1915000, DEL);// SIRPA
        int count = data.get(ampDel1);
        assertEquals(1, count);

        ampDel1 = rad(_20, 2692001, 2693000, AMP);// EBF4
        count = data.get(ampDel1);
        assertEquals(1, count);

        ampDel1 = rad(_20, 58692001, 58693000, DEL);// NPEPL1,STX16-NPEPL1
        count = data.get(ampDel1);
        assertEquals(2, count);

        ampDel1 = rad(_20, 58981001, 58982000, DEL);// NELFCD
        count = data.get(ampDel1);
        assertEquals(2, count);

        ChromosomeRegionCounts counts21 = counts.data().get(_21);
        data = counts21.counts();
        assertEquals(1, data.size());

        ampDel1 = rad(_21, 7929001, 8550172, DEL);// MIR6724-1
        count = data.get(ampDel1);
        assertEquals(2, count);

        ChromosomeRegionCounts counts22 = counts.data().get(_22);
        data = counts22.counts();
        assertEquals(1, data.size());

        ampDel1 = rad(_22, 10699413, 12694000, AMP);// ACTR3BP7
        count = data.get(ampDel1);
        assertEquals(2, count);
    }

    private EventCounts runForSamples(String... samples)
    {
        SampleGermlineGeneTask task = new SampleGermlineGeneTask(1, mPurpleDataDir);
        Arrays.stream(samples).forEach(task.getSampleIds()::add);
        task.call();
        return task.getEventCounts();
    }
}
