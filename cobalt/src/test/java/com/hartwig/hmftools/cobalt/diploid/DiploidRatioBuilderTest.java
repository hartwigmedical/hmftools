package com.hartwig.hmftools.cobalt.diploid;

import static junit.framework.TestCase.assertEquals;

import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.samtools.util.Locatable;

public class DiploidRatioBuilderTest {

    @Test
    public void testBuildRatios() {
        DiploidRatioBuilder victim = new DiploidRatioBuilder();
        victim.accept(locatable("1", 1001, 3000));
        victim.accept(locatable("1", 5001, 6000));
        victim.accept(locatable("2", 1001, 3000));

        ListMultimap<Chromosome, ReadRatio> result =  victim.build();
        assertEquals(5, result.size());
        assertReadRatio("1", 1001, result.get(HumanChromosome._1).get(0));
        assertReadRatio("1", 2001, result.get(HumanChromosome._1).get(1));
        assertReadRatio("1", 5001, result.get(HumanChromosome._1).get(2));
        assertReadRatio("2", 1001, result.get(HumanChromosome._2).get(0));
        assertReadRatio("2", 2001, result.get(HumanChromosome._2).get(1));
    }

    private void assertReadRatio(@NotNull String contig, long position, @NotNull ReadRatio victim) {
        assertEquals(contig, victim.chromosome());
        assertEquals(position, victim.position());
        assertEquals(1, victim.ratio(), 0.01);
    }

    @NotNull
    private static Locatable locatable(String contig, int start, int end) {
        return new Locatable() {
            @Override
            public String getContig() {
                return contig;
            }

            @Override
            public int getStart() {
                return start;
            }

            @Override
            public int getEnd() {
                return end;
            }
        };
    }
}
