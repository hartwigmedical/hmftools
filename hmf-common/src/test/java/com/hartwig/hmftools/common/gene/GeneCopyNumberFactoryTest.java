package com.hartwig.hmftools.common.gene;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.purple.copynumber.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.segment.StructuralVariantSupport;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import com.hartwig.hmftools.common.region.hmfslicer.ImmutableHmfGenomeRegion;

import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

public class GeneCopyNumberFactoryTest {

    private static final double EPSILON = 1e-10;

    private List<HmfGenomeRegion> genes;
    private List<PurpleCopyNumber> copyNumbers;

    @Before
    public void setup() {
        genes = Lists.newArrayList();
        genes.add(create("1", 300, 600, "First"));
        genes.add(create("1", 800, 1000, "Second"));
        genes.add(create("2", 200, 700, "Third"));

        copyNumbers = Lists.newArrayList();
        copyNumbers.add(createCopyNumber("1", 100, 500, 3));
        copyNumbers.add(createCopyNumber("1", 501, 700, 4));
        copyNumbers.add(createCopyNumber("1", 701, 799, 5));
        copyNumbers.add(createCopyNumber("1", 800, 1000, 6));
        copyNumbers.add(createCopyNumber("2", 1, 1000, 7));
    }

    @Test
    public void testGeneCopies() throws IOException, EmptyFileException {
        final List<GeneCopyNumber> geneCopies = GeneCopyNumberFactory.geneCopyNumbers(genes, copyNumbers);
        assertEquals(3, geneCopies.size());

        assertGeneCopy(geneCopies.get(0), "First", 3, 4, 2);
        assertGeneCopy(geneCopies.get(1), "Second", 6, 6, 1);
        assertGeneCopy(geneCopies.get(2), "Third", 7, 7, 1);
    }

    private void assertGeneCopy(GeneCopyNumber victim, String gene, double min, double max, int count) {
        assertEquals(victim.gene(), gene);
        assertEquals(victim.minCopyNumber(), min, EPSILON);
        assertEquals(victim.maxCopyNumber(), max, EPSILON);
        assertEquals(victim.regions(), count);
    }

    private static HmfGenomeRegion create(String chromosome, long start, long end, String name) {
        return new ImmutableHmfGenomeRegion(chromosome, start, end, name, 1, name, name, start, end, chromosome, name);
    }

    @NotNull
    private static PurpleCopyNumber createCopyNumber(String chromosome, long start, long end, double copyNumber) {
        return ImmutablePurpleCopyNumber.builder()
                .chromosome(chromosome)
                .start(start)
                .end(end)
                .averageTumorCopyNumber(copyNumber)
                .bafCount(0)
                .averageObservedBAF(0.5)
                .averageActualBAF(0.5)
                .ratioSupport(true)
                .structuralVariantSupport(StructuralVariantSupport.NONE)
                .build();
    }

}
