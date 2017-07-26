package com.hartwig.hmftools.common.gene;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.purple.PurpleDatamodelTest;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
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
        copyNumbers = Lists.newArrayList();
    }

    @Test
    public void oldFreecTest() {
        genes.add(create("1", 101, 200, "First"));
        genes.add(create("1", 301, 400, "Second"));
        genes.add(create("1", 601, 700, "Third"));

        copyNumbers.add(createCopyNumber("1", 11, 20, 1));
        copyNumbers.add(createCopyNumber("1", 21, 120, 2));
        copyNumbers.add(createCopyNumber("1", 121, 170, 4));
        copyNumbers.add(createCopyNumber("1", 171, 190, 2));
        copyNumbers.add(createCopyNumber("1", 191, 260, 5));
        copyNumbers.add(createCopyNumber("1", 261, 290, 2));
        copyNumbers.add(createCopyNumber("1", 291, 500, 0));
        copyNumbers.add(createCopyNumber("1", 501, 700, 2));

        final List<GeneCopyNumber> geneCopies = GeneCopyNumberFactory.geneCopyNumbers(genes, copyNumbers);
        assertEquals(3, geneCopies.size());

        assertGeneCopy(geneCopies.get(0), "First", 2, 5, 3.3);
        assertGeneCopy(geneCopies.get(1), "Second", 0, 0, 0.0);
        assertGeneCopy(geneCopies.get(2), "Third", 2, 2, 2.0);
    }

    @Test
    public void testGeneCopies() throws IOException, EmptyFileException {

        genes.add(create("1", 300, 600, "First"));
        genes.add(create("1", 800, 1000, "Second"));
        genes.add(create("2", 200, 700, "Third"));

        copyNumbers.add(createCopyNumber("1", 100, 500, 3));
        copyNumbers.add(createCopyNumber("1", 501, 700, 4));
        copyNumbers.add(createCopyNumber("1", 701, 799, 5));
        copyNumbers.add(createCopyNumber("1", 800, 1000, 6));
        copyNumbers.add(createCopyNumber("2", 1, 1000, 7));

        final List<GeneCopyNumber> geneCopies = GeneCopyNumberFactory.geneCopyNumbers(genes, copyNumbers);
        assertEquals(3, geneCopies.size());

        assertGeneCopy(geneCopies.get(0), "First", 3, 4, 2);
        assertGeneCopy(geneCopies.get(1), "Second", 6, 6, 1);
        assertGeneCopy(geneCopies.get(2), "Third", 7, 7, 1);
    }

    private void assertGeneCopy(GeneCopyNumber victim, String gene, double min, double max, int count) {
        assertEquals(gene, victim.gene());
        assertEquals(min, victim.minCopyNumber(), EPSILON);
        assertEquals(max, victim.maxCopyNumber(), EPSILON);
        assertEquals(count, victim.regions());
    }

    private void assertGeneCopy(GeneCopyNumber victim, String gene, double min, double max, double mean) {
        assertEquals(gene, victim.gene());
        assertEquals(min, victim.minCopyNumber(), EPSILON);
        assertEquals(max, victim.maxCopyNumber(), EPSILON);
        assertEquals(mean, victim.meanCopyNumber(), EPSILON);
    }

    private static HmfGenomeRegion create(String chromosome, long start, long end, String name) {
        return new ImmutableHmfGenomeRegion(chromosome, start, end, name, 1, name, name, start, end, chromosome, name);
    }

    @NotNull
    private static PurpleCopyNumber createCopyNumber(String chromosome, long start, long end, double copyNumber) {
        return PurpleDatamodelTest.createCopyNumber(chromosome, start, end, copyNumber).build();
    }

}
