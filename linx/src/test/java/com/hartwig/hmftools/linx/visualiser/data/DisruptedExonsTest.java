package com.hartwig.hmftools.linx.visualiser.data;

import static com.hartwig.hmftools.linx.visualiser.file.VisGeneAnnotationType.EXON_LOST;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.linx.visualiser.file.VisFusion;
import com.hartwig.hmftools.linx.visualiser.file.VisGeneExon;

import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

public class DisruptedExonsTest
{
    private static final String UP = "UP";
    private static final String DOWN = "DOWN";

    private List<VisGeneExon> normalUp = Lists.newArrayList();
    private List<VisGeneExon> reverseUp = Lists.newArrayList();
    private List<VisGeneExon> normalDown = Lists.newArrayList();
    private List<VisGeneExon> reverseDown = Lists.newArrayList();

    @Before
    public void setup()
    {
        normalUp.add(create(UP, 1, 100, 1));
        normalUp.add(create(UP, 201, 300, 2));

        normalDown.add(create(DOWN, 1, 100, 1));
        normalDown.add(create(DOWN, 201, 300, 2));


        reverseUp.add(create(UP, 1, 100, 2));
        reverseUp.add(create(UP, 201, 300, 1));

        reverseDown.add(create(DOWN, 1, 100, 2));
        reverseDown.add(create(DOWN, 201, 300, 1));
    }

    @Test
    public void testNormalOrientation()
    {
        final List<VisGeneExon> exons = Lists.newArrayList();
        exons.addAll(normalUp);
        exons.addAll(normalDown);

        final VisFusion fusion = create(1, 150, 1, 1, 150, 2);
        final List<GenomeRegion> victim = DisruptedExons.disruptedGeneRegions(fusion, exons);

        assertEquals(2, victim.size());
        assertRegion(UP, 100, 300, victim.get(0));
        assertRegion(DOWN, 1, 201, victim.get(1));
    }

    @Test
    public void testNormalOrientationPartialExon()
    {
        final List<VisGeneExon> exons = Lists.newArrayList();
        exons.addAll(normalUp);
        exons.addAll(normalDown);

        final VisFusion fusion = create(1, 50, 1, 1, 250, 2);
        final List<GenomeRegion> victim = DisruptedExons.disruptedGeneRegions(fusion, exons);

        assertEquals(2, victim.size());
        assertRegion(UP, 50, 300, victim.get(0));
        assertRegion(DOWN, 1, 250, victim.get(1));
    }

    @Test
    public void testReverseOrientation()
    {
        final List<VisGeneExon> exons = Lists.newArrayList();
        exons.addAll(reverseUp);
        exons.addAll(reverseDown);

        final VisFusion fusion = create(-1, 150, 1, -1, 150, 2);
        final List<GenomeRegion> victim = DisruptedExons.disruptedGeneRegions(fusion, exons);

        assertEquals(2, victim.size());
        assertRegion(UP, 1, 201, victim.get(0));
        assertRegion(DOWN, 100, 300, victim.get(1));
    }

    @Test
    public void testReverseOrientationPartialExon()
    {
        final List<VisGeneExon> exons = Lists.newArrayList();
        exons.addAll(reverseUp);
        exons.addAll(reverseDown);

        final VisFusion fusion = create(-1, 250, 1, -1, 50, 2);
        final List<GenomeRegion> victim = DisruptedExons.disruptedGeneRegions(fusion, exons);

        assertEquals(2, victim.size());
        assertRegion(UP, 1, 250, victim.get(0));
        assertRegion(DOWN, 50, 300, victim.get(1));
    }


    private static void assertRegion(String chr, long start, long end, GenomeRegion victim)
    {
        assertEquals(chr, victim.chromosome());
        assertEquals(start, victim.start());
        assertEquals(end, victim.end());
    }

    @NotNull
    static VisGeneExon create(String chr, int start, int end, int rank)
    {
        return new VisGeneExon("sample", 1, chr, chr, chr, EXON_LOST, rank, start, end);
    }

    @NotNull
    static VisFusion create(int strandUp, int positionUp, int fusedExonUp, int strandDown, int positionDown, int fusedExonDown)
    {
        return new VisFusion(
                "sample", 1, true, UP, UP, UP, positionUp, strandUp, UP, fusedExonUp,
                DOWN, DOWN, DOWN, positionDown, strandDown, DOWN, fusedExonDown);
    }

}
