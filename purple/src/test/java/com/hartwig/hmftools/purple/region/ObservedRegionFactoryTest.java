package com.hartwig.hmftools.purple.region;

import static com.hartwig.hmftools.purple.FittingTestUtils.buildCobaltChromosomes;
import static com.hartwig.hmftools.purple.FittingTestUtils.createCobaltRatio;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.purple.segment.PurpleSupportSegment;

import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

public class ObservedRegionFactoryTest
{
    private final int segmentLength = 1000;
    private final Chromosome chrA = HumanChromosome._1;
    private final Chromosome chrB = HumanChromosome._2;
    ObservedRegionFactory factory = new ObservedRegionFactory(segmentLength, buildCobaltChromosomes());
    private final List<AmberBAF> allBafs = new LinkedList<>();
    private final List<CobaltRatio> allCobaltRatios = new ArrayList<>();

    @Before
    public void setup()
    {
        allBafs.clear();
        allCobaltRatios.clear();
    }

    @Test
    public void noSupportSegments()
    {
        Multimap<Chromosome, AmberBAF> amberMap = ArrayListMultimap.create();
        amberMap.put(chrA, baf(5000));
        Map<Chromosome, List<CobaltRatio>> ratios = new HashMap<>();
        ratios.put(chrA, List.of(createCobaltRatio(chrA, 6000, 0.6, 0.5)));
        List<ObservedRegion> lor = factory.formObservedRegions(List.of(), amberMap, ratios);
        assertTrue(lor.isEmpty());
    }

    @Test
    public void measurementsOutsideSupportSegment()
    {
        Map<Chromosome, List<CobaltRatio>> ratios = new HashMap<>();
        ratios.put(chrA, List.of(createCobaltRatio(chrA, 5000, 0.6, 0.7)));
        List<ObservedRegion> lor = factory.formObservedRegions(List.of(purpleSupportSegment(10_000)), bafMap(), ratios);
        assertEquals(1, lor.size());
        assertEquals(0, lor.get(0).bafCount());
        assertEquals(0, lor.get(0).observedTumorRatio(), 0.001);
        assertEquals(0.0, lor.get(0).gcContent(), 0.001);
    }

    @Test
    public void measurementsInsideSupportSegment()
    {
        Map<Chromosome, List<CobaltRatio>> ratios = new HashMap<>();
        ratios.put(chrA, List.of(createCobaltRatio(chrA, 5001, 0.6, 0.7)));
        PurpleSupportSegment segment = purpleSupportSegment(5000);
        List<ObservedRegion> lor = factory.formObservedRegions(List.of(segment), bafMap(), ratios);
        assertEquals(1, lor.size());
        checkRegion(lor.get(0), segment, 1, 0.6, 0.7);
    }

    @Test
    public void measurementsOnDifferentChromosome()
    {
        Map<Chromosome, List<CobaltRatio>> ratios = new HashMap<>();
        ratios.put(chrA, List.of(createCobaltRatio(chrA, 5001, 0.6, 0.7)));
        PurpleSupportSegment segment = purpleSupportSegment(chrB, 5000);
        List<ObservedRegion> lor = factory.formObservedRegions(List.of(segment), bafMap(), ratios);
        assertEquals(1, lor.size());
        checkRegion(lor.get(0), segment, 0, 0.0, 0.0);
    }

    @Test
    public void dataAreAssignedToSegments()
    {
        tumorOnlyBAF(5001, 0.40, 40);
        tumorOnlyBAF(5101, 0.41, 41);
        tumorOnlyBAF(5201, 0.42, 42);
        tumorOnlyBAF(15001, 0.50, 50);
        tumorOnlyBAF(15101, 0.51, 51);
        tumorOnlyBAF(15201, 0.52, 52);
        tumorOnlyBAF(25001, 0.60, 60);
        tumorOnlyBAF(25101, 0.61, 61);
        tumorOnlyBAF(25201, 0.62, 62);
        Multimap<Chromosome, AmberBAF> amberMap = ArrayListMultimap.create();
        allBafs.forEach(amberBAF -> amberMap.put(amberBAF.chr(), amberBAF));

        cobaltRatio(5002, 0.40, 0.44);
        cobaltRatio(5102, 0.41, 0.46);
        cobaltRatio(15002, 0.50, 0.54);
        cobaltRatio(15102, 0.51, 0.56);
        cobaltRatio(25002, 0.60, 0.64);
        cobaltRatio(25102, 0.61, 0.66);
        Map<Chromosome, List<CobaltRatio>> cobaltRatios = new HashMap<>();
        cobaltRatios.put(chrA, allCobaltRatios);

        PurpleSupportSegment segmentA = purpleSupportSegment(5000);
        PurpleSupportSegment segmentB = purpleSupportSegment(15000);
        PurpleSupportSegment segmentC = purpleSupportSegment(25000);

        List<ObservedRegion> lor = factory.formObservedRegions(List.of(segmentA, segmentB, segmentC), amberMap, cobaltRatios);
        assertEquals(3, lor.size());
        checkRegion(lor.get(0), segmentA, 3, 0.405, 0.45);
        checkRegion(lor.get(1), segmentB, 3, 0.505, 0.55);
        checkRegion(lor.get(2), segmentC, 3, 0.605, 0.65);
    }

    @Test
    public void negativeGcRatiosAreIgnored()
    {
        Map<Chromosome, List<CobaltRatio>> ratios = new HashMap<>();
        ratios.put(chrA, List.of(createCobaltRatio(chrA, 5000, -1.0, 0.4)));
        PurpleSupportSegment segment = purpleSupportSegment(4900);
        List<ObservedRegion> lor = factory.formObservedRegions(List.of(segment), bafMap(), ratios);
        assertEquals(1, lor.size());
        assertEquals(0.0, lor.get(0).observedTumorRatio(), 0.001);
        assertEquals(0.0, lor.get(0).gcContent(), 0.001); // If the ratio is -1 then the gc content value is not used.
    }

    @Test
    public void windowsWithNegativeReferenceDiploidGcRatiosAreIgnored()
    {
        Map<Chromosome, List<CobaltRatio>> ratios = new HashMap<>();
        CobaltRatio ratio = new CobaltRatio(chrA.toString(), 5000, 0.5, 0.4, 0.4, -1.0, 0.4, 0.5, 0.5);
        ratios.put(chrA, List.of(ratio));
        PurpleSupportSegment segment = purpleSupportSegment(4900);
        List<ObservedRegion> lor = factory.formObservedRegions(List.of(segment), bafMap(), ratios);
        assertEquals(1, lor.size());
        assertEquals(0.0, lor.get(0).observedTumorRatio(), 0.001);
        assertEquals(0.0, lor.get(0).gcContent(), 0.001);
    }

    @Test
    public void germlineStatus()
    {
        tumorOnlyBAF(5001, 0.4, 40);
        tumorOnlyBAF(15001, 0.5, 50);
        tumorOnlyBAF(25001, 0.6, 60);
        Multimap<Chromosome, AmberBAF> amberMap = ArrayListMultimap.create();
        allBafs.forEach(amberBAF -> amberMap.put(amberBAF.chr(), amberBAF));

        cobaltRatio(5002, 0.40, 0.4, 0.0);
        cobaltRatio(15002, 0.05, 0.5, 0.09);
        cobaltRatio(25102, 0.60, 0.6, 1.0);
        Map<Chromosome, List<CobaltRatio>> cobaltRatios = new HashMap<>();
        cobaltRatios.put(chrA, allCobaltRatios);

        PurpleSupportSegment segmentA = purpleSupportSegment(5000);
        PurpleSupportSegment segmentB = purpleSupportSegment(15000);
        PurpleSupportSegment segmentC = purpleSupportSegment(25000);

        List<ObservedRegion> lor = factory.formObservedRegions(List.of(segmentA, segmentB, segmentC), amberMap, cobaltRatios);
        assertEquals(3, lor.size());
        assertEquals(GermlineStatus.UNKNOWN, lor.get(0).germlineStatus());
        assertEquals(GermlineStatus.HOM_DELETION, lor.get(1).germlineStatus());
        assertEquals(GermlineStatus.DIPLOID, lor.get(2).germlineStatus());
    }

    @Test
    public void gcContentIsMeanAndGcRatioIsMedian()
    {
        tumorOnlyBAF(5001, 0.40, 40);
        tumorOnlyBAF(5101, 0.41, 41);
        tumorOnlyBAF(5201, 0.42, 42);
        Multimap<Chromosome, AmberBAF> amberMap = ArrayListMultimap.create();
        allBafs.forEach(amberBAF -> amberMap.put(amberBAF.chr(), amberBAF));

        cobaltRatio(5002, 0.3, 0.3);
        cobaltRatio(5102, 0.6, 0.6);
        cobaltRatio(5202, 0.6, 0.6);
        Map<Chromosome, List<CobaltRatio>> cobaltRatios = new HashMap<>();
        cobaltRatios.put(chrA, allCobaltRatios);

        PurpleSupportSegment segmentA = purpleSupportSegment(5000);

        List<ObservedRegion> lor = factory.formObservedRegions(List.of(segmentA), amberMap, cobaltRatios);
        assertEquals(1, lor.size());
        checkRegion(lor.get(0), segmentA, 3, 0.6, 0.5);
    }

    @Test
    public void minSupportExtension()
    {
        tumorOnlyBAF(5001, 0.4, 40);
        tumorOnlyBAF(15001, 0.5, 50);
        tumorOnlyBAF(25001, 0.6, 60);
        Multimap<Chromosome, AmberBAF> amberMap = ArrayListMultimap.create();
        allBafs.forEach(amberBAF -> amberMap.put(amberBAF.chr(), amberBAF));

        cobaltRatio(5002, 0.40, 0.4, 0.0);
        cobaltRatio(15002, 0.05, 0.5, 0.0);
        cobaltRatio(25102, 0.60, 0.6, 1.0);
        Map<Chromosome, List<CobaltRatio>> cobaltRatios = new HashMap<>();
        cobaltRatios.put(chrA, allCobaltRatios);

        PurpleSupportSegment segmentA = purpleSupportSegment(5000);
        PurpleSupportSegment segmentB = purpleSupportSegment(15000);
        PurpleSupportSegment segmentC = purpleSupportSegment(25000);

        List<ObservedRegion> lor = factory.formObservedRegions(List.of(segmentA, segmentB, segmentC), amberMap, cobaltRatios);
        assertEquals(3, lor.size());
        assertEquals(GermlineStatus.UNKNOWN, lor.get(0).germlineStatus());
        assertEquals(GermlineStatus.UNKNOWN, lor.get(1).germlineStatus());
        assertEquals(GermlineStatus.DIPLOID, lor.get(2).germlineStatus());
        assertEquals(5000, lor.get(2).minStart());
    }

    @NotNull
    private Multimap<Chromosome, AmberBAF> bafMap()
    {
        Multimap<Chromosome, AmberBAF> amberMap = ArrayListMultimap.create();
        amberMap.put(chrA, baf(5001));
        return amberMap;
    }

    private void cobaltRatio(int position, double tumorRatio, double tumorGcContent)
    {
        CobaltRatio cobaltRatio = createCobaltRatio(chrA, position, tumorRatio, tumorGcContent);
        allCobaltRatios.add(cobaltRatio);
    }

    private void cobaltRatio(int position, double tumorRatio, double tumorGcContent, double referenceGcDiploidRatio)
    {
        CobaltRatio cobaltRatio = new CobaltRatio(
                chrA.toString(),
                position,
                0,
                0,
                1,
                referenceGcDiploidRatio,
                0.5,
                tumorRatio,
                tumorGcContent);
        allCobaltRatios.add(cobaltRatio);
    }

    private void tumorOnlyBAF(int position, double tumorBAF, int tumorDepth)
    {
        AmberBAF baf = new AmberBAF(chrA.toString(), position, tumorBAF, tumorDepth, -1, 0);
        allBafs.add(baf);
    }

    private AmberBAF baf(int position)
    {
        return new AmberBAF(chrA.toString(), position, 0.5, 50, -1, 0);
    }

    private void checkRegion(ObservedRegion observedRegion, PurpleSupportSegment segment, int bafCount, double observedTumorRatio,
            double gcContent)
    {
        assertEquals(segment.Start, observedRegion.start());
        assertEquals(segment.End, observedRegion.end());
        assertEquals(bafCount, observedRegion.bafCount());
        assertEquals(observedTumorRatio, observedRegion.observedTumorRatio(), 0.001);
        assertEquals(gcContent, observedRegion.gcContent(), 0.001);

    }

    private PurpleSupportSegment purpleSupportSegment(int start)
    {
        return purpleSupportSegment(chrA, start);
    }

    private PurpleSupportSegment purpleSupportSegment(Chromosome chr, int start)
    {
        return new PurpleSupportSegment(chr.toString(), start, start + 1000, true, SegmentSupport.NONE, false, start, start);
    }
}
