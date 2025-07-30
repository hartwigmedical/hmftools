package com.hartwig.hmftools.purple.segment;

import static java.util.Collections.emptyList;

import static com.hartwig.hmftools.common.purple.PurpleTestUtils.createStructuralVariant;
import static com.hartwig.hmftools.common.purple.SegmentSupport.CENTROMERE;
import static com.hartwig.hmftools.common.purple.SegmentSupport.NONE;
import static com.hartwig.hmftools.common.purple.SegmentSupport.TELOMERE;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.purple.MiscTestUtils.gcProfile;

import static org.immutables.value.internal.$guava$.collect.$ImmutableList.of;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.PurpleTestUtils;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.utils.pcf.PCFPosition;
import com.hartwig.hmftools.common.utils.pcf.PCFSource;
import com.hartwig.hmftools.purple.AmberData;
import com.hartwig.hmftools.purple.CobaltData;
import com.hartwig.hmftools.purple.region.ObservedRegion;

import org.apache.commons.lang3.RandomUtils;
import org.junit.Before;
import org.junit.Test;
import org.mockito.Mockito;

public class SegmentationTest
{
    private static final int SV_LENGTH = 10_000;
    private final Chromosome chr1 = HumanChromosome._1;
    private final Chromosome chr2 = HumanChromosome._2;
    private final GenomePosition CHR1_LENGTH = GenomePositions.create(chr1.toString(), 2_000_000);
    private final GenomePosition CHR1_CENTROMERE = GenomePositions.create(chr1.toString(), 1_000_000);
    private final GenomePosition CHR2_LENGTH = GenomePositions.create(chr2.toString(), 1_000_000);
    private final GenomePosition CHR2_CENTROMERE = GenomePositions.create(chr2.toString(), 500_000);
    private CobaltChromosomes cobaltChromosomes;
    private Segmentation segmentation;

    @Before
    public void setUp()
    {
        Map<Chromosome, GenomePosition> centromeres = new HashMap<>();
        centromeres.put(chr1, CHR1_CENTROMERE);
        centromeres.put(chr2, CHR2_CENTROMERE);
        Map<Chromosome, GenomePosition> chromosomeLengths = new HashMap<>();
        chromosomeLengths.put(chr1, CHR1_LENGTH);
        chromosomeLengths.put(chr2, CHR2_LENGTH);
        Multimap<Chromosome, GCProfile> gcProfiles = HashMultimap.create();
        gcProfiles.put(chr1, gcProfile(chr1, 1, 500_000, 0.95, 0.53));
        gcProfiles.put(chr1, gcProfile(chr1, 500_001, 500_000, 0.96, 0.55));
        gcProfiles.put(chr1, gcProfile(chr1, 1_000_001, 500_000, 0.97, 0.57));
        gcProfiles.put(chr1, gcProfile(chr1, 1_500_001, 500_000, 0.98, 0.59));
        gcProfiles.put(chr2, gcProfile(chr2, 1, 500_000, 0.90, 0.60));
        gcProfiles.put(chr2, gcProfile(chr2, 500_001, 500_000, 0.91, 0.62));
        SegmentationReferenceData referenceData = Mockito.mock(SegmentationReferenceData.class);
        Mockito.when(referenceData.centromeres()).thenReturn(centromeres);
        Mockito.when(referenceData.chromosomeLengths()).thenReturn(chromosomeLengths);

        Set<MedianRatio> ratios = new HashSet<>();
        ratios.add(new MedianRatio(chr1.toString(), 0.6, 60));
        ratios.add(new MedianRatio(chr2.toString(), 0.5, 50));
        cobaltChromosomes = new CobaltChromosomes(ratios);

        //        AmberData amberData
        segmentation = new Segmentation(referenceData);
    }

    @Test
    public void noStructuralVariantsNoAmberDataNoCobaltData()
    {
        List<ObservedRegion> segments = segmentation.createObservedRegions(emptyList(), emptyAmberData(), emptyCobaltData());
        assertEquals(0, segments.size());
    }

    @Test
    public void oneStructuralVariantNoAmberChromosomesNoCobaltData()
    {
        StructuralVariant sv = sv(chr1, 100_000);
        List<ObservedRegion> segments = segmentation.createObservedRegions(of(sv), emptyAmberData(), emptyCobaltData());
        // The list of chromosomes in a key internal calculation comes from the amber data, which is empty.
        assertEquals(0, segments.size());
    }

    @Test
    public void splitByChromosomeArms()
    {
        List<ObservedRegion> segments = segmentation.createObservedRegions(of(), amberPcfPerChromosomeArm(), emptyCobaltData());
        assertEquals(4, segments.size());

        assertEquals(1, segments.get(0).start());
        assertEquals(CHR1_CENTROMERE.position() - 1, segments.get(0).end());
        assertEquals(chr1.toString(), segments.get(0).chromosome());
        assertEquals(TELOMERE, segments.get(0).support());

        assertEquals(CHR1_CENTROMERE.position(), segments.get(1).start());
        assertEquals(CHR1_LENGTH.position(), segments.get(1).end());
        assertEquals(chr1.toString(), segments.get(1).chromosome());
        assertEquals(CENTROMERE, segments.get(1).support());

        assertEquals(1, segments.get(2).start());
        assertEquals(CHR2_CENTROMERE.position() - 1, segments.get(2).end());
        assertEquals(chr2.toString(), segments.get(2).chromosome());
        assertEquals(TELOMERE, segments.get(2).support());

        assertEquals(CHR2_CENTROMERE.position(), segments.get(3).start());
        assertEquals(CHR2_LENGTH.position(), segments.get(3).end());
        assertEquals(chr2.toString(), segments.get(3).chromosome());
        assertEquals(CENTROMERE, segments.get(3).support());
    }

    @Test
    public void amberPcfRegionsWithinAChromosomeArm()
    {
        int amberPos1 = 100_000;
        int amberPos2 = 800_000;
        AmberData amberData = amberPcfPerChromosomeArm();
        amberData.TumorSegments.put(chr1, new PCFPosition(PCFSource.TUMOR_BAF, chr1.toString(), amberPos1));
        amberData.TumorSegments.put(chr1, new PCFPosition(PCFSource.TUMOR_BAF, chr1.toString(), amberPos2));

        List<ObservedRegion> segments = segmentation.createObservedRegions(of(), amberData, emptyCobaltData());
        assertEquals(6, segments.size());

        assertEquals(1, segments.get(0).start());
        assertEquals(amberPos1 - 1, segments.get(0).end());
        assertEquals(TELOMERE, segments.get(0).support());

        assertEquals(amberPos1, segments.get(1).start());
        assertEquals(amberPos2 - 1, segments.get(1).end());
        assertEquals(NONE, segments.get(1).support());

        assertEquals(amberPos2, segments.get(2).start());
        assertEquals(CHR1_CENTROMERE.position() - 1, segments.get(2).end());
        assertEquals(NONE, segments.get(2).support());
    }

    @Test
    public void inversionsInAChromosomeArm()
    {
        int sv1Start = 100_000;
        int sv2Start = 200_000;

        StructuralVariant sv1 = sv(chr1, sv1Start, INV);
        StructuralVariant sv2 = sv(chr1, sv2Start, INV);
        List<ObservedRegion> segments = segmentation.createObservedRegions(of(sv1, sv2), amberPcfPerChromosomeArm(), emptyCobaltData());
        assertEquals(8, segments.size());

        assertEquals(1, segments.get(0).start());
        assertEquals(sv1Start, segments.get(0).end());
        assertEquals(chr1.toString(), segments.get(0).chromosome());
        assertEquals(TELOMERE, segments.get(0).support());

        assertEquals(sv1Start + 1, segments.get(1).start());
        assertEquals(sv1Start + SV_LENGTH, segments.get(1).end());
        assertEquals(chr1.toString(), segments.get(1).chromosome());
        assertEquals(SegmentSupport.INV, segments.get(1).support());

        assertEquals(sv1Start + SV_LENGTH + 1, segments.get(2).start());
        assertEquals(sv2Start, segments.get(2).end());
        assertEquals(chr1.toString(), segments.get(2).chromosome());
        assertEquals(SegmentSupport.INV, segments.get(2).support());

        assertEquals(sv2Start + 1, segments.get(3).start());
        assertEquals(sv2Start + SV_LENGTH, segments.get(3).end());
        assertEquals(chr1.toString(), segments.get(3).chromosome());
        assertEquals(SegmentSupport.INV, segments.get(3).support());

        assertEquals(sv2Start + SV_LENGTH + 1, segments.get(4).start());
        assertEquals(CHR1_CENTROMERE.position() - 1, segments.get(4).end());
        assertEquals(chr1.toString(), segments.get(4).chromosome());
        assertEquals(SegmentSupport.INV, segments.get(4).support());
    }

    @Test
    public void deletionsInAChromosomeArm()
    {
        int del1Start = 100_000;
        int del2Start = 200_000;

        StructuralVariant sv1 = sv(chr1, del1Start, DEL);
        StructuralVariant sv2 = sv(chr1, del2Start, DEL);
        List<ObservedRegion> segments = segmentation.createObservedRegions(of(sv1, sv2), amberPcfPerChromosomeArm(), emptyCobaltData());
        assertEquals(8, segments.size());

        assertEquals(1, segments.get(0).start());
        assertEquals(del1Start, segments.get(0).end());
        assertEquals(TELOMERE, segments.get(0).support());

        assertEquals(del1Start + 1, segments.get(1).start());
        assertEquals(del1Start + SV_LENGTH - 1, segments.get(1).end());
        assertEquals(SegmentSupport.DEL, segments.get(1).support());

        assertEquals(del1Start + SV_LENGTH, segments.get(2).start());
        assertEquals(del2Start, segments.get(2).end());
        assertEquals(SegmentSupport.DEL, segments.get(2).support());

        assertEquals(del2Start + 1, segments.get(3).start());
        assertEquals(del2Start + SV_LENGTH - 1, segments.get(3).end());
        assertEquals(SegmentSupport.DEL, segments.get(3).support());

        assertEquals(del2Start + SV_LENGTH, segments.get(4).start());
        assertEquals(CHR1_CENTROMERE.position() - 1, segments.get(4).end());
        assertEquals(SegmentSupport.DEL, segments.get(4).support());
    }

    @Test
    public void splitByCobaltRegionsWithinChromosomeArm()
    {
        int cobaltTumorPos = 100_000;
        int cobaltReferencePos = 200_000;

        CobaltData cobaltData = new CobaltData(cobaltChromosomes);
        cobaltData.TumorSegments.put(chr1, new PCFPosition(PCFSource.TUMOR_RATIO, chr1.toString(), cobaltTumorPos));
        cobaltData.TumorSegments.put(chr1, new PCFPosition(PCFSource.REFERENCE_RATIO, chr1.toString(), cobaltReferencePos));

        List<ObservedRegion> segments = segmentation.createObservedRegions(of(), amberPcfPerChromosomeArm(), cobaltData);
        assertEquals(6, segments.size());

        assertEquals(1, segments.get(0).start());
        assertEquals(cobaltTumorPos - 1, segments.get(0).end());
        assertEquals(chr1.toString(), segments.get(0).chromosome());
        assertEquals(TELOMERE, segments.get(0).support());

        assertEquals(cobaltTumorPos, segments.get(1).start());
        assertEquals(cobaltReferencePos - 1, segments.get(1).end());
        assertEquals(chr1.toString(), segments.get(1).chromosome());
        assertEquals(SegmentSupport.NONE, segments.get(1).support());

        assertEquals(cobaltReferencePos, segments.get(2).start());
        assertEquals(CHR1_CENTROMERE.position() - 1, segments.get(2).end());
        assertEquals(chr1.toString(), segments.get(2).chromosome());
        assertEquals(SegmentSupport.NONE, segments.get(2).support());
    }

    @Test
    public void cobaltRatiosTumorRatio()
    {
        int sv1Start = 100_500;
        int sv2Start = 200_500;
        StructuralVariant sv1 = sv(chr1, sv1Start, INV);
        StructuralVariant sv2 = sv(chr1, sv2Start, INV);

        LinkedHashMap<Integer, Double> positionToApproximateRatio = new LinkedHashMap<>();
        positionToApproximateRatio.put(105_001, 2.0);
        positionToApproximateRatio.put(205_000, 0.5);
        positionToApproximateRatio.put(800_000, 1.0);
        List<CobaltRatio> cobaltRatios = createCobaltRatios(chr1, positionToApproximateRatio,1, CHR1_CENTROMERE.position());
        CobaltData cobaltData = new CobaltData(cobaltChromosomes);
        cobaltData.Ratios.put(chr1, cobaltRatios);

        List<ObservedRegion> segments = segmentation.createObservedRegions(of(sv1, sv2), amberPcfPerChromosomeArm(), cobaltData);
        assertEquals(8, segments.size());

        assertEquals(1, segments.get(0).start());
        assertEquals(sv1Start, segments.get(0).end());
        assertEquals(TELOMERE, segments.get(0).support());
        assertEquals(2.0, segments.get(0).observedTumorRatio(), 0.1);

        assertEquals(sv1Start + 1, segments.get(1).start());
        assertEquals(sv1Start + SV_LENGTH, segments.get(1).end());
        assertEquals(SegmentSupport.INV, segments.get(1).support());
        assertEquals(2.0, segments.get(1).observedTumorRatio(), 0.1);

        assertEquals(sv1Start + SV_LENGTH + 1, segments.get(2).start());
        assertEquals(sv2Start, segments.get(2).end());
        assertEquals(SegmentSupport.INV, segments.get(2).support());
        assertEquals(0.5, segments.get(2).observedTumorRatio(), 0.1);

        assertEquals(sv2Start + 1, segments.get(3).start());
        assertEquals(sv2Start + SV_LENGTH, segments.get(3).end());
        assertEquals(SegmentSupport.INV, segments.get(3).support());
        assertEquals(1.0, segments.get(3).observedTumorRatio(), 0.1);

        assertEquals(sv2Start + SV_LENGTH + 1, segments.get(4).start());
        assertEquals(CHR1_CENTROMERE.position() - 1, segments.get(4).end());
        assertEquals(SegmentSupport.INV, segments.get(4).support());
        assertEquals(1.0, segments.get(4).observedTumorRatio(), 0.1);

        // Now with the same structural variants but slight changes to the cobalt ratio change positions.
        // The observed tumor ratios for corresponding regions are quite different because they are the
        // medians of the values from the CobaltRatios that are in the region.
        positionToApproximateRatio = new LinkedHashMap<>();
        positionToApproximateRatio.put(105_000, 2.0);
        positionToApproximateRatio.put(205_001, 0.5);
        positionToApproximateRatio.put(800_000, 1.0);
        cobaltRatios = createCobaltRatios(chr1, positionToApproximateRatio,1, CHR1_CENTROMERE.position());
        cobaltData = new CobaltData(cobaltChromosomes);
        cobaltData.Ratios.put(chr1, cobaltRatios);

        segments = segmentation.createObservedRegions(of(sv1, sv2), amberPcfPerChromosomeArm(), cobaltData);
        assertEquals(8, segments.size());

        assertEquals(1, segments.get(0).start());
        assertEquals(sv1Start, segments.get(0).end());
        assertEquals(TELOMERE, segments.get(0).support());
        assertEquals(2.0, segments.get(0).observedTumorRatio(), 0.1);

        assertEquals(sv1Start + 1, segments.get(1).start());
        assertEquals(sv1Start + SV_LENGTH, segments.get(1).end());
        assertEquals(SegmentSupport.INV, segments.get(1).support());
        assertEquals(0.5, segments.get(1).observedTumorRatio(), 0.1);

        assertEquals(sv1Start + SV_LENGTH + 1, segments.get(2).start());
        assertEquals(sv2Start, segments.get(2).end());
        assertEquals(SegmentSupport.INV, segments.get(2).support());
        assertEquals(0.5, segments.get(2).observedTumorRatio(), 0.1);

        assertEquals(sv2Start + 1, segments.get(3).start());
        assertEquals(sv2Start + SV_LENGTH, segments.get(3).end());
        assertEquals(SegmentSupport.INV, segments.get(3).support());
        assertEquals(0.5, segments.get(3).observedTumorRatio(), 0.1);

        assertEquals(sv2Start + SV_LENGTH + 1, segments.get(4).start());
        assertEquals(CHR1_CENTROMERE.position() - 1, segments.get(4).end());
        assertEquals(SegmentSupport.INV, segments.get(4).support());
        assertEquals(1.0, segments.get(4).observedTumorRatio(), 0.1);
    }

    private List<CobaltRatio> createCobaltRatios(Chromosome chromosome, Map<Integer, Double> positionToApproximateRatio, int armStart, int armEnd)
    {
        List<CobaltRatio> cobaltRatios = new ArrayList<>();
        Iterator<Integer> limitsIterator = positionToApproximateRatio.keySet().iterator();
        int currentRegionEnd = limitsIterator.next();
        Double currentApproximateRatio = positionToApproximateRatio.get(currentRegionEnd);
        assertTrue(currentRegionEnd > armStart);
        int windowStart = armStart;
        while (windowStart < armEnd)
        {
            if (windowStart > currentRegionEnd){
                if (limitsIterator.hasNext()){
                    int nextRegionEnd = limitsIterator.next();
                    assertTrue(nextRegionEnd > currentRegionEnd);
                    currentRegionEnd = nextRegionEnd;
                    currentApproximateRatio = positionToApproximateRatio.get(currentRegionEnd);
                }
                else
                {
                    currentRegionEnd = armEnd;
                    currentApproximateRatio = null;
                }
            }
            cobaltRatios.add(cobaltRatio(chromosome, windowStart, currentApproximateRatio));
            windowStart += 1000;
        }
        return cobaltRatios;
    }

    private CobaltRatio cobaltRatio(Chromosome chromosome, int position, Double approximateRatio)
    {
        if (approximateRatio == null)
        {
            return PurpleTestUtils.cobalt(chromosome.toString(), position, -1.0);
        }
        double ratio = RandomUtils.nextDouble(0.98 * approximateRatio, 1.02 * approximateRatio);
        return PurpleTestUtils.cobalt(chromosome.toString(), position, ratio);
    }

    private StructuralVariant sv(Chromosome chromosome, int position)
    {
        return sv(chromosome, position, INV);
    }

    private StructuralVariant sv(Chromosome chromosome, int position, StructuralVariantType type)
    {
        return createStructuralVariant(chromosome.toString(), position, chromosome.toString(), position + SV_LENGTH, type).build();
    }

    private AmberData amberPcfPerChromosomeArm()
    {
        AmberData data = new AmberData(50, Gender.MALE);
        data.TumorSegments.put(chr1, new PCFPosition(PCFSource.TUMOR_BAF, chr1.toString(), 1));
        data.TumorSegments.put(chr1, new PCFPosition(PCFSource.TUMOR_BAF, chr1.toString(), CHR1_CENTROMERE.position()));
        data.TumorSegments.put(chr2, new PCFPosition(PCFSource.TUMOR_BAF, chr2.toString(), 1));
        data.TumorSegments.put(chr2, new PCFPosition(PCFSource.TUMOR_BAF, chr2.toString(), CHR2_CENTROMERE.position()));
        return data;
    }

    private AmberData emptyAmberData()
    {
        return new AmberData(50, Gender.MALE);
    }

    private CobaltData emptyCobaltData()
    {
        return new CobaltData(cobaltChromosomes);
    }
}