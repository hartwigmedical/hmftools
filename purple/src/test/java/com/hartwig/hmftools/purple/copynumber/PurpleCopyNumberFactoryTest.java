package com.hartwig.hmftools.purple.copynumber;

import static com.hartwig.hmftools.purple.copynumber.PurpleCopyNumberFactory.calculateDeletedDepthWindows;

import static org.junit.Assert.assertEquals;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleTestUtils;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.purple.FittingTestUtils;
import com.hartwig.hmftools.purple.fitting.PurityAdjuster;
import com.hartwig.hmftools.purple.region.ObservedRegion;

import org.junit.Before;
import org.junit.Test;

public class PurpleCopyNumberFactoryTest
{
    private final Chromosome chr1 = HumanChromosome._1;
    private final Chromosome chr2 = HumanChromosome._2;
    private CobaltChromosomes cobaltChromosomes;
    PurityAdjuster purityAdjuster;
    PurpleCopyNumberFactory factory;
    List<PurpleCopyNumber> copyNumbers;
    List<ObservedRegion> observedRegions;
    List<PurpleCopyNumber> result;
    List<StructuralVariant> structuralVariants;

    @Before
    public void setup()
    {
        Set<MedianRatio> ratios = new HashSet<>();
        ratios.add(new MedianRatio(chr1.toString(), 0.5, 50));
        ratios.add(new MedianRatio(chr2.toString(), 0.5, 50));
        cobaltChromosomes = new CobaltChromosomes(ratios);
        purityAdjuster = new PurityAdjuster(1.0, 1.0, cobaltChromosomes);
        factory = new PurpleCopyNumberFactory(1, 1, 100, 1.0, purityAdjuster, cobaltChromosomes);

        copyNumbers = Lists.newArrayList();
        observedRegions = Lists.newArrayList();
        structuralVariants = Lists.newArrayList();
    }

    @Test
    public void mergedRegionCoversGapBetweenConstituentRegions()
    {
        createObservedRegion(10_000, 20_000, 2, GermlineStatus.DIPLOID);
        createObservedRegion(30_000, 40_000, 20, GermlineStatus.UNKNOWN);
        doBuild();

        assertEquals(1, result.size());
        PurpleCopyNumber cn = result.get(0);
        assertEquals(10_000, cn.start());
        assertEquals(40_000, cn.end());
    }

    @Test
    public void copyNumberIsExtendedFromLeft()
    {
        createObservedRegion(10_000, 20_000, 2, GermlineStatus.DIPLOID);
        createObservedRegion(30_000, 40_000, 20, GermlineStatus.UNKNOWN);
        doBuild();

        PurpleCopyNumber cn = result.get(0);
        assertEquals(2.0, cn.averageTumorCopyNumber(), 0.0001);
    }

    @Test
    public void noMergeAcrossCentromere()
    {
        createObservedRegion(10_000, 20_000, 2, SegmentSupport.TELOMERE);
        createObservedRegion(30_000, 40_000, 20, SegmentSupport.CENTROMERE);
        doBuild();

        assertEquals(2, result.size());
        PurpleCopyNumber cn = result.get(0);
        assertEquals(SegmentSupport.TELOMERE, cn.segmentStartSupport());
        assertEquals(SegmentSupport.CENTROMERE, cn.segmentEndSupport());

        cn = result.get(1);
        assertEquals(SegmentSupport.CENTROMERE, cn.segmentStartSupport());
        assertEquals(SegmentSupport.TELOMERE, cn.segmentEndSupport());
    }

    @Test
    public void testDeletedDepthWindowsCalc()
    {
        createCopyNumber("1", 1, 10000, 2, 10);
        createCopyNumber("1", 10001, 20000, 0.4, 10);
        createCopyNumber("1", 20001, 30000, 2, 10);
        createCopyNumber("1", 30001, 40000, 0.1, 10);
        createCopyNumber("1", 40001, 50000, 2, 60);

        createCopyNumber("9", 8000000, 10000000, 0.4, 40);
        createCopyNumber("9", 10000001, 11000000, 2, 40);
        createCopyNumber("9", 11000001, 13000000, 0.4, 20);

        // total = 200

        double deletedPercent = calculateDeletedDepthWindows(copyNumbers);
        assertEquals(0.25, deletedPercent, 0.01);
    }

    private void doBuild()
    {
        factory.buildCopyNumbers(observedRegions, structuralVariants);
        result = factory.copyNumbers();
    }

    private void createCopyNumber(String chromosome, int start, int end, double copyNumber, int depthWindowCount)
    {
        copyNumbers.add(PurpleTestUtils.createCopyNumber(chromosome, start, end, copyNumber)
                .depthWindowCount(depthWindowCount)
                .build());
    }

    private void createObservedRegion(int start, int end, double copyNumber, SegmentSupport support)
    {
        ObservedRegion or =
                FittingTestUtils.createObservedRegion(chr1.toString(), start, end, support, 0.5, 0.5, GermlineStatus.UNKNOWN, copyNumber);
        observedRegions.add(or);
    }

    private void createObservedRegion(int start, int end, double copyNumber, GermlineStatus status)
    {
        ObservedRegion or =
                FittingTestUtils.createObservedRegion(chr1.toString(), start, end, SegmentSupport.NONE, 0.5, 0.5, status, copyNumber);
        observedRegions.add(or);
    }

    private void createObservedRegion(int start, int end, double copyNumber, SegmentSupport segmentSupport, GermlineStatus status)
    {
        ObservedRegion or =
                FittingTestUtils.createObservedRegion(chr1.toString(), start, end, segmentSupport, 0.5, 0.5, status, copyNumber);
        observedRegions.add(or);
    }

    private void createStructuralVariant(int start, int end)
    {
        structuralVariants.add(PurpleTestUtils.createStructuralVariant(chr1.toString(), start, chr1.toString(), end, StructuralVariantType.INV)
                .build());
    }
}
