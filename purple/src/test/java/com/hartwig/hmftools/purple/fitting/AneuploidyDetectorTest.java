package com.hartwig.hmftools.purple.fitting;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._10;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._3;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._4;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._5;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._6;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._7;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._8;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._9;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.purple.GermlineStatus.LIKELY_DIPLOID;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.chromosome.PerChromosomeData;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.purple.region.FittingRegion;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

public class AneuploidyDetectorTest
{
    private static final double JUST_ABOVE_CUTOFF = 0.571;
    private static final HumanChromosome[] FIRST_10 = { _1, _2, _3, _4, _5, _6, _7, _8, _9, _10 };
    private static final double VERY_HIGH = 0.701;
    private List<FittingRegion> Regions;
    private ListMultimap<Chromosome, AmberBAF> AmberData;

    @Before
    public void setup()
    {
        Regions = new ArrayList<>();
        AmberData = ArrayListMultimap.create();
    }

    @Test
    public void regionUsedForAneuploidyDetectionTest()
    {
        PerChromosomeData pcd = new PCD(_1, _2, _3);

        FittingRegion fr1 = fr(_1, 100, 200, GermlineStatus.DIPLOID, 10, 0.5);
        Assert.assertTrue(AneuploidyDetector.useRegionToDetectAneuploidy(pcd, fr1));

        FittingRegion fr2 = fr(_10, 100, 200, GermlineStatus.DIPLOID, 10, 0.5);
        assertFalse(AneuploidyDetector.useRegionToDetectAneuploidy(pcd, fr2));

        FittingRegion fr3 = fr(_1, 100, 200, LIKELY_DIPLOID, 10, 0.5);
        assertFalse(AneuploidyDetector.useRegionToDetectAneuploidy(pcd, fr3));

        FittingRegion fr4 = fr(_1, 100, 200, GermlineStatus.DIPLOID, -1, 0.5);
        assertFalse(AneuploidyDetector.useRegionToDetectAneuploidy(pcd, fr4));

        FittingRegion fr5 = fr(_1, 100, 200, GermlineStatus.DIPLOID, 10, -0.5);
        assertFalse(AneuploidyDetector.useRegionToDetectAneuploidy(pcd, fr5));
    }

    @Test
    public void noRegionsElevatedTest()
    {
        for(HumanChromosome chromosome : FIRST_10)
        {
            int position = 1000;
            for(int i = 0; i < 100; i++)
            {
                createAndAddRegion(chromosome, position, position + 1000, 10, 0.5);
                position += 2000;
            }
        }
        AneuploidyDetector detector = new AneuploidyDetector(Regions, AmberData, pcdFirst10Chromosomes());
        assertFalse(detector.hasAneuploidy());
    }

    @Test
    public void onePercentOfRegionsElevatedTest()
    {
        for(HumanChromosome chromosome : FIRST_10)
        {
            int position = 1000;
            for(int i = 0; i < 100; i++)
            {
                if(i == 10)
                {
                    createAndAddRegion(chromosome, position, position + 1000, 10, JUST_ABOVE_CUTOFF);
                }
                else
                {
                    createAndAddRegion(chromosome, position, position + 1000, 10, 0.5);
                }
                position += 2000;
            }
        }
        AneuploidyDetector detector = new AneuploidyDetector(Regions, AmberData, pcdFirst10Chromosomes());
        assertTrue(detector.hasAneuploidy());
    }

    @Test
    public void minimalPercentageOfBafsElevatedTest()
    {
        assertFalse(checkForElevateBafCountOutOfAbout5000(40));
        setup();
        assertTrue(checkForElevateBafCountOutOfAbout5000(41));
    }

    @Test
    public void elevatedRegionsAreTooShortToCountTest()
    {
        assertFalse(checkFor4PercentElevated(1)); // too short
        setup();
        assertTrue(checkFor4PercentElevated(2));
    }

    @Test
    public void inputRegionsAreFilteredTest()
    {
        for(HumanChromosome chromosome : FIRST_10)
        {
            int position = 1000;
            for(int i = 0; i < 100; i++)
            {
                if(i == 10)
                {
                    Regions.add(new FR(chromosome, position, position + 1000, LIKELY_DIPLOID, 10, JUST_ABOVE_CUTOFF, 0.5, 0.5));
                }
                else
                {
                    createAndAddRegion(chromosome, position, position + 1000, 10, 0.5);
                }
                position += 2000;
            }
        }
        AneuploidyDetector detector = new AneuploidyDetector(Regions, AmberData, pcdFirst10Chromosomes());
        assertFalse(detector.hasAneuploidy());
    }

    @Test
    public void singleSignificantlyElevatedRegionTest()
    {
        PerChromosomeData pcd = new PCD(_1);

        createUniformlyDiploidRegionsPlusSingleSignificantlyElevatedRegion(10, 99, VERY_HIGH, 10, 0.2);
        AneuploidyDetector detector = new AneuploidyDetector(Regions, AmberData, pcd);
        assertTrue(detector.hasAneuploidy());

        // Cutoff % for aneuploidy from "fairly high" regions is 0.008 ~ 1/125, so with more than 125 regions
        // we exercise the code that looks for a single highly aneuploidic region.
        setup();
        createUniformlyDiploidRegionsPlusSingleSignificantlyElevatedRegion(4, 124, VERY_HIGH, 4, 0.2);
        detector = new AneuploidyDetector(Regions, AmberData, pcd);
        assertTrue(detector.hasAneuploidy());

        setup();
        createUniformlyDiploidRegionsPlusSingleSignificantlyElevatedRegion(4, 126, VERY_HIGH, 4, 0.2);
        detector = new AneuploidyDetector(Regions, AmberData, pcd);
        assertFalse(detector.hasAneuploidy()); // less than 1% of vaf in the elevated region

        setup();
        createUniformlyDiploidRegionsPlusSingleSignificantlyElevatedRegion(10, 165, VERY_HIGH, 5, 0.2);
        detector = new AneuploidyDetector(Regions, AmberData, pcd);
        assertTrue(detector.hasAneuploidy()); // 5 points in the elevated region

        setup();
        double vaf = VERY_HIGH - 0.002;
        createUniformlyDiploidRegionsPlusSingleSignificantlyElevatedRegion(5, 126, vaf, 5, 0.2);
        detector = new AneuploidyDetector(Regions, AmberData, pcd);
        assertFalse(detector.hasAneuploidy()); // region vaf not quite high enough

        setup();
        createUniformlyDiploidRegionsPlusSingleSignificantlyElevatedRegion(5, 126, VERY_HIGH, 5, 0.301);
        detector = new AneuploidyDetector(Regions, AmberData, pcd);
        assertFalse(detector.hasAneuploidy()); // min vaf in elevated region too high
    }

    private void createUniformlyDiploidRegionsPlusSingleSignificantlyElevatedRegion(
            final int pointsPerRegion,
            final int numberOfDiploidRegions,
            final double aneuploidicRegionVaf,
            final int pointsInElevatedRegion,
            final double minVafInElevatedRegion)
    {
        int position = 1000;
        // Create regions with no aneuploidy
        for(int i = 0; i < numberOfDiploidRegions; i++)
        {
            createAndAddRegion(_1, position, position + 1000, pointsPerRegion, 0.5);
            position += 2000;
        }
        // Add another region, of the same length as the other, that has significantly elevated vaf
        Regions.add(frDiploid(_1, position, position + 1000, pointsInElevatedRegion, aneuploidicRegionVaf));
        int amberPosition = createAmberPointsForRegion(_1, position, pointsInElevatedRegion, aneuploidicRegionVaf);
        // Add some low points to this region
        createAmberPointsForRegion(_1, amberPosition, 1, minVafInElevatedRegion);
    }

    private boolean checkForElevateBafCountOutOfAbout5000(final int elevatedBafPointCount)
    {
        for(HumanChromosome chromosome : FIRST_10)
        {
            int position = 1000;
            for(int i = 0; i < 1000; i++)
            {
                if(i == 10)
                {
                    createAndAddRegion(chromosome, position, position + 1000, elevatedBafPointCount, JUST_ABOVE_CUTOFF);
                }
                else
                {
                    createAndAddRegion(chromosome, position, position + 1000, 5, 0.5);
                }
                position += 2000;
            }
        }
        // total baf count = 999 * 5 + 55 = 5050

        return new AneuploidyDetector(Regions, AmberData, pcdFirst10Chromosomes()).hasAneuploidy();
    }

    private boolean checkFor4PercentElevated(final int multiplier)
    {
        for(HumanChromosome chromosome : FIRST_10)
        {
            int position = 1000;
            for(int i = 0; i < 100; i++)
            {
                if(i % 5 == 0)
                {
                    createAndAddRegion(chromosome, position, position + 1000, multiplier, JUST_ABOVE_CUTOFF);
                }
                else
                {
                    createAndAddRegion(chromosome, position, position + 1000, 5 * multiplier, 0.5);
                }
                position += 2000;
            }
        }
        return new AneuploidyDetector(Regions, AmberData, pcdFirst10Chromosomes()).hasAneuploidy();
    }

    private PerChromosomeData pcdFirst10Chromosomes()
    {
        return new PCD(FIRST_10);
    }

    private FittingRegion fr(HumanChromosome chromosome, int start, int end, GermlineStatus status, int bafCount, double tumorRatio)
    {
        return new FR(chromosome, start, end, status, bafCount, 0.5, 0.5, tumorRatio);
    }

    private void createAndAddRegion(HumanChromosome chromosome, int start, int end, int bafCount, double observedBAF)
    {
        Regions.add(frDiploid(chromosome, start, end, bafCount, observedBAF));
        createAmberPointsForRegion(chromosome, start, bafCount, observedBAF);
    }

    private int createAmberPointsForRegion(final HumanChromosome chromosome, final int start, final int bafCount, final double observedBAF)
    {
        int bafPosition = start + 10;
        for(int i = 0; i < bafCount; i++)
        {
            AmberBAF baf = new AmberBAF(V38.versionedChromosome(chromosome), bafPosition, observedBAF, 10, -1.0, 0);
            AmberData.put(chromosome, baf);
            bafPosition += 10;
        }
        return bafPosition;
    }

    private FittingRegion frDiploid(HumanChromosome chromosome, int start, int end, int bafCount, double observedBAF)
    {
        return new FR(chromosome, start, end, GermlineStatus.DIPLOID, bafCount, observedBAF, 0.5, 0.5);
    }
}

record FR(HumanChromosome chr, int start, int end,
          GermlineStatus germlineStatus, int bafCount, double observedBAF,
          double observedNormalRatio, double observedTumorRatio) implements FittingRegion
{
    @Override
    public String chromosome()
    {
        return V38.versionedChromosome(chr);
    }
}

class PCD implements PerChromosomeData
{
    private final Set<String> chrNames = new HashSet<>();

    PCD(HumanChromosome... chromosomes)
    {
        for(HumanChromosome chromosome : chromosomes)
        {
            chrNames.add(V38.versionedChromosome(chromosome));
        }
    }

    @Override
    public boolean hasChromosome(final String chromosome)
    {
        return chrNames.contains(chromosome);
    }
}
