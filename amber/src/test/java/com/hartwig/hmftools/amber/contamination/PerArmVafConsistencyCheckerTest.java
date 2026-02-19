package com.hartwig.hmftools.amber.contamination;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._3;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._4;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import java.util.ArrayList;
import java.util.List;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.amber.BaseDepthData;
import com.hartwig.hmftools.common.amber.ImmutableBaseDepthData;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.segmentation.Arm;
import com.hartwig.hmftools.common.segmentation.ChrArm;
import com.hartwig.hmftools.common.segmentation.ChrArmLocator;

import org.junit.Assert;
import org.junit.Test;

public class PerArmVafConsistencyCheckerTest
{
    private final static int CENTRE = 10_000_000;
    private final ChrArmLocator chrArmLocatorP = cobaltRatio -> new ChrArm(cobaltRatio.chr(), Arm.P);
    private final ChrArmLocator chrArmLocatorPQ = new ChrArmLocator()
    {
        @Override
        public ChrArm map(final GenomePosition position)
        {
            Arm arm = position.position() < CENTRE ? Arm.P : Arm.Q;
            return new ChrArm(position.chr(), arm);
        }
    };
    private final VafPredicate classifier = contamination -> contamination.Tumor.refSupport() % 10 == 3;

    @Test
    public void giniTest()
    {
        // ChrArm TotalPoints EvidencePoints
        // 1P 10 2
        // 2P 10 1
        // 3P 5 1
        // 4P 15 5
        // Order by slope, from lowest to highest, calculate the area under the curve of cumulative values.
        // 2P -> 10 * 1 * 0.5
        // 1P + 3P -> 15 * 3 * 0.5 + 15 * 1 (these have the same slope so can be merged)
        // 4P -> 15 * 5 * 0.5 + 15 * 4
        // Total area is: 5 + 37.5 + 97.5 = 140
        // Area of triangle from (0, 0) to (40, 9) = 180
        // Gini = 1.0 - (140 / 180) = 1.0 - 0.7778 = 0.2222
        PerArmVafConsistencyChecker perArmVafConsistencyChecker = new PerArmVafConsistencyChecker(classifier, chrArmLocatorP);

        // 1P
        for(int i = 1; i <= 8; i++)
        {
            perArmVafConsistencyChecker.offer(tc(_1, i * 1000, 5));
        }
        perArmVafConsistencyChecker.offer(tc(_1, 9000, 3));
        perArmVafConsistencyChecker.offer(tc(_1, 10000, 3));

        // 2P
        for(int i = 1; i <= 9; i++)
        {
            perArmVafConsistencyChecker.offer(tc(_2, i * 1000, 4));
        }
        perArmVafConsistencyChecker.offer(tc(_2, 10000, 3));

        // 3P
        for(int i = 1; i <= 4; i++)
        {
            perArmVafConsistencyChecker.offer(tc(_3, i * 1000, 9));
        }
        perArmVafConsistencyChecker.offer(tc(_3, 5000, 3));

        // 4P
        for(int i = 1; i <= 10; i++)
        {
            perArmVafConsistencyChecker.offer(tc(_4, i * 1000, 2));
        }
        for(int i = 11; i <= 15; i++)
        {
            perArmVafConsistencyChecker.offer(tc(_4, i * 1000, 3));
        }

        Assert.assertEquals(0.2222, perArmVafConsistencyChecker.gini().gini(), 0.0001);
    }

    @Test
    public void sameContaminationOnAllChromosomeArms()
    {
        List<TumorContamination> data = new ArrayList<>();
        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            // P arm
            data.addAll(createEvenlySpacedDataPointsWithContamination(chromosome, 1000, 50, 10, 100, 20, 100));
            // Q arm
            data.addAll(createEvenlySpacedDataPointsWithContamination(chromosome, CENTRE + 1000, 30, 10, 100, 20, 100));
        }
        double factor = PerArmVafConsistencyChecker.calculateConfirmationFactor(chrArmLocatorPQ, 0.2, data).gini();
        Assert.assertEquals(0.0, factor, 0.0001);
    }

    @Test
    public void allContaminationIn1P()
    {
        List<TumorContamination> data = new ArrayList<>();
        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            if(chromosome == HumanChromosome._1)
            {
                data.addAll(createEvenlySpacedDataPointsWithContamination(chromosome, 1000, 50, 10, 100, 20, 100));
                data.addAll(createEvenlySpacedDataPointsWithContamination(chromosome, CENTRE + 1000, 30, 0, 100, 20, 100));
            }
            data.addAll(createEvenlySpacedDataPointsWithContamination(chromosome, 1000, 50, 0, 100, 20, 100));
            data.addAll(createEvenlySpacedDataPointsWithContamination(chromosome, CENTRE + 1000, 30, 0, 100, 20, 100));
        }
        double factor = PerArmVafConsistencyChecker.calculateConfirmationFactor(chrArmLocatorPQ, 0.2, data).gini();
        Assert.assertEquals(0.0, factor, 0.0001);
    }

    @Test
    public void highContaminationInPArmsOnly()
    {
        // Same number of data points in each arm.
        // P arms have 20% of data points with contamination.
        // Q arms have 0% of data points with contamination.
        List<TumorContamination> data = new ArrayList<>();
        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            data.addAll(createEvenlySpacedDataPointsWithContamination(chromosome, 1000, 100, 20, 100, 10, 50));
            data.addAll(createEvenlySpacedDataPointsWithContamination(chromosome, CENTRE + 1000, 100, 0, 100, 0, 40));
        }
        double factor = PerArmVafConsistencyChecker.calculateConfirmationFactor(chrArmLocatorPQ, 0.2, data).gini();
        Assert.assertEquals(0.0, factor, 0.0001);
    }

    @Test
    public void noContamination()
    {
        List<TumorContamination> data = new ArrayList<>();
        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            data.addAll(createEvenlySpacedDataPointsWithContamination(chromosome, 1000, 50, 0, 100, 20, 100));
            data.addAll(createEvenlySpacedDataPointsWithContamination(chromosome, CENTRE + 1000, 30, 0, 100, 20, 100));
        }
        double factor = PerArmVafConsistencyChecker.calculateConfirmationFactor(chrArmLocatorPQ, 0.2, data).gini();
        Assert.assertEquals(0.0, factor, 0.0001);
    }

    @Test
    public void noContaminationIn5Arms()
    {
        // what about when not all arms have any data at all?
        // how many arms do we need?
    }

    private List<TumorContamination> createEvenlySpacedDataPointsWithContamination(
            HumanChromosome chromosome,
            int startPosition,
            int numberOfBlocks,
            int numberOfContaminationPointsPerBlock,
            int totalNumberOfPointsPerBlock,
            int altReadsInContaminatedPoints,
            int totalReadsInContaminatedPoints
    )
    {
        Preconditions.checkArgument(altReadsInContaminatedPoints <= totalReadsInContaminatedPoints);
        List<TumorContamination> result = new ArrayList<>();
        int currentPosition = startPosition;
        for(int block = 0; block < numberOfBlocks; block++)
        {
            for(int point = 0; point < totalNumberOfPointsPerBlock; point++)
            {
                if(point < numberOfContaminationPointsPerBlock)
                {
                    int refSupport = totalReadsInContaminatedPoints - altReadsInContaminatedPoints;
                    TumorContamination tc =
                            tc(chromosome, currentPosition, refSupport, altReadsInContaminatedPoints, totalReadsInContaminatedPoints);
                    result.add(tc);
                }
                else
                {
                    int altSupport = 0;
                    TumorContamination tc =
                            tc(chromosome, currentPosition, totalReadsInContaminatedPoints, altSupport, totalReadsInContaminatedPoints);
                    result.add(tc);
                }
                currentPosition += 1000;
            }
        }
        return result;
    }

    private TumorContamination tc(HumanChromosome chromosome, int position, int refSupport)
    {
        BaseDepthData bdd = ImmutableBaseDepthData.builder()
                .ref(BaseDepthData.Base.A)
                .alt(BaseDepthData.Base.C)
                .readDepth(10)
                .refSupport(refSupport)
                .altSupport(10 - refSupport)
                .indelCount(0)
                .build();
        return new TumorContamination(V38.versionedChromosome(chromosome), position, null, bdd);
    }

    private TumorContamination tc(HumanChromosome chromosome, int position, int refSupport, int altSupport, int readDepth)
    {
        BaseDepthData bdd = ImmutableBaseDepthData.builder()
                .ref(BaseDepthData.Base.A)
                .alt(BaseDepthData.Base.C)
                .readDepth(readDepth)
                .refSupport(refSupport)
                .altSupport(altSupport)
                .indelCount(0)
                .build();
        return new TumorContamination(V38.versionedChromosome(chromosome), position, null, bdd);
    }
}
