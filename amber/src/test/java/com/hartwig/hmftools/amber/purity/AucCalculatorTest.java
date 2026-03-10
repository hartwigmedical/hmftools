package com.hartwig.hmftools.amber.purity;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._3;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._4;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.amber.PositionEvidence;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.segmentation.Arm;
import com.hartwig.hmftools.common.segmentation.ChrArm;
import com.hartwig.hmftools.common.segmentation.ChrArmLocator;

import org.apache.commons.lang3.tuple.Pair;
import org.junit.Assert;
import org.junit.Test;

public class AucCalculatorTest
{
    private final static int CENTRE = 10_000_000;
    private final ChrArmLocator chrArmLocatorPQ = position ->
    {
        Arm arm = position.position() < CENTRE ? Arm.P : Arm.Q;
        return new ChrArm(position.chr(), arm);
    };
    Function<PositionEvidence, ChrArm> chrArmLocatorPQClassifier =
            evidence -> chrArmLocatorPQ.map(evidence.chromosome(), evidence.position());

    @Test
    public void unevenDistributionCostTest()
    {
        // ChrArm BaselinePoints EvidencePoints
        // 1P 8 2
        // 2P 9 1
        // 3P 4 1
        // 4P 10 5
        // Order by slope, from lowest to highest, calculate the area under the curve of cumulative values.
        // 2P -> 9 * 1 * 0.5
        // 1P + 3P -> 12 * 3 * 0.5 + 12 * 1 (these have the same slope so can be merged)
        // 4P -> 10 * 5 * 0.5 + 10 * 4
        // Total area is: 4.5 + 30.0 + 65.0 = 99.5
        // Area of triangle from (0, 0) to (31, 9) = 139.5

        List<PositionEvidence> baselineValues = new ArrayList<>();
        List<PositionEvidence> hitValues = new ArrayList<>();
        // 1P
        for(int i = 1; i <= 8; i++)
        {
            baselineValues.add(peDepth10(_1, i * 1000, 5));
        }
        hitValues.add(peDepth10(_1, 9000, 3));
        hitValues.add(peDepth10(_1, 10000, 3));

        // 2P
        for(int i = 1; i <= 9; i++)
        {
            baselineValues.add(peDepth10(_2, i * 1000, 4));
        }
        hitValues.add(peDepth10(_2, 10000, 3));

        // 3P
        for(int i = 1; i <= 4; i++)
        {
            baselineValues.add(peDepth10(_3, i * 1000, 9));
        }
        hitValues.add(peDepth10(_3, 5000, 3));

        // 4P
        for(int i = 1; i <= 10; i++)
        {
            baselineValues.add(peDepth10(_4, i * 1000, 2));
        }
        for(int i = 11; i <= 15; i++)
        {
            hitValues.add(peDepth10(_4, i * 1000, 3));
        }

        Function<PositionEvidence, ChrArm> chrArmVafClassifier = positionEvidence -> new ChrArm(positionEvidence.chr(), Arm.P);
        AucCalculator<PositionEvidence, ChrArm> calculator = new AucCalculator<>(chrArmVafClassifier, baselineValues);

        assertEquals(99.5 / 139.5, calculator.calculateAuc(hitValues), 0.0001);
    }

    @Test
    public void sameContaminationOnAllChromosomeArms()
    {
        List<PositionEvidence> baselinePoints = new ArrayList<>();
        List<PositionEvidence> contaminationPoints = new ArrayList<>();
        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            // P arm
            final Pair<List<PositionEvidence>, List<PositionEvidence>> pData =
                    createEvenlySpacedDataPointsWithContamination(chromosome, 1000, 50, 10, 100, 20, 100);
            baselinePoints.addAll(pData.getLeft());
            contaminationPoints.addAll(pData.getRight());
            // Q arm
            final Pair<List<PositionEvidence>, List<PositionEvidence>> qData =
                    createEvenlySpacedDataPointsWithContamination(chromosome, CENTRE + 1000, 30, 10, 100, 20, 100);
            baselinePoints.addAll(qData.getLeft());
            contaminationPoints.addAll(qData.getRight());
        }
        AucCalculator<PositionEvidence, ChrArm> calculator = new AucCalculator<>(chrArmLocatorPQClassifier, baselinePoints);
        Assert.assertEquals(1.0, calculator.calculateAuc(contaminationPoints), 0.0001);
    }

    @Test
    public void allContaminationIn1P()
    {
        List<PositionEvidence> baselinePoints = new ArrayList<>();
        List<PositionEvidence> contaminationPoints = new ArrayList<>();
        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            if(chromosome == HumanChromosome._1)
            {
                final Pair<List<PositionEvidence>, List<PositionEvidence>> p1Data =
                        createEvenlySpacedDataPointsWithContamination(chromosome, 1000, 50, 10, 100, 20, 100);
                baselinePoints.addAll(p1Data.getLeft());
                contaminationPoints.addAll(p1Data.getRight());
                final Pair<List<PositionEvidence>, List<PositionEvidence>> q1Data =
                        createEvenlySpacedDataPointsWithContamination(chromosome, CENTRE + 1000, 30, 0, 100, 20, 100);
                baselinePoints.addAll(q1Data.getLeft());
            }
            else
            {
                final Pair<List<PositionEvidence>, List<PositionEvidence>> pData =
                        createEvenlySpacedDataPointsWithContamination(chromosome, 1000, 50, 0, 100, 20, 100);
                baselinePoints.addAll(pData.getLeft());
                final Pair<List<PositionEvidence>, List<PositionEvidence>> qData =
                        createEvenlySpacedDataPointsWithContamination(chromosome, CENTRE + 1000, 30, 0, 100, 20, 100);
                baselinePoints.addAll(qData.getLeft());
            }
        }
        // chromosomes 2-Y have 5000 in P, 3000 in Q, no contamination
        // 1P has 4500 and 500 contamination, 1Q has 3000, no contamination
        // total area = (23 * 8000 + 3000 + 4500)*500*0.5, 1P area = 4500*500*0.5
        AucCalculator<PositionEvidence, ChrArm> calculator = new AucCalculator<>(chrArmLocatorPQClassifier, baselinePoints);
        Assert.assertEquals(4500.0 / 191500.0, calculator.calculateAuc(contaminationPoints), 0.0001);
    }

    @Test
    public void highContaminationInPArmsOnly()
    {
        // Same number of data points in each arm.
        // P arms have 25% of data points with contamination.
        // Q arms have 0% of data points with contamination.
        List<PositionEvidence> baselinePoints = new ArrayList<>();
        List<PositionEvidence> contaminationPoints = new ArrayList<>();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            final Pair<List<PositionEvidence>, List<PositionEvidence>> pData =
                    createEvenlySpacedDataPointsWithContamination(chromosome, 1000, 10, 20, 100, 10, 50);
            baselinePoints.addAll(pData.getLeft());
            contaminationPoints.addAll(pData.getRight());
            final Pair<List<PositionEvidence>, List<PositionEvidence>> qData =
                    createEvenlySpacedDataPointsWithContamination(chromosome, CENTRE + 1000, 10, 0, 100, 0, 40);
            baselinePoints.addAll(qData.getLeft());
            contaminationPoints.addAll(qData.getRight());
        }
        // All P arms have 800 baseline, 200 contamination.
        // All Q arms have 1000 baseline, 0 contamination.
        // Total height = 24 * 200.
        // Contamination base = 24*800
        // Total base = 24000 + contamination base.
        // factor = (24*200*24*800) / (24000 + 24*800)*24*200 = 800/(1000 + 800)
        AucCalculator<PositionEvidence, ChrArm> calculator = new AucCalculator<>(chrArmLocatorPQClassifier, baselinePoints);
        Assert.assertEquals(4.0 / 9.0, calculator.calculateAuc(contaminationPoints), 0.0001);
    }

    @Test
    public void noContamination()
    {
        List<PositionEvidence> baselinePoints = new ArrayList<>();
        List<PositionEvidence> contaminationPoints = new ArrayList<>();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            final Pair<List<PositionEvidence>, List<PositionEvidence>> pData =
                    createEvenlySpacedDataPointsWithContamination(chromosome, 1000, 50, 0, 100, 20, 100);
            baselinePoints.addAll(pData.getLeft());
            contaminationPoints.addAll(pData.getRight());

            final Pair<List<PositionEvidence>, List<PositionEvidence>> qData =
                    createEvenlySpacedDataPointsWithContamination(chromosome, CENTRE + 1000, 30, 0, 100, 20, 100);
            baselinePoints.addAll(qData.getLeft());
            contaminationPoints.addAll(qData.getRight());

        }
        // Height is 0
        AucCalculator<PositionEvidence, ChrArm> calculator = new AucCalculator<>(chrArmLocatorPQClassifier, baselinePoints);
        Assert.assertEquals(Double.NaN, calculator.calculateAuc(contaminationPoints), 0.0001);
    }

    private Pair<List<PositionEvidence>, List<PositionEvidence>> createEvenlySpacedDataPointsWithContamination(
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
        List<PositionEvidence> baselinePoints = new ArrayList<>();
        List<PositionEvidence> contaminationPoints = new ArrayList<>();
        int currentPosition = startPosition;
        for(int block = 0; block < numberOfBlocks; block++)
        {
            for(int point = 0; point < totalNumberOfPointsPerBlock; point++)
            {
                if(point < numberOfContaminationPointsPerBlock)
                {
                    int refSupport = totalReadsInContaminatedPoints - altReadsInContaminatedPoints;
                    contaminationPoints.add(pe(chromosome, currentPosition, refSupport, altReadsInContaminatedPoints, totalReadsInContaminatedPoints));
                }
                else
                {
                    int altSupport = 0;
                    baselinePoints.add(pe(chromosome, currentPosition, totalReadsInContaminatedPoints, altSupport, totalReadsInContaminatedPoints));
                }
                currentPosition += 1000;
            }
        }
        return Pair.of(baselinePoints, contaminationPoints);
    }

    private PositionEvidence peDepth10(HumanChromosome chromosome, int position, int refSupport)
    {
        final PositionEvidence positionEvidence = new PositionEvidence(V38.versionedChromosome(chromosome), position, "A", "C");
        positionEvidence.ReadDepth = 10;
        positionEvidence.RefSupport = refSupport;
        positionEvidence.AltSupport = 10 - refSupport;
        return positionEvidence;
    }

    private PositionEvidence pe(HumanChromosome chromosome, int position, int refSupport, int altSupport, int readDepth)
    {
        final PositionEvidence positionEvidence = new PositionEvidence(V38.versionedChromosome(chromosome), position, "A", "C");
        positionEvidence.ReadDepth = readDepth;
        positionEvidence.RefSupport = refSupport;
        positionEvidence.AltSupport = altSupport;
        return positionEvidence;
    }
}
