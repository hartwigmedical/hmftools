package com.hartwig.hmftools.amber.contamination;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._3;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._4;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import com.hartwig.hmftools.common.amber.BaseDepthData;
import com.hartwig.hmftools.common.amber.ImmutableBaseDepthData;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.segmentation.Arm;
import com.hartwig.hmftools.common.segmentation.ChrArm;
import com.hartwig.hmftools.common.segmentation.ChrArmLocator;

import org.junit.Assert;
import org.junit.Test;

public class ContaminationConfirmationTest
{
    private ChrArmLocator chrArmLocator = cobaltRatio -> new ChrArm(cobaltRatio.chr(), Arm.P);
    private ContaminationPredicate classifier = contamination -> contamination.Tumor.altSupport() % 10 == 3;

    @Test
    public void giniTest()
    {
        // ChrArm TotalPoints EvidencePoints
        // 1P 10 2
        // 2P 10 1
        // 3P 5 1
        // 4P 15 5
        // Order by slope, from lowest to highest, calculate area under the curve of cumulative values.
        // 2P -> 10 * 1 * 0.5
        // 1P + 3P -> 15 * 3 * 0.5 + 15 * 1 (these have the same slope so can be merged)
        // 4P -> 15 * 5 * 0.5 + 15 * 4
        // Total area is: 5 + 37.5 + 97.5 = 140
        // Area of triangle from (0, 0) to (40, 9) = 180
        // Gini = (180 - 140) / 180 = 0.2222
        ContaminationConfirmation contaminationConfirmation = new ContaminationConfirmation(classifier, chrArmLocator);

        // 1P
        for(int i = 1; i <= 8; i++)
        {
            contaminationConfirmation.offer(tc(_1, i * 1000, 10, 5));
        }
        contaminationConfirmation.offer(tc(_1, 9000, 10, 3));
        contaminationConfirmation.offer(tc(_1, 10000, 10, 3));

        // 2P
        for(int i = 1; i <= 9; i++)
        {
            contaminationConfirmation.offer(tc(_2, i * 1000, 10, 4));
        }
        contaminationConfirmation.offer(tc(_2, 10000, 10, 3));

        // 3P
        for(int i = 1; i <= 4; i++)
        {
            contaminationConfirmation.offer(tc(_3, i * 1000, 10, 9));
        }
        contaminationConfirmation.offer(tc(_3, 5000, 10, 3));

        // 4P
        for(int i = 1; i <= 10; i++)
        {
            contaminationConfirmation.offer(tc(_4, i * 1000, 10, 2));
        }
        for(int i = 11; i <= 15; i++)
        {
            contaminationConfirmation.offer(tc(_4, i * 1000, 10, 3));
        }

        Assert.assertEquals(0.2222, contaminationConfirmation.gini(), 0.0001);
    }

    private TumorContamination tc(HumanChromosome chromosome, int position, int readDepth, int refSupport)
    {
        BaseDepthData bdd = ImmutableBaseDepthData.builder()
                .ref(BaseDepthData.Base.A)
                .alt(BaseDepthData.Base.C)
                .readDepth(readDepth)
                .refSupport(refSupport)
                .altSupport(readDepth - refSupport)
                .indelCount(0)
                .build();
        return new TumorContamination(V38.versionedChromosome(chromosome), position, null, bdd);
    }
}
