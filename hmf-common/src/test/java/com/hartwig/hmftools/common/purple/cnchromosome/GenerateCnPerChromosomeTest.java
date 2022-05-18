package com.hartwig.hmftools.common.purple.cnchromosome;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.purple.PurpleTestUtils;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.segment.ChromosomeArm;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GenerateCnPerChromosomeTest
{
    private static final double EPSILON = 1.0E-3;

    private static final String PURPLE_COPYNUMBER_TSV = Resources.getResource("purple/sample.purple.cnv.somatic.tsv").getPath();

    @Test
    public void canExtractCopyNumberPerChromosomeArmFromFile() throws IOException
    {
        List<CnPerChromosomeArmData> cnPerChromosomeArmData =
                GenerateCnPerChromosome.fromPurpleSomaticCopynumberTsv(PURPLE_COPYNUMBER_TSV, RefGenomeCoordinates.COORDS_37);

        assertNotNull(cnPerChromosomeArmData);
    }

    @Test
    public void canDetermineCnPerChromosomeArm()
    {
        List<PurpleCopyNumber> copyNumbers = Lists.newArrayList();
        // Chromosome 1: 1-123035434-249250621
        copyNumbers.add(PurpleTestUtils.createCopyNumber("1", 1, 123035434, 2).build());
        copyNumbers.add(PurpleTestUtils.createCopyNumber("1", 123035435, 124035434, 300).build());
        copyNumbers.add(PurpleTestUtils.createCopyNumber("1", 124035435, 249250621, 3).build());

        List<CnPerChromosomeArmData> cnPerChromosomeArm =
                GenerateCnPerChromosome.extractCnPerChromosomeArm(copyNumbers, RefGenomeCoordinates.COORDS_37);

        assertEquals(2, cnPerChromosomeArm.size());
        CnPerChromosomeArmData cnPerChromosomeArmData1 =
                findByChromosomeAndArm(cnPerChromosomeArm, HumanChromosome._1, ChromosomeArm.P_ARM);
        assertEquals(2D, cnPerChromosomeArmData1.copyNumber(), EPSILON);

        CnPerChromosomeArmData cnPerChromosomeArmData2 =
                findByChromosomeAndArm(cnPerChromosomeArm, HumanChromosome._1, ChromosomeArm.Q_ARM);
        assertEquals(5.35312, cnPerChromosomeArmData2.copyNumber(), EPSILON);
    }

    @NotNull
    private static CnPerChromosomeArmData findByChromosomeAndArm(@NotNull List<CnPerChromosomeArmData> dataList,
            @NotNull HumanChromosome chromosome, @NotNull ChromosomeArm arm)
    {
        for(CnPerChromosomeArmData data : dataList)
        {
            if(data.chromosome() == chromosome && data.chromosomeArm() == arm)
            {
                return data;
            }
        }

        throw new IllegalStateException("Could not find data with chromosome " + chromosome + " and arm " + arm);
    }
}