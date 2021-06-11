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
import org.junit.Ignore;
import org.junit.Test;

public class CnPerChromosomeFactoryTest {

    private static final String PURPLE_COPYNUMBER_TSV = Resources.getResource("purple/sample_purple.cnv.somatic.tsv").getPath();

    private static final double EPSILON = 1.0E-3;

    @Test
    public void extractCopyNumberPerChromosomeArm() throws IOException {
        assertNotNull(CnPerChromosomeFactory.fromPurpleSomaticCopynumberTsv(PURPLE_COPYNUMBER_TSV));
    }

    @Test
    public void canDetermineCnPerChromosomeArm() {
        List<PurpleCopyNumber> copyNumbers = Lists.newArrayList();
        // Chromosome 1: 1-123035434-249250621
        copyNumbers.add(PurpleTestUtils.createCopyNumber("1", 1, 123035434, 2).build());
        copyNumbers.add(PurpleTestUtils.createCopyNumber("1", 123035435, 124035434, 300).build());
        copyNumbers.add(PurpleTestUtils.createCopyNumber("1", 124035435, 249250621, 3).build());

        List<CnPerChromosomeArmData> cnPerChromosomeArm =
                CnPerChromosomeFactory.extractCnPerChromosomeArm(copyNumbers, RefGenomeCoordinates.COORDS_37);

        assertEquals(2, cnPerChromosomeArm.size());
        CnPerChromosomeArmData cnPerChromosomeArmData1 = cnPerChromosomeArm.get(0);
        assertEquals(HumanChromosome.fromString("1"), cnPerChromosomeArmData1.chromosome());
        assertEquals(ChromosomeArm.P_ARM, cnPerChromosomeArmData1.chromosomeArm());
        assertEquals(2D, cnPerChromosomeArmData1.copyNumber(), EPSILON);


        CnPerChromosomeArmData cnPerChromosomeArmData2 = cnPerChromosomeArm.get(1);
        assertEquals(HumanChromosome.fromString("1"), cnPerChromosomeArmData2.chromosome());
        assertEquals(ChromosomeArm.Q_ARM, cnPerChromosomeArmData2.chromosomeArm());
        assertEquals(5.35312, cnPerChromosomeArmData2.copyNumber(), EPSILON);
    }
}