package com.hartwig.hmftools.protect.cnchromosome;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.PurpleTestUtils;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.segment.ChromosomeArm;

import org.apache.commons.compress.utils.Lists;
import org.junit.Ignore;
import org.junit.Test;

public class CnPerChromosomeFactoryTest {

    private static final String PURPLE_COPYNUMBER_TSV = Resources.getResource("cnchromosome/sample_purple.cnv.somatic.tsv").getPath();

    private static final double EPSILON = 1.0E-3;

    @Test
    public void extractCopyNumberPerChromosomeArm() throws IOException {
        CnPerChromosomeFactory.fromPurpleSomaticCopynumberTsv(PURPLE_COPYNUMBER_TSV);
    }

    @Test
    @Ignore
    public void canDetermineCnPerChromosomeArm() {
        List<PurpleCopyNumber> copyNumbers = Lists.newArrayList();
        // Chromosome 1: 1-123035434-249250621
        copyNumbers.add(PurpleTestUtils.createCopyNumber("1", 1, 123000000, 2).build());
        copyNumbers.add(PurpleTestUtils.createCopyNumber("1", 124000000, 125000000, 300).build());
        copyNumbers.add(PurpleTestUtils.createCopyNumber("1", 125000001, 250000000, 3).build());

        Map<ChromosomeArmKey, Double> cnPerChromosomeArm = CnPerChromosomeFactory.extractCnPerChromosomeArm(copyNumbers);
        assertEquals(2D, cnPerChromosomeArm.get(new ChromosomeArmKey(HumanChromosome._1, ChromosomeArm.P_ARM)), EPSILON);
        // Below is not exactly 6 but close enough.
        assertEquals(6, cnPerChromosomeArm.get(new ChromosomeArmKey(HumanChromosome._1, ChromosomeArm.Q_ARM)), EPSILON);
    }

}