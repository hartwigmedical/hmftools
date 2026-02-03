package com.hartwig.hmftools.purple.copynumber;

import static org.junit.Assert.assertEquals;
import java.io.IOException;
import java.util.List;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.ChromosomeArm;
import com.hartwig.hmftools.common.segmentation.Arm;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

public class ChromosomeArmCopyNumbersFileTest
{
    @Rule
    public TemporaryFolder temp = new TemporaryFolder();

    @Test
    public void roundTripTest() throws IOException
    {
        List<ChromosomeArmCopyNumber> expected = Lists.newArrayList();
        expected.add(new ChromosomeArmCopyNumber(HumanChromosome._1, Arm.P, 2.0, 2.1, 1.9, 2.2));
        expected.add(new ChromosomeArmCopyNumber(HumanChromosome._2, Arm.Q, 3.0, 3.1, 2.9, 3.2));

        String filename = temp.newFile("test" + ChromosomeArmCopyNumbersFile.EXTENSION).getAbsolutePath();
        ChromosomeArmCopyNumbersFile.write(filename, expected);

        List<ChromosomeArmCopyNumber> actual = ChromosomeArmCopyNumbersFile.read(filename);
        assertEquals(expected.size(), actual.size());
        for (int i = 0; i < expected.size(); i++) {
            assertEquals(expected.get(i), actual.get(i));
        }
    }
}
