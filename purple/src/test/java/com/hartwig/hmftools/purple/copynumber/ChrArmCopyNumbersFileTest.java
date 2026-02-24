package com.hartwig.hmftools.purple.copynumber;

import static org.junit.Assert.assertEquals;
import java.io.IOException;
import java.util.List;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.ChrArmCopyNumber;
import com.hartwig.hmftools.common.purple.ChrArmCopyNumbersFile;
import com.hartwig.hmftools.common.segmentation.Arm;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

public class ChrArmCopyNumbersFileTest
{
    @Rule
    public TemporaryFolder temp = new TemporaryFolder();

    @Test
    public void roundTripTest() throws IOException
    {
        List<ChrArmCopyNumber> expected = Lists.newArrayList();
        expected.add(new ChrArmCopyNumber(HumanChromosome._1, Arm.P, 2.0, 2.1, 1.9, 2.2));
        expected.add(new ChrArmCopyNumber(HumanChromosome._2, Arm.Q, 3.0, 3.1, 2.9, 3.2));

        String filename = temp.newFile("test" + ChrArmCopyNumbersFile.EXTENSION).getAbsolutePath();
        ChrArmCopyNumbersFile.write(filename, expected);

        List<ChrArmCopyNumber> actual = ChrArmCopyNumbersFile.read(filename);
        assertEquals(expected.size(), actual.size());

        for (int i = 0; i < expected.size(); i++)
        {
            assertEquals(expected.get(i), actual.get(i));
        }
    }
}
