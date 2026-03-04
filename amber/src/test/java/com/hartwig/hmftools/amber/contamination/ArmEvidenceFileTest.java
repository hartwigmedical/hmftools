package com.hartwig.hmftools.amber.contamination;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.segmentation.Arm;
import com.hartwig.hmftools.common.segmentation.ChrArm;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

public class ArmEvidenceFileTest
{
    private List<CategoryEvidence<ChrArm>> Evidence;

    @Before
    public void setup()
    {
        Evidence = new ArrayList<>();
    }

    @Test
    public void writeReadTest() throws IOException
    {
        addEvidence(HumanChromosome._1, Arm.P, 12, 24);
        addEvidence(HumanChromosome._1, Arm.Q, 0, 20);
        addEvidence(HumanChromosome._X, Arm.P, 60, 120);
        addEvidence(HumanChromosome._X, Arm.Q, 100, 120);
        File tempDir = Files.createTempDirectory("amber").toFile();
        tempDir.deleteOnExit();
        File destination = new File(tempDir, "arms.tsv");
        ArmEvidenceFile.write(destination.getAbsolutePath(), Evidence);
        List<CategoryEvidence<ChrArm>> read = ArmEvidenceFile.read(destination.getAbsolutePath());
        Assert.assertEquals(Evidence.stream().sorted().toList(), read);
    }

    private void addEvidence(HumanChromosome chromosome, Arm arm, int evidencePoints, int totalPoints)
    {
        Preconditions.checkArgument(totalPoints >= evidencePoints);
        CategoryEvidence<ChrArm> e = new CategoryEvidence<>(new ChrArm(chromosome, arm));
        for(int j = 0; j < totalPoints; j++)
        {
            e.register(j < evidencePoints);
        }
        Evidence.add(e);
    }
}
