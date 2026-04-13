package com.hartwig.hmftools.panelbuilder.probequality;

import java.nio.file.Files;
import java.nio.file.Paths;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;

public class Utils
{
    public static BwaMemAligner createBwaMemAligner(String bwaIndexImageFile, int threads)
    {
        if(!Files.exists(Paths.get(bwaIndexImageFile)) || bwaIndexImageFile.isEmpty())
        {
            throw new RuntimeException("Reference genome file is missing or empty");
        }

        BwaMemIndex index = new BwaMemIndex(bwaIndexImageFile);
        BwaMemAligner aligner = new BwaMemAligner(index);

        aligner.setNThreadsOption(threads);

        return aligner;
    }
}
