package com.hartwig.hmftools.common.hla;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class HlaFiles
{
    @NotNull
    public static List<HlaTypeDetails> typeDetails(final String lilacFile) throws IOException
    {
        // never used in production, may revisit
        return Lists.newArrayList();
        //return Files.readAllLines(new File(lilacFile).toPath()).stream().skip(1).map(HlaFiles::fromString).collect(Collectors.toList());
    }

    @NotNull
    private static HlaTypeDetails fromString(String line)
    {
        String[] split = line.split("\t");
        return ImmutableHlaTypeDetails.builder()
                .type(split[0])
                .referenceUniqueCoverage(Integer.parseInt(split[2]))
                .referenceSharedCoverage(Integer.parseInt(split[3]))
                .referenceWildcardCoverage(Integer.parseInt(split[4]))
                .tumorUniqueCoverage(Integer.parseInt(split[6]))
                .tumorSharedCoverage(Integer.parseInt(split[7]))
                .tumorWildcardCoverage(Integer.parseInt(split[8]))
                .tumorCopyNumber(Double.parseDouble(split[9]))
                .somaticMissense(Double.parseDouble(split[10]))
                .somaticNonsenseOrFrameshift(Double.parseDouble(split[11]))
                .somaticSplice(Double.parseDouble(split[12]))
                .somaticSynonymous(Double.parseDouble(split[13]))
                .somaticInframeIndel(Double.parseDouble(split[14]))
                .build();
    }

    @NotNull
    public static HlaType type(final String lilacFile, final String lilacQc) throws IOException
    {
        String[] qcLines = Files.readAllLines(new File(lilacQc).toPath()).get(1).split("\t");
        String qcStatus = qcLines[0];
        int qcVariants = Integer.parseInt(qcLines[18]);

        List<String> alleles =
                Files.readAllLines(new File(lilacFile).toPath()).stream().skip(1).map(x -> x.split("\t")[0]).collect(Collectors.toList());
        if(alleles.size() != 6)
        {
            throw new IllegalStateException();
        }

        return ImmutableHlaType.builder()
                .typeA1(alleles.get(0))
                .typeA2(alleles.get(1))
                .typeB1(alleles.get(2))
                .typeB2(alleles.get(3))
                .typeC1(alleles.get(4))
                .typeC2(alleles.get(5))
                .somaticVariants(qcVariants)
                .status(qcStatus)
                .build();
    }
}
