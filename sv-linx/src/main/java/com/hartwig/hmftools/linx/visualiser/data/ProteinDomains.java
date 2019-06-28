package com.hartwig.hmftools.linx.visualiser.data;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.linx.visualiser.file.VisProteinDomainFile;

import org.jetbrains.annotations.NotNull;

public class ProteinDomains
{

    @NotNull
    public static List<ProteinDomain> readProteinDomains(@NotNull final String fileName) throws IOException
    {
        return VisProteinDomainFile.read(fileName).stream().map(ProteinDomains::fromFile).collect(Collectors.toList());
    }

    @NotNull
    private static ProteinDomain fromFile(@NotNull final VisProteinDomainFile file)
    {
        return ImmutableProteinDomain.builder()
                .sampleId(file.SampleId)
                .clusterId(file.ClusterId)
                .chromosome(file.Chromosome)
                .start(file.Start)
                .end(file.End)
                .name(file.Info)
                .build();

    }

}
