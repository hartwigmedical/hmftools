package com.hartwig.hmftools.linx.visualiser.data;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.visualiser.file.VisFusionFile;

import org.jetbrains.annotations.NotNull;

public class VisFusions
{
    @NotNull
    public static List<Fusion> fromFile(@NotNull final String fileName) throws IOException
    {
        if(!Files.exists(Paths.get(fileName)))
            return Lists.newArrayList();

        final List<VisFusionFile> visFusions = VisFusionFile.read(fileName);

        return visFusions.stream().map(x -> ImmutableFusion.builder()
                .sampleId(x.SampleId)
                .reportable(x.Reportable)
                .clusterId(x.ClusterId)
                .geneUp(x.GeneNameUp)
                .transcriptUp(x.TranscriptUp)
                .chromosomeUp(x.ChrUp)
                .positionUp(x.PosUp)
                .strandUp(x.StrandUp)
                .regionTypeUp(x.RegionTypeUp)
                .fusedExonUp(x.FusedExonUp)
                .geneDown(x.GeneNameDown)
                .transcriptDown(x.TranscriptDown)
                .chromosomeDown(x.ChrDown)
                .positionDown(x.PosDown)
                .strandDown(x.StrandDown)
                .regionTypeDown(x.RegionTypeDown)
                .fusedExonDown(x.FusedExonDown)
                .build()).collect(Collectors.toList());
    }
}
