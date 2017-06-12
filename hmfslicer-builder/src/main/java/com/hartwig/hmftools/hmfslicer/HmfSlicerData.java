package com.hartwig.hmftools.hmfslicer;

import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public class HmfSlicerData {
    @NotNull
    private final GenomeRegion genomeRegion;
    @NotNull
    private final String transcriptID;
    private final int transcriptVersion;
    @NotNull
    private final String gene;
    @NotNull
    private final String chromosomeBand;
    @NotNull
    private final String entrezId;

    public HmfSlicerData(@NotNull final GenomeRegion genomeRegion, @NotNull final String transcriptID,
            final int transcriptVersion, @NotNull final String gene, @NotNull final String chromosomeBand,
            @NotNull final String entrezId) {
        this.genomeRegion = genomeRegion;
        this.transcriptID = transcriptID;
        this.transcriptVersion = transcriptVersion;
        this.gene = gene;
        this.chromosomeBand = chromosomeBand;
        this.entrezId = entrezId;
    }
}
