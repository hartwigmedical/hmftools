package com.hartwig.hmftools.serve.extraction;

import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.extraction.characteristic.TumorCharacteristicExtractor;
import com.hartwig.hmftools.serve.extraction.codon.CodonExtractor;
import com.hartwig.hmftools.serve.extraction.copynumber.CopyNumberExtractor;
import com.hartwig.hmftools.serve.extraction.exon.ExonExtractor;
import com.hartwig.hmftools.serve.extraction.fusion.FusionExtractor;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelExtractor;
import com.hartwig.hmftools.serve.extraction.hotspot.HotspotExtractor;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class EventExtractor {

    @NotNull
    private final HotspotExtractor hotspotExtractor;
    @NotNull
    private final CodonExtractor codonExtractor;
    @NotNull
    private final ExonExtractor exonExtractor;
    @NotNull
    private final GeneLevelExtractor geneLevelExtractor;
    @NotNull
    private final CopyNumberExtractor copyNumberExtractor;
    @NotNull
    private final FusionExtractor fusionExtractor;
    @NotNull
    private final TumorCharacteristicExtractor tumorCharacteristicExtractor;

    public EventExtractor(@NotNull final HotspotExtractor hotspotExtractor, @NotNull final CodonExtractor codonExtractor,
            @NotNull final ExonExtractor exonExtractor, @NotNull final GeneLevelExtractor geneLevelExtractor,
            @NotNull final CopyNumberExtractor copyNumberExtractor, @NotNull final FusionExtractor fusionExtractor,
            @NotNull final TumorCharacteristicExtractor tumorCharacteristicExtractor) {
        this.hotspotExtractor = hotspotExtractor;
        this.codonExtractor = codonExtractor;
        this.exonExtractor = exonExtractor;
        this.geneLevelExtractor = geneLevelExtractor;
        this.copyNumberExtractor = copyNumberExtractor;
        this.fusionExtractor = fusionExtractor;
        this.tumorCharacteristicExtractor = tumorCharacteristicExtractor;
    }

    @NotNull
    public EventExtractorOutput extract(@NotNull String gene, @Nullable String transcriptId, @NotNull EventType type,
            @NotNull String event) {
        return ImmutableEventExtractorOutput.builder()
                .hotspots(hotspotExtractor.extract(gene, transcriptId, type, event))
                .codons(codonExtractor.extract(gene, transcriptId, type, event))
                .exons(exonExtractor.extract(gene, transcriptId, type, event))
                .geneLevelEvent(geneLevelExtractor.extract(gene, type, event))
                .knownCopyNumber(copyNumberExtractor.extract(gene, type))
                .knownFusionPair(fusionExtractor.extract(gene, type, event))
                .characteristic(tumorCharacteristicExtractor.extract(type, event))
                .build();
    }
}
