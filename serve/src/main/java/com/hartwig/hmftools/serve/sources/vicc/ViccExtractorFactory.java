package com.hartwig.hmftools.serve.sources.vicc;

import java.util.Map;

import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.serve.hotspot.ProteinResolver;
import com.hartwig.hmftools.serve.sources.vicc.extractor.CopyNumberExtractor;
import com.hartwig.hmftools.serve.sources.vicc.extractor.FusionExtractor;
import com.hartwig.hmftools.serve.sources.vicc.extractor.GeneLevelEventExtractor;
import com.hartwig.hmftools.serve.sources.vicc.extractor.GeneRangeExtractor;
import com.hartwig.hmftools.serve.sources.vicc.extractor.HotspotExtractor;
import com.hartwig.hmftools.serve.sources.vicc.extractor.SignaturesExtractor;

import org.jetbrains.annotations.NotNull;

public final class ViccExtractorFactory {
    @NotNull
    private final Map<String, HmfTranscriptRegion> transcriptPerGeneMap;

    private ViccExtractorFactory(@NotNull Map<String, HmfTranscriptRegion> transcriptPerGeneMap) {
        this.transcriptPerGeneMap = transcriptPerGeneMap;
    }

    @NotNull
    public static ViccExtractor buildViccExtractor(@NotNull ProteinResolver proteinResolver) {
        Map<String, HmfTranscriptRegion> transcriptPerGeneMap = HmfGenePanelSupplier.allGenesMap37();

        return new ViccExtractor(new HotspotExtractor(proteinResolver),
                new CopyNumberExtractor(transcriptPerGeneMap),
                new FusionExtractor(transcriptPerGeneMap),
                new GeneLevelEventExtractor(),
                new GeneRangeExtractor(transcriptPerGeneMap),
                new SignaturesExtractor());
    }
}
