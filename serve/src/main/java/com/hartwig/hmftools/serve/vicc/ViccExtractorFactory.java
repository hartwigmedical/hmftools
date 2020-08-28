package com.hartwig.hmftools.serve.vicc;

import java.util.Map;

import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.serve.hotspot.HotspotGenerator;
import com.hartwig.hmftools.serve.vicc.copynumber.CopyNumberExtractor;
import com.hartwig.hmftools.serve.vicc.fusion.FusionExtractor;
import com.hartwig.hmftools.serve.vicc.hotspot.HotspotExtractor;
import com.hartwig.hmftools.serve.vicc.range.GeneLevelEventExtractor;
import com.hartwig.hmftools.serve.vicc.range.GeneRangeExtractor;
import com.hartwig.hmftools.serve.vicc.signatures.SignaturesExtractor;

import org.jetbrains.annotations.NotNull;

public final class ViccExtractorFactory {

    private ViccExtractorFactory() {
    }

    @NotNull
    public static ViccExtractor buildViccExtractor(@NotNull HotspotGenerator hotspotGenerator) {
        Map<String, HmfTranscriptRegion> transcriptPerGeneMap = HmfGenePanelSupplier.allGenesMap37();

        return new ViccExtractor(new HotspotExtractor(hotspotGenerator),
                new CopyNumberExtractor(),
                new FusionExtractor(),
                new GeneLevelEventExtractor(),
                new GeneRangeExtractor(transcriptPerGeneMap),
                new SignaturesExtractor());
    }
}
