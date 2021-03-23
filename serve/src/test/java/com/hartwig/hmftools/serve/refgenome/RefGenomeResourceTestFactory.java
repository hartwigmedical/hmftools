package com.hartwig.hmftools.serve.refgenome;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.extraction.hotspot.ProteinResolver;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class RefGenomeResourceTestFactory {

    private RefGenomeResourceTestFactory() {
    }

    @NotNull
    public static RefGenomeResource buildTest37() {
        return ImmutableRefGenomeResource.builder()
                .fastaFile(Strings.EMPTY)
                .canonicalTranscriptPerGeneMap(HmfGenePanelSupplier.allGenesMap37())
                .proteinResolver(new TestProteinResolver())
                .build();
    }

    private static class TestProteinResolver implements ProteinResolver {

        @NotNull
        @Override
        public List<VariantHotspot> resolve(@NotNull final String gene, @Nullable final String specificTranscript,
                @NotNull final String proteinAnnotation) {
            return Lists.newArrayList(ImmutableVariantHotspotImpl.builder().chromosome("1").position(10).ref("A").alt("T").build());
        }

        @NotNull
        @Override
        public Set<String> unresolvedProteinAnnotations() {
            return Sets.newHashSet();
        }
    }
}
