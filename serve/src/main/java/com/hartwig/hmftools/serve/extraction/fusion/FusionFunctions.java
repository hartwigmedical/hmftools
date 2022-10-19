package com.hartwig.hmftools.serve.extraction.fusion;

import java.util.Map;
import java.util.Objects;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.datamodel.fusion.FusionPair;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class FusionFunctions {

    private FusionFunctions() {
    }

    @NotNull
    public static Set<KnownFusionPair> consolidate(@NotNull Iterable<KnownFusionPair> fusionPairs) {
        Map<FusionKey, Set<Knowledgebase>> sourcesPerFusion = Maps.newHashMap();
        for (KnownFusionPair fusionPair : fusionPairs) {
            FusionKey key = new FusionKey(fusionPair);
            Set<Knowledgebase> sources = sourcesPerFusion.get(key);
            if (sources == null) {
                sources = Sets.newHashSet();
            }
            sources.addAll(fusionPair.sources());
            sourcesPerFusion.put(key, sources);
        }

        Set<KnownFusionPair> consolidated = Sets.newHashSet();
        for (Map.Entry<FusionKey, Set<Knowledgebase>> entry : sourcesPerFusion.entrySet()) {
            consolidated.add(ImmutableKnownFusionPair.builder()
                    .geneUp(entry.getKey().geneUp())
                    .minExonUp(entry.getKey().minExonUp())
                    .maxExonUp(entry.getKey().maxExonUp())
                    .geneDown(entry.getKey().geneDown())
                    .minExonDown(entry.getKey().minExonDown())
                    .maxExonDown(entry.getKey().maxExonDown())
                    .sources(entry.getValue())
                    .build());
        }
        return consolidated;
    }

    private static class FusionKey implements FusionPair {

        @NotNull
        private final String geneUp;
        @Nullable
        private final Integer minExonUp;
        @Nullable
        private final Integer maxExonUp;
        @NotNull
        private final String geneDown;
        @Nullable
        private final Integer minExonDown;
        @Nullable
        private final Integer maxExonDown;

        public FusionKey(@NotNull KnownFusionPair fusionPair) {
            this.geneUp = fusionPair.geneUp();
            this.minExonUp = fusionPair.minExonUp();
            this.maxExonUp = fusionPair.maxExonUp();
            this.geneDown = fusionPair.geneDown();
            this.minExonDown = fusionPair.minExonDown();
            this.maxExonDown = fusionPair.maxExonDown();
        }

        @NotNull
        @Override
        public String geneUp() {
            return geneUp;
        }

        @Nullable
        @Override
        public Integer minExonUp() {
            return minExonUp;
        }

        @Nullable
        @Override
        public Integer maxExonUp() {
            return maxExonUp;
        }

        @NotNull
        @Override
        public String geneDown() {
            return geneDown;
        }

        @Nullable
        @Override
        public Integer minExonDown() {
            return minExonDown;
        }

        @Nullable
        @Override
        public Integer maxExonDown() {
            return maxExonDown;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) {
                return true;
            }
            if (o == null || getClass() != o.getClass()) {
                return false;
            }
            final FusionKey fusionKey = (FusionKey) o;
            return geneUp.equals(fusionKey.geneUp) && Objects.equals(minExonUp, fusionKey.minExonUp) && Objects.equals(maxExonUp,
                    fusionKey.maxExonUp) && geneDown.equals(fusionKey.geneDown) && Objects.equals(minExonDown, fusionKey.minExonDown)
                    && Objects.equals(maxExonDown, fusionKey.maxExonDown);
        }

        @Override
        public int hashCode() {
            return Objects.hash(geneUp, minExonUp, maxExonUp, geneDown, minExonDown, maxExonDown);
        }
    }
}
