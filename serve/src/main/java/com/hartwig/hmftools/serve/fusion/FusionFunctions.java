package com.hartwig.hmftools.serve.fusion;

import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class FusionFunctions {

    private FusionFunctions() {
    }

    @NotNull
    public static List<KnownFusionPair> consolidate(@NotNull List<KnownFusionPair> fusionPairs) {
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

        List<KnownFusionPair> consolidated = Lists.newArrayList();
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

    private static class FusionKey {

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
        public String geneUp() {
            return geneUp;
        }

        @Nullable
        public Integer minExonUp() {
            return minExonUp;
        }

        @Nullable
        public Integer maxExonUp() {
            return maxExonUp;
        }

        @NotNull
        public String geneDown() {
            return geneDown;
        }

        @Nullable
        public Integer minExonDown() {
            return minExonDown;
        }

        @Nullable
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
