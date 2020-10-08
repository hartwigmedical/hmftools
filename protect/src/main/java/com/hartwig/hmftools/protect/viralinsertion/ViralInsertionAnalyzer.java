package com.hartwig.hmftools.protect.viralinsertion;

import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.structural.linx.LinxViralInsertion;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ViralInsertionAnalyzer {

    private static final Set<String> EXCLUDED_VIRAL_INSERTIONS =
            Sets.newHashSet("Human immunodeficiency virus", "Human immunodeficiency virus 1", "Human immunodeficiency virus 2");

    private ViralInsertionAnalyzer() {
    }

    @Nullable
    public static List<ViralInsertion> analyzeViralInsertions(@NotNull List<LinxViralInsertion> linxViralInsertions) {
        Map<ViralInsertionAnalyzer.VirusKey, List<LinxViralInsertion>> itemsPerKey = Maps.newHashMap();
        for (LinxViralInsertion viralInsertion : linxViralInsertions) {
            ViralInsertionAnalyzer.VirusKey key = new ViralInsertionAnalyzer.VirusKey(viralInsertion.VirusId);
            List<LinxViralInsertion> items = itemsPerKey.get(key);

            if (items == null) {
                items = Lists.newArrayList();
            }
            items.add(viralInsertion);
            itemsPerKey.put(key, items);
        }

        List<ViralInsertion> viralInsertions = Lists.newArrayList();
        for (Map.Entry<ViralInsertionAnalyzer.VirusKey, List<LinxViralInsertion>> entry : itemsPerKey.entrySet()) {
            List<LinxViralInsertion> itemsForKey = entry.getValue();

            int count = itemsForKey.size();
            assert count > 0;
            String virusName = itemsForKey.get(0).VirusName;
            if (!EXCLUDED_VIRAL_INSERTIONS.contains(virusName)) {
                viralInsertions.add(ImmutableViralInsertion.builder().virus(virusName).viralInsertionCount(count).build());
            }
        }

        return viralInsertions;
    }

    private static class VirusKey {

        @NotNull
        private final String virusId;

        private VirusKey(@NotNull final String virusId) {
            this.virusId = virusId;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) {
                return true;
            }
            if (o == null || getClass() != o.getClass()) {
                return false;
            }
            final ViralInsertionAnalyzer.VirusKey key = (ViralInsertionAnalyzer.VirusKey) o;
            return Objects.equals(virusId, key.virusId);
        }

        @Override
        public int hashCode() {
            return Objects.hash(virusId);
        }
    }
}
