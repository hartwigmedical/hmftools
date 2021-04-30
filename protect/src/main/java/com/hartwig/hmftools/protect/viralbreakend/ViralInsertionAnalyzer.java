package com.hartwig.hmftools.protect.viralbreakend;

import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

public final class ViralInsertionAnalyzer {

    private static final Set<String> EXCLUDED_VIRAL_INSERTIONS =
            Sets.newHashSet("Human immunodeficiency virus", "Human immunodeficiency virus 1", "Human immunodeficiency virus 2");

    private ViralInsertionAnalyzer() {
    }

    @NotNull
    public static List<Viralbreakend> analyzeViralInsertions(@NotNull List<Viralbreakend> viralbreakends) {
        Map<ViralInsertionAnalyzer.VirusKey, List<Viralbreakend>> itemsPerKey = Maps.newHashMap();
        for (Viralbreakend viralBreakend : viralbreakends) {
            ViralInsertionAnalyzer.VirusKey key = new ViralInsertionAnalyzer.VirusKey(viralBreakend.Reference());
            List<Viralbreakend> items = itemsPerKey.get(key);

            if (items == null) {
                items = Lists.newArrayList();
            }
            items.add(viralBreakend);
            itemsPerKey.put(key, items);
        }

        List<Viralbreakend> viralInsertions = Lists.newArrayList();
        for (Map.Entry<ViralInsertionAnalyzer.VirusKey, List<Viralbreakend>> entry : itemsPerKey.entrySet()) {
            List<Viralbreakend> itemsForKey = entry.getValue();

            int count = itemsForKey.size();
            assert count > 0;
            String virusName = itemsForKey.get(0).nameAssigned();
            if (!EXCLUDED_VIRAL_INSERTIONS.contains(virusName)) {
                viralInsertions.add(ImmutableViralbreakend.builder()
                        .from(itemsForKey.get(0))
                        .integrations(Integer.toString(count)) // TODO: can we use integrations of file?
                        .build());
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
