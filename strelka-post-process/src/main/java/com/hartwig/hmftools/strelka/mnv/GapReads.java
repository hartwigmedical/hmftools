package com.hartwig.hmftools.strelka.mnv;

import java.util.Comparator;
import java.util.Map;
import java.util.Optional;

import com.google.common.collect.Maps;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
abstract class GapReads {
    private static final Logger LOGGER = LogManager.getLogger(GapReads.class);

    abstract Map<Character, Integer> countPerBases();

    Character mostFrequentRead() {
        final Optional<Map.Entry<Character, Integer>> maxCountEntry =
                countPerBases().entrySet().stream().max(Comparator.comparing(Map.Entry::getValue));
        if (maxCountEntry.isPresent()) {
            return maxCountEntry.get().getKey();
        } else {
            LOGGER.warn("Couldn't find a frequent read for this position.");
            return 'N';
        }
    }

    static GapReads empty() {
        final Map<Character, Integer> reads = Maps.newHashMap();
        reads.put('G', 0);
        reads.put('A', 0);
        reads.put('T', 0);
        reads.put('C', 0);
        reads.put('N', 0);
        return ImmutableGapReads.of(reads);
    }

    static GapReads addRead(@NotNull final GapReads reads, final char read) {
        final Map<Character, Integer> updatedReads = Maps.newHashMap();
        updatedReads.putAll(reads.countPerBases());
        updatedReads.put(read, reads.countPerBases().get(read) + 1);
        return ImmutableGapReads.of(updatedReads);
    }
}
