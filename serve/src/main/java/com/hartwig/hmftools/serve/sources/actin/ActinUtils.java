package com.hartwig.hmftools.serve.sources.actin;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinRule;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ActinUtils {

    private static final Logger LOGGER = LogManager.getLogger(ActinUtils.class);
    private static final String FIELD_DELIMITER = "\t";

    private ActinUtils() {
    }

    public static void writeActinEventTypes(@NotNull String actinMutationTsv, @NotNull List<ActinEntry> entries) throws IOException {
        List<String> lines = Lists.newArrayList();
        String header = new StringJoiner(FIELD_DELIMITER).add("study").add("rule").add("gene").add("mutation").add("type").toString();
        lines.add(header);

        Set<ActinEntryTrial> mutationEntries = Sets.newHashSet();

        for (ActinEntry entry : entries) {
            for (EventType type : entry.type()) {
                mutationEntries.add(new ActinEntryTrial(entry.trial(), entry.rule(), entry.gene(), entry.mutation(), type));
            }
        }

        for (ActinEntryTrial entry : mutationEntries) {
            lines.add(new StringJoiner(FIELD_DELIMITER).add(entry.trial)
                    .add(entry.actinRule.toString())
                    .add(entry.gene)
                    .add(entry.mutation)
                    .add(entry.type.toString())
                    .toString());
        }

        LOGGER.info("Writing {} unique ACTIN mutations to {}", lines.size() - 1, actinMutationTsv);
        Files.write(new File(actinMutationTsv).toPath(), lines);
    }

    private static class ActinEntryTrial {

        @Nullable
        private final String trial;
        @NotNull
        private final ActinRule actinRule;
        @NotNull
        private final String gene;
        @Nullable
        private final String mutation;
        @NotNull
        private final EventType type;

        public ActinEntryTrial(@Nullable final String trial, @NotNull final ActinRule actinRule, @NotNull final String gene,
                @Nullable final String mutation, @NotNull final EventType type) {
            this.trial = trial;
            this.actinRule = actinRule;
            this.gene = gene;
            this.mutation = mutation;
            this.type = type;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) {
                return true;
            }
            if (o == null || getClass() != o.getClass()) {
                return false;
            }
            final ActinEntryTrial that = (ActinEntryTrial) o;
            return Objects.equals(trial, that.trial) && actinRule == that.actinRule && gene.equals(that.gene) && Objects.equals(mutation,
                    that.mutation) && type == that.type;
        }

        @Override
        public int hashCode() {
            return Objects.hash(trial, actinRule, gene, mutation, type);
        }
    }
}
