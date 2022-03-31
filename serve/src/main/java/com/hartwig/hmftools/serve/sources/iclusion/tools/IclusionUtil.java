package com.hartwig.hmftools.serve.sources.iclusion.tools;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.iclusion.datamodel.IclusionMutation;
import com.hartwig.hmftools.iclusion.datamodel.IclusionMutationCondition;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTrial;
import com.hartwig.hmftools.serve.extraction.ExtractionResult;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class IclusionUtil {

    private static final Logger LOGGER = LogManager.getLogger(IclusionUtil.class);

    private static final String FIELD_DELIMITER = "\t";

    private IclusionUtil() {
    }

    public static void printIclusionResult(@NotNull ExtractionResult result) {
        LOGGER.info("Analysis performed on iClusion extraction result");
        LOGGER.info(" {} actionable hotspot records generated", result.actionableHotspots().size());
        LOGGER.info(" {} actionable range records generated", result.actionableRanges().size());
        LOGGER.info(" {} actionable gene records generated", result.actionableGenes().size());
        LOGGER.info(" {} actionable fusion records generated", result.actionableFusions().size());
        LOGGER.info(" {} actionable tumor characteristics records generated", result.actionableCharacteristics().size());
    }

    public static void writeIclusionMutationTypes(@NotNull String iClusionMutationTsv, @NotNull List<IclusionTrial> trials)
            throws IOException {
        List<String> lines = Lists.newArrayList();
        String header = new StringJoiner(FIELD_DELIMITER).add("gene").add("name").add("type").toString();
        lines.add(header);

        Set<MutationEntry> mutationEntries = Sets.newHashSet();

        for (IclusionTrial trial : trials) {
            for (IclusionMutationCondition condition : trial.mutationConditions()) {
                for (IclusionMutation mutation : condition.mutations()) {
                    mutationEntries.add(new MutationEntry(mutation.gene(), mutation.name(), mutation.type().toString()));
                }
            }
        }

        for (MutationEntry entry : mutationEntries) {
            lines.add(new StringJoiner(FIELD_DELIMITER).add(entry.gene()).add(entry.name()).add(entry.type()).toString());
        }

        LOGGER.info("Writing {} unique iClusion mutations to {}", lines.size() - 1, iClusionMutationTsv);
        Files.write(new File(iClusionMutationTsv).toPath(), lines);
    }

    private static class MutationEntry {

        @Nullable
        private final String gene;
        @NotNull
        private final String name;
        @NotNull
        private final String type;

        public MutationEntry(@Nullable final String gene, @NotNull final String name, @NotNull final String type) {
            this.gene = gene;
            this.name = name;
            this.type = type;
        }

        @Nullable
        public String gene() {
            return gene;
        }

        @NotNull
        public String name() {
            return name;
        }

        @NotNull
        public String type() {
            return type;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) {
                return true;
            }
            if (o == null || getClass() != o.getClass()) {
                return false;
            }
            final MutationEntry that = (MutationEntry) o;
            return Objects.equals(gene, that.gene) && name.equals(that.name) && type.equals(that.type);
        }

        @Override
        public int hashCode() {
            return Objects.hash(gene, name, type);
        }
    }
}
