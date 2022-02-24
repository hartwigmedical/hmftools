package com.hartwig.hmftools.serve.actionability;

import java.util.Objects;
import java.util.Set;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.serve.sources.Sources;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ActionabilityTestUtil {

    public static final String TEST_SERVE_OUTPUT_DIR = Resources.getResource("serve_output").getPath();

    private ActionabilityTestUtil() {
    }

    @NotNull
    public static ActionableEvent create(@NotNull Sources source, @NotNull String treatment,
            @NotNull String cancerType, @NotNull String doid, @NotNull String tumorLocationBlacklisting,
            @NotNull EvidenceLevel level, @NotNull EvidenceDirection direction, @Nullable Set<String> sourceUrls,
            @NotNull Set<String> evidenceUrls) {
        return new ActionableEventImpl(
                source,
                treatment,
                cancerType,
                doid,
                tumorLocationBlacklisting,
                level,
                direction,
                sourceUrls,
                evidenceUrls);
    }

    private static class ActionableEventImpl implements ActionableEvent {

        @NotNull
        private final Sources source;
        @NotNull
        private final String treatment;
        @NotNull
        private final String cancerType;
        @NotNull
        private final String doid;
        @NotNull
        private final String tumorLocationBlacklisting;
        @NotNull
        private final EvidenceLevel level;
        @NotNull
        private final EvidenceDirection direction;
        @Nullable
        private final Set<String> sourceUrls;
        @NotNull
        private final Set<String> evidenceUrls;

        public ActionableEventImpl(@NotNull final Sources source, @NotNull final String treatment,
                @NotNull final String cancerType, @NotNull final String doid, @NotNull final String tumorLocationBlacklisting,
                @NotNull final EvidenceLevel level, @NotNull final EvidenceDirection direction,
                @Nullable Set<String> sourceUrls, @NotNull final Set<String> evidenceUrls) {
            this.source = source;
            this.treatment = treatment;
            this.cancerType = cancerType;
            this.doid = doid;
            this.tumorLocationBlacklisting = tumorLocationBlacklisting;
            this.level = level;
            this.direction = direction;
            this.sourceUrls = sourceUrls;
            this.evidenceUrls = evidenceUrls;
        }

        @NotNull
        @Override
        public Sources source() {
            return source;
        }

        @NotNull
        @Override
        public String treatment() {
            return treatment;
        }

        @NotNull
        @Override
        public String cancerType() {
            return cancerType;
        }

        @NotNull
        @Override
        public String doid() {
            return doid;
        }

        @NotNull
        @Override
        public String tumorLocationBlacklisting() {
            return tumorLocationBlacklisting;
        }

        @NotNull
        @Override
        public EvidenceLevel level() {
            return level;
        }

        @NotNull
        @Override
        public EvidenceDirection direction() {
            return direction;
        }

        @NotNull
        @Override
        public Set<String> sourceUrls() {
            return sourceUrls;
        }

        @NotNull
        @Override
        public Set<String> evidenceUrls() {
            return evidenceUrls;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) {
                return true;
            }
            if (o == null || getClass() != o.getClass()) {
                return false;
            }
            final ActionableEventImpl that = (ActionableEventImpl) o;
            return source == that.source && treatment.equals(that.treatment)
                    && cancerType.equals(that.cancerType) && doid.equals(that.doid) && tumorLocationBlacklisting.equals(that.tumorLocationBlacklisting)
                    && level == that.level && direction == that.direction
                    && sourceUrls.equals(that.sourceUrls) && evidenceUrls.equals(that.evidenceUrls);
        }

        @Override
        public int hashCode() {
            return Objects.hash(
                    source,
                    treatment,
                    cancerType,
                    doid,
                    tumorLocationBlacklisting,
                    level,
                    direction,
                    sourceUrls,
                    evidenceUrls);
        }

        @Override
        public String toString() {
            return "ActionableEventImpl{" + "source=" + source + ", treatment='" + treatment + '\''
                    + ", cancerType='" + cancerType + '\'' + ", doid='" + doid + '\''
                    + ", tumorLocationBlacklisting='" + tumorLocationBlacklisting + '\'' + ", level=" + level + ", direction=" + direction + ", sourceUrls="
                    + sourceUrls + ", evidenceUrls=" + evidenceUrls + '}';
        }
    }
}