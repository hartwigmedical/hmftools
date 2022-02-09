package com.hartwig.hmftools.serve.actionability;

import java.util.Objects;
import java.util.Set;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ActionabilityTestUtil {

    public static final String TEST_SERVE_OUTPUT_DIR = Resources.getResource("serve_output").getPath();

    private ActionabilityTestUtil() {
    }

    @NotNull
    public static ActionableEvent create(@NotNull String rawInput, @NotNull Knowledgebase source, @NotNull String treatment,
            @NotNull String cancerType, @NotNull String doid, @NotNull String blacklistCancerType, @NotNull String blacklistedDoid,
            @NotNull EvidenceLevel level, @NotNull EvidenceDirection direction, @Nullable Set<String> urlSource,
            @NotNull Set<String> urls) {
        return new ActionableEventImpl(rawInput,
                source,
                treatment,
                cancerType,
                doid,
                blacklistCancerType,
                blacklistedDoid,
                level,
                direction,
                urlSource,
                urls);
    }

    private static class ActionableEventImpl implements ActionableEvent {

        @NotNull
        private final String rawInput;
        @NotNull
        private final Knowledgebase source;
        @NotNull
        private final String treatment;
        @NotNull
        private final String cancerType;
        @NotNull
        private final String doid;
        @NotNull
        private final String blacklistCancerType;
        @NotNull
        private final String blacklistedDoid;
        @NotNull
        private final EvidenceLevel level;
        @NotNull
        private final EvidenceDirection direction;
        @Nullable
        private final Set<String> urlSource;
        @NotNull
        private final Set<String> urls;

        public ActionableEventImpl(@NotNull String rawInput, @NotNull final Knowledgebase source, @NotNull final String treatment,
                @NotNull final String cancerType, @NotNull final String doid, @NotNull final String blacklistCancerType,
                @NotNull final String blacklistedDoid, @NotNull final EvidenceLevel level, @NotNull final EvidenceDirection direction,
                @Nullable Set<String> urlSource, @NotNull final Set<String> urls) {
            this.rawInput = rawInput;
            this.source = source;
            this.treatment = treatment;
            this.cancerType = cancerType;
            this.doid = doid;
            this.blacklistCancerType = blacklistCancerType;
            this.blacklistedDoid = blacklistedDoid;
            this.level = level;
            this.direction = direction;
            this.urlSource = urlSource;
            this.urls = urls;
        }

        @NotNull
        @Override
        public String rawInput() {
            return rawInput;
        }

        @NotNull
        @Override
        public Knowledgebase source() {
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
        public String blacklistCancerType() {
            return blacklistCancerType;
        }

        @NotNull
        @Override
        public String blacklistedDoid() {
            return blacklistedDoid;
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
        public Set<String> urlSource() {
            return urlSource;
        }

        @NotNull
        @Override
        public Set<String> urls() {
            return urls;
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
            return rawInput.equals(that.rawInput()) && source == that.source && treatment.equals(that.treatment)
                    && cancerType.equals(that.cancerType) && doid.equals(that.doid) && blacklistCancerType.equals(that.blacklistCancerType)
                    && blacklistedDoid.equals(that.blacklistedDoid) && level == that.level && direction == that.direction
                    && urlSource.equals(that.urlSource) && urls.equals(that.urls);
        }

        @Override
        public int hashCode() {
            return Objects.hash(rawInput,
                    source,
                    treatment,
                    cancerType,
                    doid,
                    blacklistCancerType,
                    blacklistedDoid,
                    level,
                    direction,
                    urlSource,
                    urls);
        }

        @Override
        public String toString() {
            return "ActionableEventImpl{" + "rawInput=" + rawInput + ",source=" + source + ", treatment='" + treatment + '\''
                    + ", cancerType='" + cancerType + '\'' + ", doid='" + doid + ", blacklistCancerType='" + blacklistCancerType + '\''
                    + ", blacklistedDoid='" + blacklistedDoid + '\'' + ", level=" + level + ", direction=" + direction + ", urlSource="
                    + urlSource + ", urls=" + urls + '}';
        }
    }
}
