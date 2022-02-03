package com.hartwig.hmftools.serve.actionability;

import java.util.Objects;
import java.util.Set;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;

import org.jetbrains.annotations.NotNull;

public final class ActionabilityTestUtil {

    public static final String TEST_SERVE_OUTPUT_DIR = Resources.getResource("serve_output").getPath();

    private ActionabilityTestUtil() {
    }

    @NotNull
    public static ActionableEvent create(@NotNull String rawInput, @NotNull Knowledgebase source, @NotNull String treatment, @NotNull String cancerType,
            @NotNull String doid, @NotNull EvidenceLevel level, @NotNull EvidenceDirection direction, @NotNull Set<String> urls) {
        return new ActionableEventImpl(rawInput, source, treatment, cancerType, doid, level, direction, urls);
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
        private final EvidenceLevel level;
        @NotNull
        private final EvidenceDirection direction;
        @NotNull
        private final Set<String> urls;

        public ActionableEventImpl(@NotNull String rawInput, @NotNull final Knowledgebase source, @NotNull final String treatment, @NotNull final String cancerType,
                @NotNull final String doid, @NotNull final EvidenceLevel level, @NotNull final EvidenceDirection direction,
                @NotNull final Set<String> urls) {
            this.rawInput = rawInput;
            this.source = source;
            this.treatment = treatment;
            this.cancerType = cancerType;
            this.doid = doid;
            this.level = level;
            this.direction = direction;
            this.urls = urls;
        }

        @NotNull
        @Override
        public String rawInput() { return  rawInput; }
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
            return rawInput.equals(that.rawInput()) && source == that.source && treatment.equals(that.treatment) && cancerType.equals(that.cancerType) && doid.equals(that.doid)
                    && level == that.level && direction == that.direction && urls.equals(that.urls);
        }

        @Override
        public int hashCode() {
            return Objects.hash(rawInput, source, treatment, cancerType, doid, level, direction, urls);
        }

        @Override
        public String toString() {
            return "ActionableEventImpl{" + "rawInput=" + rawInput + ",source=" + source + ", treatment='" + treatment + '\'' + ", cancerType='" + cancerType + '\''
                    + ", doid='" + doid + '\'' + ", level=" + level + ", direction=" + direction + ", urls=" + urls + '}';
        }
    }
}
