package com.hartwig.hmftools.serve.actionability;

import java.util.Objects;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;

import org.jetbrains.annotations.NotNull;

public final class ActionabilityTestUtil {

    public static final String TEST_ACTIONABILITY_DIR = Resources.getResource("actionability").getPath();

    private ActionabilityTestUtil() {
    }

    @NotNull
    public static ActionableEvent create(@NotNull Knowledgebase source, @NotNull String treatment, @NotNull String cancerType,
            @NotNull String doid, @NotNull EvidenceLevel level, @NotNull EvidenceDirection direction, @NotNull String url) {
        return new ActionableEventImpl(source, treatment, cancerType, doid, level, direction, url);
    }

    private static class ActionableEventImpl implements ActionableEvent {

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
        private final String url;

        public ActionableEventImpl(@NotNull final Knowledgebase source, @NotNull final String treatment, @NotNull final String cancerType,
                @NotNull final String doid, @NotNull final EvidenceLevel level, @NotNull final EvidenceDirection direction,
                @NotNull final String url) {
            this.source = source;
            this.treatment = treatment;
            this.cancerType = cancerType;
            this.doid = doid;
            this.level = level;
            this.direction = direction;
            this.url = url;
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
        public String url() {
            return url;
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
            return source == that.source && treatment.equals(that.treatment) && cancerType.equals(that.cancerType) && doid.equals(that.doid)
                    && level == that.level && direction == that.direction && url.equals(that.url);
        }

        @Override
        public int hashCode() {
            return Objects.hash(source, treatment, cancerType, doid, level, direction, url);
        }

        @Override
        public String toString() {
            return "ActionableEventImpl{" + "source=" + source + ", treatment='" + treatment + '\'' + ", cancerType='" + cancerType + '\''
                    + ", doid='" + doid + '\'' + ", level=" + level + ", direction=" + direction + ", url='" + url + '\'' + '}';
        }
    }
}
