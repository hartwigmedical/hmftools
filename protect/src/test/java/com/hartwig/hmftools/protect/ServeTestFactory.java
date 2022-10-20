package com.hartwig.hmftools.protect;

import java.util.Objects;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.common.serve.actionability.ImmutableTreatment;
import com.hartwig.hmftools.common.serve.actionability.Treatment;
import com.hartwig.hmftools.common.serve.actionability.characteristic.ActionableCharacteristic;
import com.hartwig.hmftools.common.serve.actionability.characteristic.ImmutableActionableCharacteristic;
import com.hartwig.hmftools.common.serve.actionability.fusion.ActionableFusion;
import com.hartwig.hmftools.common.serve.actionability.fusion.ImmutableActionableFusion;
import com.hartwig.hmftools.common.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.common.serve.actionability.gene.ImmutableActionableGene;
import com.hartwig.hmftools.common.serve.actionability.hotspot.ActionableHotspot;
import com.hartwig.hmftools.common.serve.actionability.hotspot.ImmutableActionableHotspot;
import com.hartwig.hmftools.common.serve.actionability.immuno.ActionableHLA;
import com.hartwig.hmftools.common.serve.actionability.immuno.ImmutableActionableHLA;
import com.hartwig.hmftools.common.serve.actionability.range.ActionableRange;
import com.hartwig.hmftools.common.serve.actionability.range.ImmutableActionableRange;
import com.hartwig.hmftools.common.serve.actionability.range.RangeType;
import com.hartwig.hmftools.common.serve.cancertype.CancerType;
import com.hartwig.hmftools.common.serve.cancertype.ImmutableCancerType;
import com.hartwig.hmftools.common.serve.datamodel.MutationTypeFilter;
import com.hartwig.hmftools.common.serve.datamodel.characteristic.TumorCharacteristicAnnotation;
import com.hartwig.hmftools.common.serve.datamodel.gene.GeneLevelEvent;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ServeTestFactory {

    private ServeTestFactory() {
    }

    @NotNull
    public static ActionableHotspot createTestActionableHotspot() {
        return ImmutableActionableHotspot.builder()
                .from(createTestBaseEvent())
                .chromosome(Strings.EMPTY)
                .position(0)
                .ref(Strings.EMPTY)
                .alt(Strings.EMPTY)
                .build();
    }

    @NotNull
    public static ActionableRange createTestActionableRange() {
        return ImmutableActionableRange.builder()
                .from(createTestBaseEvent())
                .gene(Strings.EMPTY)
                .transcript(Strings.EMPTY)
                .chromosome(Strings.EMPTY)
                .start(0)
                .end(0)
                .mutationType(MutationTypeFilter.ANY)
                .rangeType(RangeType.EXON)
                .rank(0)
                .build();
    }

    @NotNull
    public static ActionableGene createTestActionableGene() {
        return ImmutableActionableGene.builder()
                .from(createTestBaseEvent())
                .gene(Strings.EMPTY)
                .event(GeneLevelEvent.AMPLIFICATION)
                .build();
    }

    @NotNull
    public static ActionableFusion createTestActionableFusion() {
        return ImmutableActionableFusion.builder().from(createTestBaseEvent()).geneUp(Strings.EMPTY).geneDown(Strings.EMPTY).build();
    }

    @NotNull
    public static ActionableCharacteristic createTestActionableCharacteristic() {
        return ImmutableActionableCharacteristic.builder()
                .from(createTestBaseEvent())
                .name(TumorCharacteristicAnnotation.MICROSATELLITE_UNSTABLE)
                .build();
    }

    @NotNull
    public static ActionableHLA createTestActionableHLA() {
        return ImmutableActionableHLA.builder().from(createTestBaseEvent()).hlaType(Strings.EMPTY).build();
    }

    @NotNull
    private static ActionableEvent createTestBaseEvent() {
        return createTestBaseEvent(Knowledgebase.HARTWIG_CURATED);
    }

    @NotNull
    private static ActionableEvent createTestBaseEvent(@NotNull Knowledgebase source) {
        return create(source,
                "source event",
                Sets.newHashSet(),
                ImmutableTreatment.builder()
                        .treament("treatment")
                        .sourceRelevantTreatmentApproaches(Sets.newHashSet("drug classes"))
                        .relevantTreatmentApproaches(Sets.newHashSet("drug classes"))
                        .build(),
                ImmutableCancerType.builder().name("applicable name").doid("applicable doid").build(),
                Sets.newHashSet(ImmutableCancerType.builder().name("blacklist name").doid("blacklist doid").build()),
                EvidenceLevel.A,
                EvidenceDirection.RESPONSIVE,
                Sets.newHashSet());
    }

    @NotNull
    public static ActionableEvent create(@NotNull Knowledgebase source, @NotNull String sourceEvent, @NotNull Set<String> sourceUrls,
            @NotNull Treatment treatment, @NotNull CancerType applicableCancerType, @NotNull Set<CancerType> blacklistCancerTypes,
            @NotNull EvidenceLevel level, @NotNull EvidenceDirection direction, @NotNull Set<String> evidenceUrls) {
        return new ActionableEventImpl(source,
                sourceEvent,
                sourceUrls,
                treatment,
                applicableCancerType,
                blacklistCancerTypes,
                level,
                direction,
                evidenceUrls);
    }

    private static class ActionableEventImpl implements ActionableEvent {

        @NotNull
        private final Knowledgebase source;
        @NotNull
        private final String sourceEvent;
        @NotNull
        private final Set<String> sourceUrls;
        @NotNull
        private final Treatment treatment;
        @NotNull
        private final CancerType applicableCancerType;
        @NotNull
        private final Set<CancerType> blacklistCancerTypes;
        @NotNull
        private final EvidenceLevel level;
        @NotNull
        private final EvidenceDirection direction;
        @NotNull
        private final Set<String> evidenceUrls;

        public ActionableEventImpl(@NotNull Knowledgebase source, @NotNull String sourceEvent, @NotNull Set<String> sourceUrls,
                @NotNull Treatment treatment, @NotNull CancerType applicableCancerType, @NotNull Set<CancerType> blacklistCancerTypes,
                @NotNull EvidenceLevel level, @NotNull EvidenceDirection direction, @NotNull Set<String> evidenceUrls) {
            this.source = source;
            this.sourceEvent = sourceEvent;
            this.sourceUrls = sourceUrls;
            this.treatment = treatment;
            this.applicableCancerType = applicableCancerType;
            this.blacklistCancerTypes = blacklistCancerTypes;
            this.level = level;
            this.direction = direction;
            this.evidenceUrls = evidenceUrls;
        }

        @NotNull
        @Override
        public Knowledgebase source() {
            return source;
        }

        @NotNull
        @Override
        public String sourceEvent() {
            return sourceEvent;
        }

        @NotNull
        @Override
        public Set<String> sourceUrls() {
            return sourceUrls;
        }

        @NotNull
        @Override
        public Treatment treatment() {
            return treatment;
        }

        @NotNull
        @Override
        public CancerType applicableCancerType() {
            return applicableCancerType;
        }

        @NotNull
        @Override
        public Set<CancerType> blacklistCancerTypes() {
            return blacklistCancerTypes;
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
            return source == that.source && sourceEvent.equals(that.sourceEvent) && sourceUrls.equals(that.sourceUrls) && treatment.equals(
                    that.treatment) && applicableCancerType.equals(that.applicableCancerType)
                    && blacklistCancerTypes.equals(that.blacklistCancerTypes) && level == that.level && direction == that.direction
                    && evidenceUrls.equals(that.evidenceUrls);
        }

        @Override
        public int hashCode() {
            return Objects.hash(source,
                    sourceEvent,
                    sourceUrls,
                    treatment,
                    applicableCancerType,
                    blacklistCancerTypes,
                    level,
                    direction,
                    evidenceUrls);
        }

        @Override
        public String toString() {
            return "ActionableEventImpl{" + "source=" + source + ", sourceEvent='" + sourceEvent + '\'' + ", sourceUrls=" + sourceUrls
                    + ", treatment=" + treatment + ", applicableCancerType=" + applicableCancerType + ", blacklistCancerTypes="
                    + blacklistCancerTypes + ", level=" + level + ", direction=" + direction + ", evidenceUrls=" + evidenceUrls + '}';
        }
    }
}
