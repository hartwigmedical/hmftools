package com.hartwig.hmftools.common.variant.structural;

import java.util.Collection;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public interface StructuralVariant {

    @Nullable
    Integer primaryKey();

    @NotNull
    String id();

    @Nullable
    String mateId();

    @NotNull
    StructuralVariantLeg start();

    @Nullable
    StructuralVariantLeg end();

    @Nullable
    default String chromosome(boolean isStart) {
        if (isStart) {
            return start().chromosome();
        } else {
            StructuralVariantLeg endLeg = end();
            return endLeg != null ? endLeg.chromosome() : null;
        }
    }

    @Nullable
    default Long position(boolean isStart) {
        if (isStart) {
            return start().position();
        } else {
            StructuralVariantLeg endLeg = end();
            return endLeg != null ? endLeg.position() : null;
        }
    }

    @Nullable
    default Byte orientation(boolean isStart) {
        if (isStart) {
            return start().orientation();
        } else {
            StructuralVariantLeg endLeg = end();
            return endLeg != null ? endLeg.orientation() : null;
        }
    }

    @NotNull
    String insertSequence();

    @Nullable
    String insertSequenceAlignments();

    @Nullable
    String insertSequenceRepeatClass();

    @Nullable
    String insertSequenceRepeatType();

    @Nullable
    Byte insertSequenceRepeatOrientation();

    @Nullable
    Double insertSequenceRepeatCoverage();

    @NotNull
    StructuralVariantType type();

    @Nullable
    String filter();

    @Nullable
    Boolean imprecise();

    double qualityScore();

    @Nullable
    String event();

    @Nullable
    String startLinkedBy();

    @Nullable
    String endLinkedBy();

    @NotNull
    default Collection<String> startLinks() {
        String linkedBy = startLinkedBy();
        if (linkedBy == null) {
            return Collections.emptyList();
        }
        List<String> links = Lists.newArrayList();
        for (String s : linkedBy.split(",")) {
            if (s.length() > 1) {
                links.add(s);
            }
        }
        return links;
    }

    @NotNull
    default Collection<String> endLinks() {
        String linkedBy = endLinkedBy();
        if (linkedBy == null) {
            return Collections.emptyList();
        }
        List<String> links = Lists.newArrayList();
        for (String s : linkedBy.split(",")) {
            if (s.length() > 1) {
                links.add(s);
            }
        }
        return links;
    }

    default boolean isFoldBackInversion() {
        // TODO: should we add a max bounds check in here?
        // TODO: restrict to same arm?
        StructuralVariantLeg endLeg = end();
        return endLeg != null && start().chromosome().equals(endLeg.chromosome()) && start().orientation() == endLeg.orientation();
    }

    boolean recovered();
}
