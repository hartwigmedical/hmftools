package com.hartwig.hmftools.protect.evidence;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

public final class EvidenceTestFactory {

    private EvidenceTestFactory() {
    }

    @NotNull
    public static PersonalizedEvidenceFactory create(@NotNull String doid) {
        return new PersonalizedEvidenceFactory(Sets.newHashSet(doid));
    }

    @NotNull
    public static PersonalizedEvidenceFactory create() {
        return create("162");
    }
}
