package com.hartwig.hmftools.protect.evidence;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

public final class EvidenceTestFactory {

    private EvidenceTestFactory() {
    }

    @NotNull
    public static PersonalizedEvidenceFactory createTestEvidenceFactory(@NotNull String doid) {
        return new PersonalizedEvidenceFactory(Sets.newHashSet(doid));
    }

    @NotNull
    public static PersonalizedEvidenceFactory createTestEvidenceFactory() {
        return createTestEvidenceFactory("162");
    }
}
