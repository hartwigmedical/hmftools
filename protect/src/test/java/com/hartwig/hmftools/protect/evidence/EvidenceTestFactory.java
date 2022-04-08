package com.hartwig.hmftools.protect.evidence;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.doid.DoidEdge;
import com.hartwig.hmftools.common.doid.DoidParents;
import com.hartwig.hmftools.common.doid.DoidParentsTest;

import org.jetbrains.annotations.NotNull;

public final class EvidenceTestFactory {

    private EvidenceTestFactory() {
    }

    @NotNull
    public static PersonalizedEvidenceFactory createTestEvidenceFactory() {

        return new PersonalizedEvidenceFactory(Sets.newHashSet("162"));
    }
}
