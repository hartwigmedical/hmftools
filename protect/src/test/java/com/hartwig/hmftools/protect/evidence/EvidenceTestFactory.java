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

        List<DoidEdge> edges = Lists.newArrayList();
        edges.add(DoidParentsTest.createParent("299", "305"));
        edges.add(DoidParentsTest.createParent("305", "162"));
        edges.add(DoidParentsTest.createEdge("305", "has_a", "162"));

        DoidParents victim = DoidParents.fromEdges(edges);

        return new PersonalizedEvidenceFactory(Sets.newHashSet(), victim);
    }
}
