package com.hartwig.hmftools.patientdb.clinical.curators;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.doid.DoidNode;
import org.jetbrains.annotations.NotNull;

import java.util.List;

public class DoidNodesResolver {

    @NotNull
    private final List<DoidNode> doidNodes;

    public DoidNodesResolver(@NotNull List<DoidNode> doidNodes) {
        this.doidNodes = doidNodes;
    }

    @NotNull
    @VisibleForTesting
    public List<DoidNode> resolveDoidNodes(@NotNull List<String> doidsToResolve) {
        List<DoidNode> resolvedDoidNodes = Lists.newArrayList();
        for (String doid : doidsToResolve) {
            for (DoidNode doidNode : doidNodes) {
                if (doidNode.doid().equals(doid)) {
                    resolvedDoidNodes.add(doidNode);
                }
            }
        }
        return resolvedDoidNodes;
    }
}
