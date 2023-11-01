package com.hartwig.hmftools.patientdb.clinical.curators;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.doid.DoidNode;
import org.jetbrains.annotations.NotNull;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class DoidNodesResolver {

    @NotNull
    private final Map<String, DoidNode> doidNodes = new HashMap<>();

    public DoidNodesResolver(@NotNull List<DoidNode> doidNodes)
    {
        for(final DoidNode doidNode : doidNodes)
        {
            this.doidNodes.put(doidNode.doid(), doidNode);
        }
    }

    @NotNull
    public List<DoidNode> resolveDoidNodes(@NotNull List<String> doidsToResolve)
    {
        List<DoidNode> resolvedDoidNodes = Lists.newArrayList();
        for(final String doid : doidsToResolve)
        {
            DoidNode doidNode = doidNodes.get(doid);
            if(doidNode != null)
            {
                resolvedDoidNodes.add(doidNode);
            }
        }
        return resolvedDoidNodes;
    }
}
