package com.hartwig.hmftools.common.doid;

import java.util.List;
import java.util.Set;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

public class DoidParents {

    private final ListMultimap<String, String> relationship = ArrayListMultimap.create();

    public DoidParents(List<DoidEdge> edges) {
        for (DoidEdge edge : edges) {
            if (edge.predicate().equals("is_a")) {
                isA(DiseaseOntology.extractDoid(edge.subject()), DiseaseOntology.extractDoid(edge.object()));
            }
        }
    }

    public int size() {
        return relationship.size();
    }

    @NotNull
    public Set<String> parents(@NotNull String child) {
        Set<String> result = Sets.newHashSet();
        inner(child, result);
        return result;
    }

    private void isA(@NotNull String child, @NotNull String parent) {
        if (relationship.containsKey(child)) {
            List<String> parents = relationship.get(child);
            if (!parents.contains(parent)) {
                parents.add(parent);
            }
        } else {
            relationship.put(child, parent);
        }
    }

    private void inner(@NotNull String child, @NotNull final Set<String> result) {
        if (!relationship.containsKey(child)) {
            return;
        }

        for (String parent : relationship.get(child)) {
            if (result.add(parent)) {
                inner(parent, result);
            }
        }
    }

}
