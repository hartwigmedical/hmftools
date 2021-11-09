package com.hartwig.hmftools.common.doid;

import java.util.List;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

public class DoidParents {

    @NotNull
    private final ListMultimap<String, String> relationship;

    @NotNull
    public static DoidParents fromEdges(@NotNull List<DoidEdge> edges) {
        ListMultimap<String, String> relationship = ArrayListMultimap.create();
        for (DoidEdge edge : edges) {
            if (edge.predicate().equals("is_a")) {
                String child = DiseaseOntology.extractDoid(edge.subject());
                String parent = DiseaseOntology.extractDoid(edge.object());

                if (relationship.containsKey(child)) {
                    List<String> parents = relationship.get(child);
                    if (!parents.contains(parent)) {
                        parents.add(parent);
                    }
                } else {
                    relationship.put(child, parent);
                }
            }
        }
        return new DoidParents(relationship);
    }

    @VisibleForTesting
    DoidParents(@NotNull final ListMultimap<String, String> relationship) {
        this.relationship = relationship;
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

    private void inner(@NotNull String child, @NotNull Set<String> result) {
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
