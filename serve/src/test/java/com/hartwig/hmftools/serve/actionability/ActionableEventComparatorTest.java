package com.hartwig.hmftools.serve.actionability;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.serve.tumorlocation.ImmutableTumorLocation;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ActionableEventComparatorTest {

    @Test
    public void canSortActionableEvents() {
        ActionableEvent event1 = create(Knowledgebase.VICC_CGI, "DrugA", "CancerA", EvidenceLevel.A, EvidenceDirection.RESISTANT);
        ActionableEvent event2 = create(Knowledgebase.VICC_CGI, "DrugA", "CancerA", EvidenceLevel.A, EvidenceDirection.RESPONSIVE);
        ActionableEvent event3 = create(Knowledgebase.VICC_CGI, "DrugB", "CancerA", EvidenceLevel.A, EvidenceDirection.RESPONSIVE);
        ActionableEvent event4 = create(Knowledgebase.VICC_CGI, "DrugB", "CancerB", EvidenceLevel.A, EvidenceDirection.RESPONSIVE);
        ActionableEvent event5 = create(Knowledgebase.VICC_CGI, "DrugA", "CancerA", EvidenceLevel.B, EvidenceDirection.RESISTANT);
        ActionableEvent event6 = create(Knowledgebase.VICC_CIVIC, "DrugA", "CancerA", EvidenceLevel.A, EvidenceDirection.RESISTANT);

        List<ActionableEvent> events = Lists.newArrayList(event3, event5, event1, event6, event4, event2);
        events.sort(new ActionableEventComparator());

        assertEquals(event1, events.get(0));
        assertEquals(event2, events.get(1));
        assertEquals(event3, events.get(2));
        assertEquals(event4, events.get(3));
        assertEquals(event5, events.get(4));
        assertEquals(event6, events.get(5));
    }

    @NotNull
    private static ActionableEvent create(@NotNull Knowledgebase source, @NotNull String treatment, @NotNull String whiteListCancerType,
            @NotNull EvidenceLevel level, @NotNull EvidenceDirection direction) {
        return ActionabilityTestUtil.create(source,
                "rawInput",
                Sets.newHashSet(),
                treatment,
                ImmutableTumorLocation.builder()
                        .cancerType(whiteListCancerType)
                        .doid("whitlist doid")
                        .build(),
                Sets.newHashSet(ImmutableTumorLocation.builder()
                        .cancerType("blacklist cancertype")
                        .doid("blacklist doid")
                        .build()),
                level,
                direction,
                Sets.newHashSet());
    }
}