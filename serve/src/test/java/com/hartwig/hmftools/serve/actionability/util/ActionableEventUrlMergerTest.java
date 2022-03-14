package com.hartwig.hmftools.serve.actionability.util;

import static org.junit.Assert.assertEquals;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.common.serve.actionability.EvidenceLevel;
import com.hartwig.hmftools.serve.actionability.ActionabilityTestUtil;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusion;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusionUrlConsolidator;
import com.hartwig.hmftools.serve.actionability.fusion.ImmutableActionableFusion;
import com.hartwig.hmftools.serve.cancertype.ImmutableCancerType;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ActionableEventUrlMergerTest {

    private static final String GENE_1 = "gene1";
    private static final String GENE_2 = "gene2";

    @Test
    public void canMergeFusionUrls() {
        Set<ActionableFusion> fusions = Sets.newHashSet();
        fusions.add(createFusion(GENE_1, "url1"));
        fusions.add(createFusion(GENE_1, "url2"));
        fusions.add(createFusion(GENE_2, "url3"));

        Set<ActionableFusion> consolidated = ActionableEventUrlMerger.merge(fusions, new ActionableFusionUrlConsolidator());
        assertEquals(2, consolidated.size());
        assertEquals(Sets.newHashSet("url1", "url2"), findByGeneUp(consolidated, GENE_1).evidenceUrls());
        assertEquals(Sets.newHashSet("url3"), findByGeneUp(consolidated, GENE_2).evidenceUrls());
    }

    @NotNull
    private static ActionableFusion createFusion(@NotNull String gene, @NotNull String url) {
        return ImmutableActionableFusion.builder()
                .from(ActionabilityTestUtil.create(Knowledgebase.VICC_CGI,
                        "rawInput",
                        Sets.newHashSet(),
                        "treatment",
                        ImmutableCancerType.builder()
                                .cancerType("applicable cancerType")
                                .doid("applicable doid")
                                .build(),
                        Sets.newHashSet(ImmutableCancerType.builder()
                                .cancerType("blacklist cancerType")
                                .doid("blacklist doid")
                                .build()),
                        EvidenceLevel.A,
                        EvidenceDirection.RESPONSIVE,
                        Sets.newHashSet(url)))
                .geneUp(gene)
                .geneDown(gene)
                .build();
    }

    @NotNull
    private static ActionableFusion findByGeneUp(@NotNull Iterable<ActionableFusion> fusions, @NotNull String gene) {
        for (ActionableFusion fusion : fusions) {
            if (fusion.geneUp().equals(gene)) {
                return fusion;
            }
        }

        throw new IllegalStateException("Could not find gene in fusions: " + gene);
    }
}