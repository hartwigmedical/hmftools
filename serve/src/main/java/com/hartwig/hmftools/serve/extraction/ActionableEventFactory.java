package com.hartwig.hmftools.serve.extraction;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;
import com.hartwig.hmftools.serve.actionability.fusion.ActionableFusion;
import com.hartwig.hmftools.serve.actionability.fusion.ImmutableActionableFusion;
import com.hartwig.hmftools.serve.actionability.gene.ActionableGene;
import com.hartwig.hmftools.serve.actionability.gene.ImmutableActionableGene;
import com.hartwig.hmftools.serve.actionability.hotspot.ActionableHotspot;
import com.hartwig.hmftools.serve.actionability.hotspot.ImmutableActionableHotspot;
import com.hartwig.hmftools.serve.actionability.range.ActionableRange;
import com.hartwig.hmftools.serve.actionability.range.ImmutableActionableRange;
import com.hartwig.hmftools.serve.actionability.signature.ActionableSignature;
import com.hartwig.hmftools.serve.actionability.signature.ImmutableActionableSignature;
import com.hartwig.hmftools.serve.extraction.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.extraction.fusion.KnownFusionPair;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelAnnotation;
import com.hartwig.hmftools.serve.extraction.gene.GeneLevelEvent;
import com.hartwig.hmftools.serve.extraction.range.RangeAnnotation;
import com.hartwig.hmftools.serve.extraction.signature.SignatureName;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ActionableEventFactory {

    private ActionableEventFactory() {
    }

    @NotNull
    public static Set<ActionableHotspot> toActionableHotspots(@NotNull ActionableEvent actionableEvent,
            @Nullable List<VariantHotspot> hotspotList) {
        if (hotspotList == null) {
            return Sets.newHashSet();
        }

        Set<ActionableHotspot> actionableHotspots = Sets.newHashSet();
        for (VariantHotspot hotspot : hotspotList) {
            actionableHotspots.add(ImmutableActionableHotspot.builder()
                    .from(actionableEvent)
                    .chromosome(hotspot.chromosome())
                    .position(hotspot.position())
                    .ref(hotspot.ref())
                    .alt(hotspot.alt())
                    .build());
        }

        return actionableHotspots;
    }

    @NotNull
    public static Set<ActionableRange> toActionableRanges(@NotNull ActionableEvent actionableEvent,
            @Nullable List<? extends RangeAnnotation> ranges) {
        if (ranges == null) {
            return Sets.newHashSet();
        }

        Set<ActionableRange> actionableRanges = Sets.newHashSet();
        for (RangeAnnotation range : ranges) {
            actionableRanges.add(ImmutableActionableRange.builder()
                    .from(actionableEvent)
                    .chromosome(range.chromosome())
                    .start(range.start())
                    .end(range.end())
                    .gene(range.gene())
                    .mutationType(range.mutationType())
                    .build());
        }

        return actionableRanges;
    }

    @NotNull
    public static ActionableGene copyNumberToActionableGene(@NotNull ActionableEvent actionableEvent, @NotNull KnownCopyNumber copyNumber) {
        GeneLevelEvent event;
        switch (copyNumber.type()) {
            case AMPLIFICATION: {
                event = GeneLevelEvent.AMPLIFICATION;
                break;
            }
            case DELETION: {
                event = GeneLevelEvent.DELETION;
                break;
            }
            default:
                throw new IllegalStateException("Invalid copy number type: " + copyNumber.type());
        }

        return ImmutableActionableGene.builder().from(actionableEvent).gene(copyNumber.gene()).event(event).build();
    }

    @NotNull
    public static ActionableGene geneLevelEventToActionableGene(@NotNull ActionableEvent actionableEvent,
            @NotNull GeneLevelAnnotation geneLevelEvent) {
        return ImmutableActionableGene.builder().from(actionableEvent).gene(geneLevelEvent.gene()).event(geneLevelEvent.event()).build();
    }

    @NotNull
    public static ActionableFusion toActionableFusion(@NotNull ActionableEvent actionableEvent, @NotNull KnownFusionPair fusion) {
        return ImmutableActionableFusion.builder()
                .from(actionableEvent)
                .geneUp(fusion.geneUp())
                .minExonUp(fusion.minExonUp())
                .maxExonUp(fusion.maxExonUp())
                .geneDown(fusion.geneDown())
                .minExonDown(fusion.minExonDown())
                .maxExonDown(fusion.maxExonDown())
                .build();
    }

    @NotNull
    public static ActionableSignature toActionableSignature(@NotNull ActionableEvent actionableEvent, @NotNull SignatureName signature) {
        return ImmutableActionableSignature.builder().from(actionableEvent).name(signature).build();
    }
}
