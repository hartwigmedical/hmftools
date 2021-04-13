package com.hartwig.hmftools.ckb.classification;

import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.ckb.datamodel.variant.Variant;

import org.jetbrains.annotations.NotNull;

public final class EventExtractorCuration {

    private EventExtractorCuration() {
    }

    @NotNull
    public static String extractEvent(@NotNull CkbEntry entry) {
        Variant variant = entry.variants().get(0);
        String event;

        if (variant.variant().equals("fusion") && variant.impact() != null && variant.impact().equals("fusion")) {
            event = "fusion promiscuous";
        } else if (variant.impact() != null && variant.impact().equals("fusion")) {
            event = variant.variant().replaceAll("\\s+", "");
        } else if (variant.variant().contains("exon")) {
            event = variant.variant().replace("exon", "exon ");
        } else {
            event = variant.variant();
        }
        return event;
    }

    @NotNull
    public static String extractGene(@NotNull CkbEntry entry) {
        Variant variant = entry.variants().get(0);
        String geneSymbol;

        if (variant.variant().equals("fusion") && variant.impact() != null && variant.impact().equals("fusion")) {
            geneSymbol = variant.gene().geneSymbol();
        } else if (variant.impact() != null && variant.impact().equals("fusion")) {
            geneSymbol = variant.fullName();
        } else if (variant.variant().contains("exon")) {
            geneSymbol = variant.gene().geneSymbol();
        } else {
            geneSymbol = variant.gene().geneSymbol();
        }
        return geneSymbol;
    }
}
