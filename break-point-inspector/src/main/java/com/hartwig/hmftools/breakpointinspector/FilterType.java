package com.hartwig.hmftools.breakpointinspector;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.vcf.VCFFilterHeaderLine;

enum FilterType {
    BREAKPOINT_ERROR("BPI_BreakpointError", "BPI failed to determine breakpoints"),
    MIN_DEPTH("BPI_MinDepth", "The depth across one of the breakpoints is <10"),
    MIN_ANCHOR_LENGTH("BPI_MinAnchorLength", "There isn't at least one PR with >=30 bases matched in both alignments"),
    SR_SUPPORT_ZERO("BPI_SRSupportZero", "Short delete or dupe (<1000) must have SR support"),
    SR_NORMAL_SUPPORT("BPI_SRNormalSupport", "Short delete or dupe (<1000) has SR support in normal"),
    PR_NORMAL_SUPPORT("BPI_PRNormalSupport", "PR support in the normal"),
    PR_SUPPORT_ZERO("BPI_PRSupportZero", "No PR support in tumor"),
    CLIPPING_CONCORDANCE("BPI_ClippingConcordance", "At least 5 base clipped bases concordance between tumor and normal");

    @NotNull
    private final String name;
    @NotNull
    private final String description;

    FilterType(@NotNull final String name, @NotNull final String description) {
        this.name = name;
        this.description = description;
    }

    @NotNull
    VCFFilterHeaderLine toHeaderLine() {
        return new VCFFilterHeaderLine(name, description);
    }

    @Override
    public String toString() {
        return name;
    }
}
