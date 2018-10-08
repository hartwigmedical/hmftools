package com.hartwig.hmftools.svanalysis.svgraph;

import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import org.jetbrains.annotations.NotNull;

import java.util.Collection;

public class BreakendConsistency {
    private final BgSegment anchor;
    private final BgSegment refPath;
    private final EnrichedStructuralVariant sv;
    private final Collection<EnrichedStructuralVariant> alternateSvs;
    private final Collection<EnrichedStructuralVariant> oppositeOrientationSvs;
    private final double nominalPloidy;

    public BreakendConsistency(
            double ploidy,
            @NotNull BgSegment anchor, BgSegment refPath,
            @NotNull EnrichedStructuralVariant sv,
            @NotNull Collection<EnrichedStructuralVariant> alternateSvs,
            @NotNull Collection<EnrichedStructuralVariant> oppositeOrientationSvs) {
        this.nominalPloidy = ploidy;
        this.anchor = anchor;
        this.refPath = refPath;
        this.sv = sv;
        this.alternateSvs = alternateSvs;
        this.oppositeOrientationSvs = oppositeOrientationSvs;
    }

    /**
     * @return notional ploidy of event
     */
    public double ploidy() {
        return nominalPloidy;
    }

    /**
     * @return Difference in copy number between the flanking segment and the notional ploidy of the event
     */
    public double copyNumberDelta() {
        return anchorPloidy() - referencePathPloidy() - ploidy();
    }

    /**
     * Difference in ploidy between event support and notional ploidy
     * @return
     */
    public double eventDelta() {
        return ploidy() - sv.ploidy();
    }
    public double anchorPloidy() { return anchor == null ? 0 : anchor.ploidy(); }
    public double referencePathPloidy() { return refPath == null ? 0 : refPath.ploidy(); }
    public EnrichedStructuralVariant sv() { return sv; }

    /**
     * Ploidy of other SVs. SVs in the opposite orientation have negative ploidy.
     * @return
     */
    public double otherSvPloidy() {
        return alternateSvs.stream().mapToDouble(x -> x.ploidy()).sum()
            - oppositeOrientationSvs.stream().mapToDouble(x -> x.ploidy()).sum();
    }
    @Override
    public String toString() {
        return String.format("%s ploidy: %.2f ∆cn: %.2f ∆sv: %.2f ∆other: %.2f", sv.id(), ploidy(), copyNumberDelta(), eventDelta(), otherSvPloidy());
    }
}
