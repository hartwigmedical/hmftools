package com.hartwig.hmftools.svannotation;

import java.util.List;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;

public class BreakendAnnotations {

    private final StructuralVariantAnnotation parent;
    private final String chromosome;
    private final long position;
    private final int orientation;
    private final List<GeneAnnotation> genes = Lists.newArrayList();

    BreakendAnnotations(final StructuralVariantAnnotation parent) {
        this.parent = parent;
        this.chromosome = "ERR";
        this.position = 0;
        this.orientation = 0;
    }

    BreakendAnnotations(final StructuralVariantAnnotation parent, final String chromosome, final long position, final int orientation) {
        this.parent = parent;
        this.chromosome = chromosome;
        this.position = position;
        this.orientation = orientation;
    }

    void addGeneAnnotation(final GeneAnnotation gene) {
        genes.add(gene);
    }

    public StructuralVariantAnnotation getStructuralVariant() {
        return parent;
    }

    public String getChromosome() {
        return chromosome;
    }

    public long getPosition() {
        return position;
    }

    public int getOrientation() {
        return orientation;
    }

    public String getPositionString() {
        return String.format("chr%s:%d",chromosome, position);
    }

    public List<GeneAnnotation> getGenes() {
        return ImmutableList.copyOf(genes);
    }
}
