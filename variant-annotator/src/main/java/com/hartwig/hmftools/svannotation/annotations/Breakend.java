package com.hartwig.hmftools.svannotation.annotations;

import java.util.List;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;

public class Breakend {

    private final StructuralVariantAnnotation parent;
    private final String chromosome;
    private final long position;
    private final int orientation;
    private final Double alleleFrequency;
    private final List<GeneAnnotation> geneAnnotations = Lists.newArrayList();

    public Breakend(final StructuralVariantAnnotation parent) {
        this.parent = parent;
        this.chromosome = "ERR";
        this.position = 0;
        this.orientation = 0;
        this.alleleFrequency = null;
    }

    public Breakend(final StructuralVariantAnnotation parent, final String chromosome, final long position, final int orientation,
            final Double alleleFrequency) {
        this.parent = parent;
        this.chromosome = chromosome;
        this.position = position;
        this.orientation = orientation;
        this.alleleFrequency = alleleFrequency;
    }

    public void addGeneAnnotation(final GeneAnnotation geneAnnotation) {
        geneAnnotations.add(geneAnnotation);
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

    public Double getAlleleFrequency() {
        return alleleFrequency;
    }

    public String getPositionString() {
        return String.format("chr%s:%d", chromosome, position);
    }

    public List<GeneAnnotation> getGeneAnnotations() {
        return ImmutableList.copyOf(geneAnnotations);
    }
}
