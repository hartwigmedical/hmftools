package com.hartwig.hmftools.common.cosmic;

import org.jetbrains.annotations.NotNull;

public class CosmicData {
    @NotNull
    private final String description;
    @NotNull
    private final String entrezId;
    @NotNull
    private final String genomeLocation;
    @NotNull
    private final String chromosomeBand;
    @NotNull
    private final String somatic;
    @NotNull
    private final String germline;
    @NotNull
    private final String somaticTumorTypes;
    @NotNull
    private final String germlineTumorTypes;
    @NotNull
    private final String cancerSyndrome;
    @NotNull
    private final String tissueType;
    @NotNull
    private final String molecularGenetics;
    @NotNull
    private final String role;
    @NotNull
    private final String mutationTypes;
    @NotNull
    private final String translocationPartner;
    @NotNull
    private final String otherGermlineMut;
    @NotNull
    private final String otherSyndrome;

    public CosmicData(@NotNull final String description, @NotNull final String entrezId,
            @NotNull final String genomeLocation, @NotNull final String chromosomeBand, @NotNull final String somatic,
            @NotNull final String germline, @NotNull final String somaticTumorTypes,
            @NotNull final String germlineTumorTypes, @NotNull final String cancerSyndrome,
            @NotNull final String tissueType, @NotNull final String molecularGenetics, @NotNull final String role,
            @NotNull final String mutationTypes, @NotNull final String translocationPartner,
            @NotNull final String otherGermlineMut, @NotNull final String otherSyndrome) {
        this.description = description;
        this.entrezId = entrezId;
        this.genomeLocation = genomeLocation;
        this.chromosomeBand = chromosomeBand;
        this.somatic = somatic;
        this.germline = germline;
        this.somaticTumorTypes = somaticTumorTypes;
        this.germlineTumorTypes = germlineTumorTypes;
        this.cancerSyndrome = cancerSyndrome;
        this.tissueType = tissueType;
        this.molecularGenetics = molecularGenetics;
        this.role = role;
        this.mutationTypes = mutationTypes;
        this.translocationPartner = translocationPartner;
        this.otherGermlineMut = otherGermlineMut;
        this.otherSyndrome = otherSyndrome;
    }

    @NotNull
    public String description() {
        return description;
    }

    @NotNull
    public String entrezId() {
        return entrezId;
    }

    @NotNull
    public String genomeLocation() {
        return genomeLocation;
    }

    @NotNull
    public String chromosomeBand() {
        return chromosomeBand;
    }

    @NotNull
    public String somatic() {
        return somatic;
    }

    @NotNull
    public String germline() {
        return germline;
    }

    @NotNull
    public String somaticTumorTypes() {
        return somaticTumorTypes;
    }

    @NotNull
    public String germlineTumorTypes() {
        return germlineTumorTypes;
    }

    @NotNull
    public String cancerSyndrome() {
        return cancerSyndrome;
    }

    @NotNull
    public String tissueType() {
        return tissueType;
    }

    @NotNull
    public String molecularGenetics() {
        return molecularGenetics;
    }

    @NotNull
    public String role() {
        return role;
    }

    @NotNull
    public String mutationTypes() {
        return mutationTypes;
    }

    @NotNull
    public String translocationPartner() {
        return translocationPartner;
    }

    @NotNull
    public String otherGermlineMut() {
        return otherGermlineMut;
    }

    @NotNull
    public String otherSyndrome() {
        return otherSyndrome;
    }
}
