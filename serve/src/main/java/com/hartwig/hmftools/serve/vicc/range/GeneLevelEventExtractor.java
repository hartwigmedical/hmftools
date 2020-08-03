package com.hartwig.hmftools.serve.vicc.range;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class GeneLevelEventExtractor {
    private static final Logger LOGGER = LogManager.getLogger(GeneLevelEventExtractor.class);

    private static final Set<String> GENE = Sets.newHashSet("mut",
            "mutant",
            "expression",
            "gene_only",
            "EXPRESSION",
            "Protein Altering Variant",
            "Frameshift Truncation",
            "Gene Variant",
            "Transcript Variant",
            "Inframe Insertion",
            "Wild Type",
            "Transcription Variant",
            "Frameshift Variant",
            "Coding Transcript Intron Variant",
            "Nonsynonymous Variant",
            "Inframe Variant",
            "Inframe Indel", "MUTATION", "Frameshift Elongation", "wild-type");

    private static final Set<String> GENE_ACTIVATION = Sets.newHashSet("Gain-of-function Mutations",
            "act mut",
            "overexpression",
            "over exp",
            "pos",
            "positive",
            "amp over exp",
            "OVEREXPRESSION",
            "Gain Of Function Variant",
            "Stop Lost", "Missense Variant");

    private static final Set<String> GENE_INACTIVATION = Sets.newHashSet("Truncating Mutations",
            "inact mut",
            "loss",
            "biallelic inactivation",
            "undexpression",
            "dec exp",
            "negative",
            "is_deletion",
            "UNDEREXPRESSION",
            "Loss Of Function Variant",
            "Start Lost", "Loss Of Heterozygosity");

    @NotNull
    public Map<Feature, String> extractKnownGeneLevelEvents(@NotNull ViccEntry viccEntry) {
        Map<Feature, String> geneLevelEventsPerFeature = Maps.newHashMap();

        for (Feature feature : viccEntry.features()) {
            if (GENE_ACTIVATION.contains(feature.name()) || GENE_ACTIVATION.contains(feature.biomarkerType()) || GENE_ACTIVATION.contains(
                    feature.proteinAnnotation())) {
                geneLevelEventsPerFeature.put(feature, "gain of " + feature.geneSymbol());
            } else if (GENE_INACTIVATION.contains(feature.name()) || GENE_INACTIVATION.contains(feature.biomarkerType())
                    || GENE_INACTIVATION.contains(feature.provenanceRule()) || GENE_INACTIVATION.contains(feature.proteinAnnotation())) {
                geneLevelEventsPerFeature.put(feature, "loss of " + feature.geneSymbol());
            }
//            else if (GENE.contains(feature.biomarkerType()) || GENE.contains(feature.provenanceRule())
//                    || (GENE.contains(feature.proteinAnnotation())) && !feature.name().contains("+")
//                    || GENE.contains(feature.name())) { //TODO: determine gain of loss function
//                geneLevelEventsPerFeature.put(feature, "gain/loss of " + feature.geneSymbol());
//            }
        }

        return geneLevelEventsPerFeature;
    }
}
