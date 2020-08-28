package com.hartwig.hmftools.serve.vicc.copynumber;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.KbSpecificObject;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;
import com.hartwig.hmftools.vicc.datamodel.brca.Brca;
import com.hartwig.hmftools.vicc.datamodel.cgi.Cgi;
import com.hartwig.hmftools.vicc.datamodel.civic.Civic;
import com.hartwig.hmftools.vicc.datamodel.jax.Jax;
import com.hartwig.hmftools.vicc.datamodel.jaxtrials.JaxTrials;
import com.hartwig.hmftools.vicc.datamodel.molecularmatch.MolecularMatch;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.MolecularMatchTrials;
import com.hartwig.hmftools.vicc.datamodel.oncokb.OncoKb;
import com.hartwig.hmftools.vicc.datamodel.pmkb.Pmkb;
import com.hartwig.hmftools.vicc.datamodel.sage.Sage;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class CopyNumberExtractor {

    private static final Logger LOGGER = LogManager.getLogger(CopyNumberExtractor.class);

    @NotNull
    private final Set<String> uniqueAmps = Sets.newHashSet();
    @NotNull
    private final Set<String> uniqueDels = Sets.newHashSet();

    private static final Set<String> AMPLIFICATIONS = Sets.newHashSet("Amplification",
            "amplification",
            "AMPLIFICATION",
            "amp",
            "overexpression",
            "over exp",
            "amp over exp",
            "OVEREXPRESSION",
            "Overexpression");

    private static final Set<String> DELETIONS = Sets.newHashSet("Deletion",
            "deletion",
            "DELETION",
            "del",
            "undexpression",
            "dec exp",
            "UNDEREXPRESSION",
            "loss",
            "LOSS", "Copy Number Loss");

    @NotNull
    public Set<String> uniqueAmps() {
        return uniqueAmps;
    }

    @NotNull
    public Set<String> uniqueDels() {
        return uniqueDels;
    }

    private boolean isAmplification(@NotNull Feature feature) {
        String eventKeyAmplification = extractKeyAmplificationDeletion(feature);
        if (eventKeyAmplification.equals("Amplification")) {
            return true;
        } else {
            return false;
        }
    }

    private boolean isDeletion(@NotNull Feature feature) {
        String eventKeyDeletion = extractKeyAmplificationDeletion(feature);

        if (eventKeyDeletion.equals("Deletion")) {
            return true;
        } else {
            return false;
        }
    }

    private String extractKeyAmplificationDeletion(@NotNull Feature feature) {
        //TODO: fix combi events
        String featureName = feature.name();
        if (featureName.contains(" ") && !featureName.equals("Copy Number Loss")) {
            featureName = featureName.split(" ", 2)[1];
        }

        if (AMPLIFICATIONS.contains(featureName) || AMPLIFICATIONS.contains(feature.biomarkerType())) {
            return "Amplification";
        } else if (DELETIONS.contains(featureName) || DELETIONS.contains(feature.biomarkerType())) {
            return "Deletion";
        } else {
            return Strings.EMPTY;
        }
    }

    @NotNull
    public Map<Feature, KnownAmplificationDeletion> extractKnownAmplificationsDeletions(@NotNull ViccEntry viccEntry) {
        Map<Feature, KnownAmplificationDeletion> ampsDelsPerFeature = Maps.newHashMap();

        for (Feature feature : viccEntry.features()) {
            if (isAmplification(feature)) {
                ampsDelsPerFeature.put(feature, eventForGene(feature.geneSymbol(), "amp", viccEntry.source().display()));
                uniqueAmps.add(feature.geneSymbol());
            } else if (isDeletion(feature)) {
                ampsDelsPerFeature.put(feature, eventForGene(feature.geneSymbol(), "del", viccEntry.source().display()));
                uniqueDels.add(feature.geneSymbol());
            }
        }
        return ampsDelsPerFeature;
    }

    @NotNull
    private static KnownAmplificationDeletion eventForGene(@NotNull String gene, @NotNull String eventType, @NotNull String database) {
        return ImmutableKnownAmplificationDeletion.builder().gene(gene).source(database).eventType(eventType).sourceLink("link").build();
    }

    @NotNull
    public static KnownAmplificationDeletion determineKnownAmplificationDeletion(@NotNull ViccSource source, @NotNull String typeEvent,
            @NotNull String gene) {
        return knownInformation(source, typeEvent, gene);

    }

    @NotNull
    public static ActionableAmplificationDeletion determineActionableAmplificationDeletion(@NotNull ViccSource source,
            @NotNull String typeEvent, @NotNull String gene, @NotNull ViccEntry viccEntry) {
        return actionableInformation(source, typeEvent, gene, viccEntry);
    }

    @NotNull
    private static KnownAmplificationDeletion knownInformation(@NotNull ViccSource source, @NotNull String typeEvent,
            @NotNull String gene) {
        String link = Strings.EMPTY;
        switch (source) {
            case ONCOKB:
                link = "link_oncokb";
                break;
            case CGI:
                link = "link_cgi";
                break;
            case CIVIC:
                link = "link_civic";
                break;
            case JAX:
                break;
            case JAX_TRIALS:
                break;
            case BRCA:
                break;
            case SAGE:
                break;
            case PMKB:
                break;
            case MOLECULAR_MATCH:
                break;
            case MOLECULAR_MATCH_TRIALS:
                break;
            default:
                LOGGER.warn("Unknown knowledgebase");
        }
        return ImmutableKnownAmplificationDeletion.builder()
                .gene(gene)
                .eventType(typeEvent)
                .source(source.toString())
                .sourceLink(link)
                .build();
    }

    @NotNull
    private static ActionableAmplificationDeletion actionableInformation(@NotNull ViccSource source, @NotNull String typeEvent,
            @NotNull String gene, @NotNull ViccEntry viccEntry) {
        KbSpecificObject kbSpecificObject = viccEntry.kbSpecificObject();
        String drug = Strings.EMPTY;
        String drugType = Strings.EMPTY;
        String cancerType = Strings.EMPTY;
        String level = Strings.EMPTY;
        String direction = Strings.EMPTY;
        String link = Strings.EMPTY;
        switch (source) {
            case ONCOKB:
                OncoKb kbOncoKb = (OncoKb) kbSpecificObject;
                drug = Strings.EMPTY;
                drugType = Strings.EMPTY;
                cancerType = Strings.EMPTY;
                level = viccEntry.association().evidenceLabel();
                direction = viccEntry.association().responseType();
                link = "http://oncokb.org/#/gene/" + gene + "/alteration/" + "[variantName]";

                break;
            case CGI: // all events are actionable
                Cgi kbCgi = (Cgi) kbSpecificObject;
                drug = kbCgi.drug();
                drugType = kbCgi.drugFamily();
                cancerType = kbCgi.primaryTumorType();
                level = viccEntry.association().evidenceLabel();
                direction = viccEntry.association().responseType();
                link = "https://www.cancergenomeinterpreter.org/biomarkers";
                break;
            case CIVIC:
                Civic kbCivic = (Civic) kbSpecificObject;
                drug = Strings.EMPTY;
                drugType = Strings.EMPTY;
                cancerType = kbCivic.evidenceItem().disease().name();
                level = viccEntry.association().evidenceLabel();
                direction = viccEntry.association().responseType();
                link = "https://civic.genome.wustl.edu/links/variants/" + kbCivic.id();
                break;
            case JAX:
                Jax kbJax = (Jax) kbSpecificObject;
                break;
            case JAX_TRIALS:
                JaxTrials kbJaxTrials = (JaxTrials) kbSpecificObject;
                break;
            case BRCA:
                Brca kbBrca = (Brca) kbSpecificObject;
                break;
            case SAGE:
                Sage kbSage = (Sage) kbSpecificObject;
                break;
            case PMKB:
                Pmkb kbPmkb = (Pmkb) kbSpecificObject;
                break;
            case MOLECULAR_MATCH:
                MolecularMatch kbMolecularMatch = (MolecularMatch) kbSpecificObject;
                break;
            case MOLECULAR_MATCH_TRIALS:
                MolecularMatchTrials kbMolecularMatchTrials = (MolecularMatchTrials) kbSpecificObject;
                break;
            default:
                LOGGER.warn("Unknown knowledgebase");
        }
        return ImmutableActionableAmplificationDeletion.builder()
                .gene(gene)
                .eventType(typeEvent)
                .source(source.toString())
                .drug(drug)
                .drugType(drugType)
                .cancerType(cancerType)
                .level(level)
                .direction(direction)
                .sourceLink(link)
                .build();

    }
}
