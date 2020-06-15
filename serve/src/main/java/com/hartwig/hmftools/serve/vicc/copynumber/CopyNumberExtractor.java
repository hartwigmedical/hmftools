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

    private static final Set<String> ONCOKB_AMPLIFICATIONS = Sets.newHashSet("Amplification");

    private static final Set<String> ONCOKB_DELETIONS = Sets.newHashSet("Deletion");

    private static final Set<String> CGI_AMPLIFICATIONS = Sets.newHashSet("amplification");

    private static final Set<String> CGI_DELETIONS = Sets.newHashSet("deletion");

    private static final Set<String> CIVIC_AMPLIFICATIONS = Sets.newHashSet("AMPLIFICATION");

    private static final Set<String> CIVIC_DELETIONS = Sets.newHashSet("DELETION");

    @NotNull
    public Set<String> uniqueAmps() {
        return uniqueAmps;
    }

    @NotNull
    public Set<String> uniqueDels() {
        return uniqueDels;
    }

    @NotNull
    public Map<Feature, KnownAmplificationDeletion> extractKnownAmplificationsDeletions(@NotNull ViccEntry viccEntry) {
        Map<Feature, KnownAmplificationDeletion> ampsDelsPerFeature = Maps.newHashMap();
        boolean combinedEvent = false;

        if (viccEntry.source() == ViccSource.ONCOKB) {
            for (Feature feature : viccEntry.features()) {
                if (ONCOKB_AMPLIFICATIONS.contains(feature.name())) {
                    ampsDelsPerFeature.put(feature, eventForGene(feature.geneSymbol(), "amp", "OncoKB"));
                    uniqueAmps.add(feature.geneSymbol());
                } else if (ONCOKB_DELETIONS.contains(feature.name())) {
                    ampsDelsPerFeature.put(feature, eventForGene(feature.geneSymbol(), "del", "OncoKB"));
                    uniqueDels.add(feature.geneSymbol());
                }
            }
        } else if (viccEntry.source() == ViccSource.CIVIC) {
            for (Feature feature : viccEntry.features()) {
                String event = Strings.EMPTY;
                String featureName = feature.name();
                if (!featureName.contains("DEL") && !featureName.contains("Splicing alteration") && !featureName.contains("EXON")
                        && !featureName.contains("c.") && !featureName.contains("MUT") && !featureName.equals("LOSS-OF-FUNCTION")
                        && !featureName.equals("Gain-of-Function") && !featureName.contains("C.") && !featureName.equals(
                        "N-TERMINAL FRAME SHIFT") && !featureName.equals("COPY-NEUTRAL LOSS OF HETEROZYGOSITY")) {
                    //TODO
                } else if (featureName.contains("+") && !featureName.contains("c.") && !featureName.contains("C.")) {
                    //TODO

                } else {
                    ampsDelsPerFeature.put(feature, eventForGene(feature.geneSymbol(), "del", "CiViC"));
                }

                if (CIVIC_AMPLIFICATIONS.contains(event)) {
                    ampsDelsPerFeature.put(feature, eventForGene(feature.geneSymbol(), "amp", "CiViC"));
                    uniqueAmps.add(feature.geneSymbol());
                } else if (CIVIC_DELETIONS.contains(event)) {
                    ampsDelsPerFeature.put(feature, eventForGene(feature.geneSymbol(), "del", "CiViC"));
                    uniqueDels.add(feature.geneSymbol());
                }
            }
        } else if (viccEntry.source() == ViccSource.CGI) {
            for (Feature feature : viccEntry.features()) {
                String event = Strings.EMPTY;
                if (feature.name().contains("+")) {
                    LOGGER.info(feature);
                    String[] combinedEventConvertToSingleEvent = feature.name().split(" \\+ ", 2);
                    String gene = combinedEventConvertToSingleEvent[0].split(" ", 2)[0];

                    String geneCombined = combinedEventConvertToSingleEvent[1].split(" ", 2)[0];
                    String eventInfoCombined = combinedEventConvertToSingleEvent[1].split(" ", 2)[1];

                    //I assume, a combined event for actionability has 2 events. If more events, this will be not interpretated
                    if (combinedEventConvertToSingleEvent.length == 2) {
                        combinedEvent = true;

                        //TODO: fix combined events, to map

                        //                        if (eventMap.size() == 0) {
                        //                            eventMap.put(gene, Lists.newArrayList(eventInfo));
                        //                            if (eventMap.containsKey(geneCombined)) {
                        //                                eventMap.put(geneCombined, Lists.newArrayList(eventInfo, eventInfoCombined));
                        //                            } else {
                        //                                eventMap.put(gene, Lists.newArrayList(eventInfo));
                        //                                eventMap.put(geneCombined, Lists.newArrayList(eventInfoCombined));
                        //                            }
                        //                        }
                    }

                } else if (feature.name().split(" ", 2).length == 2) {
                    event = feature.name().split(" ")[1];
                } else {
                    if (feature.name().contains(":")) {
                        event = feature.name().split(":")[1];
                    }
                }

                if (CGI_AMPLIFICATIONS.contains(event)) {
                    ampsDelsPerFeature.put(feature, eventForGene(feature.geneSymbol(), "amp", "CGI"));
                    uniqueAmps.add(feature.geneSymbol());
                } else if (CGI_DELETIONS.contains(event)) {
                    ampsDelsPerFeature.put(feature, eventForGene(feature.geneSymbol(), "del", "CGI"));
                    uniqueDels.add(feature.geneSymbol());
                }
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
