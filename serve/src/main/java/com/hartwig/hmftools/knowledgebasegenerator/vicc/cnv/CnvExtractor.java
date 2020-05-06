package com.hartwig.hmftools.knowledgebasegenerator.vicc.cnv;

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

public class CnvExtractor {

    private static final Logger LOGGER = LogManager.getLogger(CnvExtractor.class);

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
        KbSpecificObject kbSpecificObject = viccEntry.KbSpecificObject();
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
