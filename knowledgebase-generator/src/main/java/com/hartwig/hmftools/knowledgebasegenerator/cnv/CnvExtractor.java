package com.hartwig.hmftools.knowledgebasegenerator.cnv;

import com.hartwig.hmftools.knowledgebasegenerator.actionability.gene.ActionableGene;
import com.hartwig.hmftools.knowledgebasegenerator.actionability.gene.ImmutableActionableGene;
import com.hartwig.hmftools.knowledgebasegenerator.sourceknowledgebase.Source;
import com.hartwig.hmftools.vicc.datamodel.KbSpecificObject;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
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
    public static KnownAmplificationDeletion determineKnownAmplificationDeletion(@NotNull Source source, @NotNull String typeEvent, @NotNull String gene){
        return knownInformation(source, typeEvent, gene);
    }

    @NotNull
    public static ActionableAmplificationDeletion determineActionableAmplificationDeletion(@NotNull Source source, @NotNull String typeEvent, @NotNull String gene){
        return actionableInformation(source, typeEvent, gene);
    }

    @NotNull
    private static KnownAmplificationDeletion knownInformation(@NotNull Source source, @NotNull String typeEvent, @NotNull String gene) {
        switch (source) {
            case ONCOKB:
                break;
            case CGI:
                break;
            case CIVIC:
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
            case MOLECULARMATCH:
                break;
            case MOLECULARMATCH_TRIALS:
                break;
            default:
                LOGGER.warn("Unknown knowledgebase");
        }
        return ImmutableKnownAmplificationDeletion.builder().gene(gene).eventType(typeEvent).source(source.toString()).build();
    }

    @NotNull
    private static ActionableAmplificationDeletion actionableInformation(@NotNull Source source, @NotNull String typeEvent, @NotNull String gene) {
        switch (source) {
            case ONCOKB:
                break;
            case CGI:
                break;
            case CIVIC:
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
            case MOLECULARMATCH:
                break;
            case MOLECULARMATCH_TRIALS:
                break;
            default:
                LOGGER.warn("Unknown knowledgebase");
        }
        return ImmutableActionableAmplificationDeletion.builder().gene(Strings.EMPTY).eventType(typeEvent).source(Strings.EMPTY).build();
    }

    @NotNull
    public static ActionableGene determineInfoOfEvent(@NotNull Source source, @NotNull String typeEvent,
            @NotNull KbSpecificObject kbSpecificObject, @NotNull String gene, @NotNull ViccEntry viccEntries) {
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
                level = Strings.EMPTY;
                direction = Strings.EMPTY;
                link = Strings.EMPTY;
                break;
            case CGI:
                Cgi kbCgi = (Cgi) kbSpecificObject;
                drug = kbCgi.drug();
                drugType = kbCgi.drugFamily();
                cancerType = kbCgi.primaryTumorType();
                level = viccEntries.association().evidenceLabel();
                direction = viccEntries.association().responseType();
                link = "https://www.cancergenomeinterpreter.org/biomarkers";
                break;
            case CIVIC:
                Civic kbCivic = (Civic) kbSpecificObject;
                drug = Strings.EMPTY;
                drugType = Strings.EMPTY;
                cancerType = Strings.EMPTY;
                level = Strings.EMPTY;
                direction = Strings.EMPTY;
                link = Strings.EMPTY;
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
            case MOLECULARMATCH:
                MolecularMatch kbMolecularMatch = (MolecularMatch) kbSpecificObject;
                break;
            case MOLECULARMATCH_TRIALS:
                MolecularMatchTrials kbMolecularMatchTrials = (MolecularMatchTrials) kbSpecificObject;
                break;
            default:
                LOGGER.warn("Unknown knowledgebase");
        }
        return ImmutableActionableGene.builder()
                .gene(gene)
                .type(typeEvent)
                .source(source.toString())
                .drug(drug)
                .drugType(drugType)
                .cancerType(cancerType)
                .level(level)
                .direction(direction)
                .link(link)
                .build();
    }
}
