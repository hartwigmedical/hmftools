package com.hartwig.hmftools.ckb.interpretation;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodelinterpretation.CkbEntryInterpretation;
import com.hartwig.hmftools.ckb.datamodelinterpretation.ImmutableCkbEntryInterpretation;
import com.hartwig.hmftools.ckb.datamodelinterpretation.clinicaltrial.ImmutableClinicalTrial;
import com.hartwig.hmftools.ckb.datamodelinterpretation.indication.ImmutableIndication;
import com.hartwig.hmftools.ckb.datamodelinterpretation.therapy.ImmutableTherapy;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.clinicaltrial.ClinicalTrial;
import com.hartwig.hmftools.ckb.json.common.ClinicalTrialInfo;
import com.hartwig.hmftools.ckb.json.common.EvidenceInfo;
import com.hartwig.hmftools.ckb.json.common.IndicationInfo;
import com.hartwig.hmftools.ckb.json.common.TherapyInfo;
import com.hartwig.hmftools.ckb.json.indication.Indication;
import com.hartwig.hmftools.ckb.json.molecularprofile.MolecularProfile;
import com.hartwig.hmftools.ckb.json.therapy.Therapy;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class InterpretationFactory {

    private static final Logger LOGGER = LogManager.getLogger(InterpretationFactory.class);

    private InterpretationFactory() {
    }

    public static List<CkbEntryInterpretation> interpretationCkbDataModel(@NotNull CkbJsonDatabase ckbEntry) {
        List<CkbEntryInterpretation> CkbEntryInterpretation = Lists.newArrayList();
        int ckbId = 0;
        for (MolecularProfile molecularProfile : ckbEntry.molecularProfile()) {
            ++ckbId;
            ImmutableCkbEntryInterpretation.Builder outputBuilder = ImmutableCkbEntryInterpretation.builder();
            outputBuilder.id(ckbId);
            LOGGER.info("molecular profile {}", molecularProfile.id());

            //extract clinical trial information
            for (ClinicalTrialInfo clinicalTrialInfo : molecularProfile.variantAssociatedClinicalTrial()) {
                ImmutableClinicalTrialInterpretation.Builder outputBuilderClinicalInterpretation = ImmutableClinicalTrialInterpretation.builder();

                for (ClinicalTrial clinicalTrial : ckbEntry.clinicalTrial()) {
                    if (clinicalTrialInfo.nctId().equals(clinicalTrial.nctId())) {
                        LOGGER.info("clinicalTrial {}", clinicalTrial.nctId());
                        outputBuilderClinicalInterpretation.clinicalTrials(ImmutableClinicalTrial.builder()
                                .nctId(clinicalTrial.nctId())
                                .title(clinicalTrial.title())
                                .phase(clinicalTrial.phase())
                                .recruitment(clinicalTrial.recruitment())
                                .ageGroups(clinicalTrial.ageGroups())
                                .gender(clinicalTrial.gender())
                                .variantRequirement(clinicalTrial.variantRequirements())
                                .sponsor(clinicalTrial.sponsors())
                                .updateDate(clinicalTrial.updateDate())
                                .clinicalTrialVariantRequirementDetails(Lists.newArrayList())
                                .locations(Lists.newArrayList())
                                .build());

                        for (IndicationInfo indicationInfo : clinicalTrial.indications()) {
                            for (Indication indication : ckbEntry.indication()) {
                                if (indicationInfo.id().equals(indication.id())) {
                                    outputBuilderClinicalInterpretation.addIndications(ImmutableIndication.builder()
                                            // TODO Switch to String for ID
                                            .id(Integer.parseInt(indication.id()))
                                            .name(indication.name())
                                            .source(indication.source())
                                            .definition(indication.definition())
                                            .currentPreferredTerm(indication.currentPreferredTerm())
                                            .lastUpdateDateFromDO(indication.lastUpdateDateFromDO())
                                            .altIds(indication.altIds())
                                            .termId(indication.termId())
                                            .build());
                                }
                            }
                        }

                        for (TherapyInfo therapyInfo : clinicalTrial.therapies()) {
                            for (Therapy therapy : ckbEntry.therapy()) {
                                if (therapyInfo.id() == therapy.id()) {
                                    outputBuilderClinicalInterpretation.addTherapies(ImmutableTherapy.builder()
                                            .id(therapy.id())
                                            .therapyName(therapy.therapyName())
                                            .synonyms(therapy.synonyms())
                                            .descriptions(Lists.newArrayList())
                                            .createDate(therapy.createDate())
                                            .updateDate(therapy.updateDate())
                                            .build());
                                }
                            }
                        }
                    }
                }
                outputBuilder.addClinicalTrialInterpretation(outputBuilderClinicalInterpretation.build());
            }

            //extract variant level information
            for (EvidenceInfo evidenceInfo : molecularProfile.variantLevelEvidence().evidence()) {
                outputBuilder.addEvidenceInterpretation();
            }
            LOGGER.info(outputBuilder.build());
            CkbEntryInterpretation.add(outputBuilder.build());
        }
        return CkbEntryInterpretation;
    }
}