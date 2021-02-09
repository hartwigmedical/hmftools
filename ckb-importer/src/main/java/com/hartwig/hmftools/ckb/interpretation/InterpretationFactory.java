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
import com.hartwig.hmftools.ckb.json.common.IndicationInfo;
import com.hartwig.hmftools.ckb.json.common.TherapyInfo;
import com.hartwig.hmftools.ckb.json.indication.Indication;
import com.hartwig.hmftools.ckb.json.molecularprofile.MolecularProfile;
import com.hartwig.hmftools.ckb.json.therapy.Therapy;
import com.hartwig.hmftools.ckb.util.DateConverter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class InterpretationFactory {

    private static final Logger LOGGER = LogManager.getLogger(InterpretationFactory.class);

    private InterpretationFactory() {
    }

    public static List<CkbEntryInterpretation> interpretationCkbDataModel(@NotNull CkbJsonDatabase ckbEntry) {
        List<CkbEntryInterpretation> CkbEntryInterpretation = Lists.newArrayList();

        for (MolecularProfile molecularProfile : ckbEntry.molecularProfile()) {
            ImmutableCkbEntryInterpretation.Builder outputBuilder = ImmutableCkbEntryInterpretation.builder();
            LOGGER.info("molecular profile {}", molecularProfile.id());
            for (ClinicalTrialInfo clinicalTrialInfo : molecularProfile.variantAssociatedClinicalTrial()) {
                for (ClinicalTrial clinicalTrial : ckbEntry.clinicalTrial()) {
                    if (clinicalTrialInfo.nctId().equals(clinicalTrial.nctId())) {
                        LOGGER.info("clinicalTrial {}", clinicalTrial.nctId());
                        outputBuilder.addClinicalTrial(ImmutableClinicalTrial.builder()
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
                                // TODO Use string comparison
                                if (Integer.parseInt(indicationInfo.id()) == indication.id()) {
                                    LOGGER.info(ImmutableIndication.builder()
                                            .id(indication.id())
                                            .name(indication.name())
                                            .source(indication.source())
                                            .definition(indication.definition())
                                            .currentPreferredTerm(indication.currentPreferredTerm())
                                            .lastUpdateDateFromDO(DateConverter.toDate(indication.lastUpdateDateFromDO()))
                                            .altIds(indication.altId())
                                            .termId(indication.termId())
                                            .build());
                                }
                            }

                        }

                        for (TherapyInfo therapyInfo : clinicalTrial.therapies()) {
                            for (Therapy therapy : ckbEntry.therapy()) {
                                if (therapyInfo.id() == therapy.id()) {
                                    LOGGER.info(ImmutableTherapy.builder()
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
            }
            LOGGER.info(outputBuilder.build());
            CkbEntryInterpretation.add(outputBuilder.build());
        }
        return CkbEntryInterpretation;
    }

}
