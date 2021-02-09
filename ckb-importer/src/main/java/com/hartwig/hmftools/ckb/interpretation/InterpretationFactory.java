package com.hartwig.hmftools.ckb.interpretation;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodel.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.datamodel.common.DescriptionInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ReferenceInfo;
import com.hartwig.hmftools.ckb.datamodel.common.VariantInfo;
import com.hartwig.hmftools.ckb.datamodel.gene.Gene;
import com.hartwig.hmftools.ckb.datamodel.molecularprofile.MolecularProfile;
import com.hartwig.hmftools.ckb.datamodel.reference.Reference;
import com.hartwig.hmftools.ckb.datamodel.variant.Variant;
import com.hartwig.hmftools.ckb.interpretation.varianttree.GeneInterpretation;
import com.hartwig.hmftools.ckb.interpretation.varianttree.ImmutableGeneInterpretation;
import com.hartwig.hmftools.ckb.interpretation.varianttree.ImmutableVariantInterpretation;
import com.hartwig.hmftools.ckb.interpretation.varianttree.ImmutableVariantTreeInterpretation;
import com.hartwig.hmftools.ckb.interpretation.varianttree.VariantInterpretation;
import com.hartwig.hmftools.ckb.interpretation.varianttree.VariantTreeInterpretation;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class InterpretationFactory {

    private InterpretationFactory() {

    }

    private static final Logger LOGGER = LogManager.getLogger(InterpretationFactory.class);

    public static List<CkbEntryInterpretation> interpretationCkbDataModel(@NotNull CkbJsonDatabase ckbEntry) {
        List<CkbEntryInterpretation> CkbEntryInterpretation = Lists.newArrayList();
        for (MolecularProfile molecularProfile : ckbEntry.molecularProfile()) {
            ImmutableCkbEntryInterpretation.Builder outputBuilder = ImmutableCkbEntryInterpretation.builder();
            outputBuilder.molecularProfile(molecularProfile);

            List<VariantTreeInterpretation> variantInterpretation = Lists.newArrayList();
            for (VariantInfo variantInfo : molecularProfile.geneVariant()) {
                variantInterpretation = matchVariantTreeInterpretation(ckbEntry, variantInfo.id(), variantInterpretation); //array
                outputBuilder.variantTreeInterpretations(variantInterpretation);

            }
            //
            //            for (TreatmentApproachInfo treatmentApproachInfo : molecularProfile.treatmentApproach()) {
            //                TreatmentInterpretation treatmentInterpretation =
            //                        matchTreatmentInterpretation(ckbEntry, treatmentApproachInfo.id()); //array
            //                outputBuilder.addTreatmentInterpretations(treatmentInterpretation);
            //            }
            //
            //            for (ClinicalTrialInfo clinicalTrialInfo: molecularProfile.variantAssociatedClinicalTrial()) {
            //                ClinicalTrialInterpretation evidenceInterpretation = matchClinicalTrialInterpretation(ckbEntry, clinicalTrialInfo.nctId()); //array
            //                outputBuilder.addClinicalTrialInterpretations(evidenceInterpretation);
            //            }
            //
            CkbEntryInterpretation.add(outputBuilder.build());
        }
        return CkbEntryInterpretation;
    }

    @NotNull
    public static VariantInterpretation extractVariantTree(@NotNull Variant variant, @NotNull CkbJsonDatabase ckbEntry) {
        ImmutableVariantInterpretation.Builder outputBuilder = ImmutableVariantInterpretation.builder();
        outputBuilder.variant(variant);
        for (DescriptionInfo descriptionInfo : variant.description()) {
            for (ReferenceInfo referenceInfo : descriptionInfo.reference()) {
                for (Reference reference : ckbEntry.reference()) {
                    if (reference.id() == referenceInfo.id()) {
                        outputBuilder.addReferences(reference);
                    }
                }
            }
        }
        return outputBuilder.build();
    }

    @NotNull
    public static GeneInterpretation extractGeneTree(@NotNull Gene gene, @NotNull CkbJsonDatabase ckbEntry) {
        ImmutableGeneInterpretation.Builder outputBuilder = ImmutableGeneInterpretation.builder();
        outputBuilder.gene(gene);
        for (DescriptionInfo descriptionInfo : gene.description()) {
            for (ReferenceInfo referenceInfo : descriptionInfo.reference()) {
                for (Reference reference : ckbEntry.reference()) {
                    if (reference.id() == referenceInfo.id()) {
                        outputBuilder.addReferences(reference);
                    }
                }
            }
        }
        return outputBuilder.build();
    }

    @NotNull
    private static List<VariantTreeInterpretation> matchVariantTreeInterpretation(@NotNull CkbJsonDatabase ckbEntry, int variantId,
            @NotNull List<VariantTreeInterpretation> variantInterpretationList) {
        ImmutableVariantTreeInterpretation.Builder outputBuilder = ImmutableVariantTreeInterpretation.builder();
        for (Variant variant : ckbEntry.variant()) {
            if (variant.id() == variantId) {
                int geneId = variant.gene().id();
                VariantInterpretation variantInterpretation = extractVariantTree(variant, ckbEntry);
                outputBuilder.variantInterpretation(variantInterpretation);

                for (Gene gene : ckbEntry.gene()) {
                    if (gene.id() == geneId) {
                        GeneInterpretation geneInterpretation = extractGeneTree(gene, ckbEntry);
                        outputBuilder.geneInterpretation(geneInterpretation);
                    }
                }
            }
            variantInterpretationList.add(outputBuilder.build());
        }
        return variantInterpretationList;
    }

    //    @NotNull
    //    private static ClinicalTrialInterpretation matchClinicalTrialInterpretation(@NotNull CkbEntry ckbEntry, @NotNull String nctId) {
    //        ImmutableClinicalTrialInterpretation.Builder outputBuilder = ImmutableClinicalTrialInterpretation.builder();
    //
    //        for (ClinicalTrial clinicalTrial: ckbEntry.clinicalTrial()) {
    //            if (clinicalTrial.nctId().equals(nctId)) {
    //                outputBuilder.clinicalTrial(clinicalTrial);
    //                for (IndicationInfo indicationInfo: clinicalTrial.indication()) {
    //                    for (Indication indication: ckbEntry.indication()) {
    //                        if (indicationInfo.id() == indication.id()) {
    //                            outputBuilder.addIndications(indication);
    //                        }
    //                    }
    //                }
    //
    //            }
    //        }
    //        return outputBuilder.build();
    //    }
    //

    //
    //    @NotNull
    //    private static TreatmentInterpretation matchTreatmentInterpretation(@NotNull CkbEntry ckbEntry, int treatmentApprochId) {
    //        ImmutableTreatmentInterpretation.Builder outputBuilder = ImmutableTreatmentInterpretation.builder();
    //        for (TreatmentApproach treatmentApproach : ckbEntry.treatmentApproach()) {
    //            if (treatmentApproach.id() == treatmentApprochId) {
    //                outputBuilder.treatmentApproach(treatmentApproach);
    //                outputBuilder.treatmentApproachInterpretation(matchTreatmentApprochInterpretation(ckbEntry,
    //                        treatmentApproach)); //array
    //            }
    //        }
    //        return outputBuilder.build();
    //    }
    //
    //    @NotNull
    //    private static TreatmentApprochInterpretation matchTreatmentApprochInterpretation(@NotNull CkbEntry ckbEntry,
    //            @NotNull TreatmentApproach treatmentApproach) {
    //        ImmutableTreatmentApprochInterpretation.Builder outputBuilder = ImmutableTreatmentApprochInterpretation.builder();
    //
    //        for (DrugClass drugClass : ckbEntry.drugClass()) {
    //            if (treatmentApproach.drugClass() != null) {
    //                outputBuilder.drugClassInterpretation(matchDrugClassInterpretation(drugClass, treatmentApproach, ckbEntry));
    //            }
    //        }
    //
    //        for (Therapy therapy : ckbEntry.therapy()) {
    //            if (treatmentApproach.therapy() != null) {
    //                if (therapy.id() == treatmentApproach.therapy().id()) {
    //                    outputBuilder.therapy(therapy); //object
    //                }
    //            }
    //        }
    //        return outputBuilder.build();
    //    }
    //
    //    @NotNull
    //    private static DrugClassInterpretation matchDrugClassInterpretation(@NotNull DrugClass drugClass,
    //            @NotNull TreatmentApproach treatmentApproach, @NotNull CkbEntry ckbEntry) {
    //        ImmutableDrugClassInterpretation.Builder outputBuilder = ImmutableDrugClassInterpretation.builder();
    //        if (drugClass.id() == treatmentApproach.drugClass().id()) {
    //            outputBuilder.addDrugClasses(drugClass);  //object
    //            for (DrugInfo drugInfo : drugClass.drug()) { //array
    //                for (Drug drug : ckbEntry.drug()) {
    //                    if (drugInfo.id() == drug.id()) {
    //                        outputBuilder.drug(drug);
    //                    }
    //                }
    //            }
    //        }
    //        return outputBuilder.build();
    //    }

}
