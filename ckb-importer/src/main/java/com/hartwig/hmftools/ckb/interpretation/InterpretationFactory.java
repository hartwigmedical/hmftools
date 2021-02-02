package com.hartwig.hmftools.ckb.interpretation;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.ckb.datamodel.common.VariantInfo;
import com.hartwig.hmftools.ckb.datamodel.gene.Gene;
import com.hartwig.hmftools.ckb.datamodel.molecularprofile.MolecularProfile;
import com.hartwig.hmftools.ckb.datamodel.variant.Variant;

import org.jetbrains.annotations.NotNull;

public class InterpretationFactory {

    private InterpretationFactory() {

    }

    public static List<CkbEntryInterpretation> interpretationCkb(@NotNull CkbEntry ckbEntry) {
        List<CkbEntryInterpretation> CkbEntryInterpretation = Lists.newArrayList();

        for (MolecularProfile molecularProfile : ckbEntry.molecularProfile()) {
            ImmutableCkbEntryInterpretation.Builder outputBuilder = ImmutableCkbEntryInterpretation.builder();
            outputBuilder.molecularProfile(molecularProfile);

            for (VariantInfo variantInfo : molecularProfile.geneVariant()) {
                VariantInterpretation variantInterpretation = matchVariantInterpretation(ckbEntry, variantInfo.id());
                outputBuilder.variantInterpretation(variantInterpretation);

            }
            //            for (TreatmentApproachInfo treatmentApproachinfo : molecularProfile.treatmentApproach()) {
            //                List<TreatmentApproach> treatmentApproachMatching = matchTreatmentAppoach(ckbEntry, treatmentApproachinfo.id());
            //                CkbEntryInterpretation.add(outputBuilder.treatmentApproach(treatmentApproachMatching).build());
            //                for (TreatmentApproach treatmentApproach : treatmentApproachMatching) {
            //                    DrugClass drugClassMatching = null;
            //                    Therapy therapyMatching = null;
            //                    if (treatmentApproach.drugClass().id() != null) {
            //                        drugClassMatching = drugClassMatching(ckbEntry, treatmentApproach.drugClass().id());
            //                        for (DrugInfo drugInfo : drugClassMatching.drug()) {
            //                            List<Drug> drugMatching = drugMatching(ckbEntry, drugInfo.id());
            //                            CkbEntryInterpretation.add(outputBuilder.drug(drugMatching).build());
            //                        }
            //                    }
            //                    if (treatmentApproach.therapy().id() != null) {
            //                        therapyMatching = therapyMatching(ckbEntry, treatmentApproach.therapy().id());
            //                    }
            //
            //                    CkbEntryInterpretation.add(outputBuilder.drugClass(drugClassMatching).therapy(therapyMatching).build());
            //                }
            //            }

            CkbEntryInterpretation.add(outputBuilder.build());
        }
        return CkbEntryInterpretation;
    }

    @NotNull
    private static VariantInterpretation matchVariantInterpretation(@NotNull CkbEntry ckbEntry, int variantId) {
        ImmutableVariantInterpretation.Builder outputBuilder = ImmutableVariantInterpretation.builder();
        for (Variant variant : ckbEntry.variant()) {
            if (variant.id() == variantId) {
                int geneId = variant.gene().id();
                outputBuilder.variant(variant);
                for (Gene gene : ckbEntry.gene()) {
                    if (gene.id() == geneId) {
                        outputBuilder.gene(gene);
                    }
                }
            }
        }
        return outputBuilder.build();
    }
    //
    //
    //
    //    @NotNull
    //    private static List<TreatmentApproach> matchTreatmentAppoach(@NotNull CkbEntry ckbEntry, int treatmentApprochId) {
    //        List<TreatmentApproach> treatmentApproaches = Lists.newArrayList();
    //        for (TreatmentApproach treatmentApproach : ckbEntry.treatmentApproach()) {
    //            if (treatmentApproach.id() == treatmentApprochId) {
    //                treatmentApproaches.add(treatmentApproach);
    //            }
    //        }
    //        return treatmentApproaches;
    //    }
    //
    //    @Nullable
    //    private static DrugClass drugClassMatching(@NotNull CkbEntry ckbEntry, int drugClassId) {
    //        for (DrugClass drugClass : ckbEntry.drugClass()) {
    //            if (drugClass.id() == drugClassId) {
    //                return drugClass;
    //            }
    //        }
    //        return null;
    //    }
    //
    //    @Nullable
    //    private static Therapy therapyMatching(@NotNull CkbEntry ckbEntry, int therapyId) {
    //        for (Therapy therapy : ckbEntry.therapy()) {
    //            if (therapy.id() == therapyId) {
    //                return therapy;
    //            }
    //        }
    //        return null;
    //    }
    //
    //    @NotNull
    //    private static List<Drug> drugMatching(@NotNull CkbEntry ckbEntry, int drugId) {
    //        List<Drug> drugs = Lists.newArrayList();
    //        for (Drug drug : ckbEntry.drug()) {
    //            if (drug.id() == drugId) {
    //                drugs.add(drug);
    //            }
    //        }
    //        return drugs;
    //    }
}
