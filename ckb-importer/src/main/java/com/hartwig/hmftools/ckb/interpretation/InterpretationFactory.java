package com.hartwig.hmftools.ckb.interpretation;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.ckb.datamodel.common.TreatmentApproachInfo;
import com.hartwig.hmftools.ckb.datamodel.common.VariantInfo;
import com.hartwig.hmftools.ckb.datamodel.drug.Drug;
import com.hartwig.hmftools.ckb.datamodel.drugclass.DrugClass;
import com.hartwig.hmftools.ckb.datamodel.gene.Gene;
import com.hartwig.hmftools.ckb.datamodel.molecularprofile.MolecularProfile;
import com.hartwig.hmftools.ckb.datamodel.therapy.Therapy;
import com.hartwig.hmftools.ckb.datamodel.treatmentapproach.TreatmentApproach;
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
                outputBuilder.addVariantInterpretation(variantInterpretation);

            }

//            for (TreatmentApproachInfo treatmentApproachInfo : molecularProfile.treatmentApproach()) {
//                TreatmentInterpretation treatmentInterpretation = matchTreatmentInterpretation(ckbEntry, treatmentApproachInfo.id());
//                outputBuilder.treatmentInterpretation(treatmentInterpretation);
//
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
                outputBuilder.variant(variant); //array
                for (Gene gene : ckbEntry.gene()) {
                    if (gene.id() == geneId) {
                        outputBuilder.gene(gene);  //object
                    }
                }
            }
        }
        return outputBuilder.build();
    }

    @NotNull
    private static TreatmentInterpretation matchTreatmentInterpretation(@NotNull CkbEntry ckbEntry, int treatmentApprochId) {
        ImmutableTreatmentInterpretation.Builder outputBuilder = ImmutableTreatmentInterpretation.builder();
        for (TreatmentApproach treatmentApproach : ckbEntry.treatmentApproach()) {
            if (treatmentApproach.id() == treatmentApprochId) {
                outputBuilder.treatmentApproach(treatmentApproach);  //Array

                for (DrugClass drugClass: ckbEntry.drugClass()) {
                    if (treatmentApproach.drugClass() != null){
                        if (drugClass.id() == treatmentApproach.drugClass().id()) {
                            outputBuilder.drugClass(drugClass);  //object
                        }
                    } else {
                        outputBuilder.drugClass(null);
                    }
                }

                for (Therapy therapy: ckbEntry.therapy()) {
                    if (treatmentApproach.therapy() != null) {
                        if (therapy.id() == treatmentApproach.therapy().id()) {
                            outputBuilder.therapy(therapy); //object
                        }
                    } else {
                        outputBuilder.therapy(null);
                    }
                }
            }
        }
        return outputBuilder.build();
    }

}
