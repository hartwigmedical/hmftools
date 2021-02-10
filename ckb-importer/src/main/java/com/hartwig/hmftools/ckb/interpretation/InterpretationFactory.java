package com.hartwig.hmftools.ckb.interpretation;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodelinterpretation.CkbEntryInterpretation;
import com.hartwig.hmftools.ckb.datamodelinterpretation.ImmutableCkbEntryInterpretation;
import com.hartwig.hmftools.ckb.datamodelinterpretation.clinicaltrial.ClinicalTrialLocation;
import com.hartwig.hmftools.ckb.datamodelinterpretation.clinicaltrial.ClinicalTrialVariantRequirementDetail;
import com.hartwig.hmftools.ckb.datamodelinterpretation.clinicaltrial.ImmutableClinicalTrial;
import com.hartwig.hmftools.ckb.datamodelinterpretation.clinicaltrial.ImmutableClinicalTrialContact;
import com.hartwig.hmftools.ckb.datamodelinterpretation.clinicaltrial.ImmutableClinicalTrialLocation;
import com.hartwig.hmftools.ckb.datamodelinterpretation.clinicaltrial.ImmutableClinicalTrialVariantRequirementDetail;
import com.hartwig.hmftools.ckb.datamodelinterpretation.common.ImmutableReferenceExtend;
import com.hartwig.hmftools.ckb.datamodelinterpretation.common.ReferenceExtend;
import com.hartwig.hmftools.ckb.datamodelinterpretation.drug.DrugDescription;
import com.hartwig.hmftools.ckb.datamodelinterpretation.drug.ImmutableDrug;
import com.hartwig.hmftools.ckb.datamodelinterpretation.drug.ImmutableDrugDescription;
import com.hartwig.hmftools.ckb.datamodelinterpretation.drugclass.ImmutableDrugClass;
import com.hartwig.hmftools.ckb.datamodelinterpretation.gene.GeneDescription;
import com.hartwig.hmftools.ckb.datamodelinterpretation.gene.ImmutableGene;
import com.hartwig.hmftools.ckb.datamodelinterpretation.gene.ImmutableGeneDescription;
import com.hartwig.hmftools.ckb.datamodelinterpretation.globaltherapyapprovalstatus.GlobalTherapyApprovalStatus;
import com.hartwig.hmftools.ckb.datamodelinterpretation.globaltherapyapprovalstatus.ImmutableGlobalTherapyApprovalStatus;
import com.hartwig.hmftools.ckb.datamodelinterpretation.indication.ImmutableIndication;
import com.hartwig.hmftools.ckb.datamodelinterpretation.therapy.ImmutableTherapy;
import com.hartwig.hmftools.ckb.datamodelinterpretation.therapy.ImmutableTherapyDescription;
import com.hartwig.hmftools.ckb.datamodelinterpretation.therapy.TherapyDescription;
import com.hartwig.hmftools.ckb.datamodelinterpretation.variant.CategoryVariantPath;
import com.hartwig.hmftools.ckb.datamodelinterpretation.variant.ImmutableCategoryVariantPath;
import com.hartwig.hmftools.ckb.datamodelinterpretation.variant.ImmutableMemberVariant;
import com.hartwig.hmftools.ckb.datamodelinterpretation.variant.ImmutableReferenceTranscriptCoordinate;
import com.hartwig.hmftools.ckb.datamodelinterpretation.variant.ImmutableVariant;
import com.hartwig.hmftools.ckb.datamodelinterpretation.variant.ImmutableVariantDescription;
import com.hartwig.hmftools.ckb.datamodelinterpretation.variant.ImmutableVariantInfo;
import com.hartwig.hmftools.ckb.datamodelinterpretation.variant.MemberVariant;
import com.hartwig.hmftools.ckb.datamodelinterpretation.variant.ReferenceTranscriptCoordinate;
import com.hartwig.hmftools.ckb.datamodelinterpretation.variant.VariantDescription;
import com.hartwig.hmftools.ckb.interpretation.clinicaltrial.ImmutableClinicalTrialInterpretation;
import com.hartwig.hmftools.ckb.interpretation.common.DrugsInterpretation;
import com.hartwig.hmftools.ckb.interpretation.common.ImmutableDrugsInterpretation;
import com.hartwig.hmftools.ckb.interpretation.common.ImmutableTherapyInterpretation;
import com.hartwig.hmftools.ckb.interpretation.common.ImmutableVariantInterpretation;
import com.hartwig.hmftools.ckb.interpretation.common.TherapyInterpretation;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.clinicaltrial.ClinicalTrial;
import com.hartwig.hmftools.ckb.json.clinicaltrial.ClinicalTrialContact;
import com.hartwig.hmftools.ckb.json.common.ClinicalTrialInfo;
import com.hartwig.hmftools.ckb.json.common.DescriptionInfo;
import com.hartwig.hmftools.ckb.json.common.DrugClassInfo;
import com.hartwig.hmftools.ckb.json.common.DrugInfo;
import com.hartwig.hmftools.ckb.json.common.EvidenceInfo;
import com.hartwig.hmftools.ckb.json.common.GlobalApprovalStatusInfo;
import com.hartwig.hmftools.ckb.json.common.IndicationInfo;
import com.hartwig.hmftools.ckb.json.common.ReferenceInfo;
import com.hartwig.hmftools.ckb.json.common.TherapyInfo;
import com.hartwig.hmftools.ckb.json.common.VariantInfo;
import com.hartwig.hmftools.ckb.json.drug.Drug;
import com.hartwig.hmftools.ckb.json.drugclass.DrugClass;
import com.hartwig.hmftools.ckb.json.gene.Gene;
import com.hartwig.hmftools.ckb.json.indication.Indication;
import com.hartwig.hmftools.ckb.json.molecularprofile.MolecularProfile;
import com.hartwig.hmftools.ckb.json.reference.Reference;
import com.hartwig.hmftools.ckb.json.therapy.Therapy;
import com.hartwig.hmftools.ckb.json.variant.Variant;
import com.hartwig.hmftools.ckb.json.variant.VariantCategoryVariantPath;
import com.hartwig.hmftools.ckb.json.variant.VariantTranscriptCoordinate;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class InterpretationFactory {

    private static final Logger LOGGER = LogManager.getLogger(InterpretationFactory.class);

    private InterpretationFactory() {
    }

    public static List<CkbEntryInterpretation> interpretationCkbDataModel(@NotNull CkbJsonDatabase ckbEntry) {
        List<CkbEntryInterpretation> CkbEntryInterpretation = Lists.newArrayList();
        int ckbId = 0;
        for (MolecularProfile molecularProfile : ckbEntry.molecularProfiles()) {
            LOGGER.info(molecularProfile.id());
            ++ckbId;
            ImmutableCkbEntryInterpretation.Builder outputBuilder = ImmutableCkbEntryInterpretation.builder();
            outputBuilder.id(ckbId);

            //extract clinical trial information
            for (ClinicalTrialInfo clinicalTrialInfo : molecularProfile.variantAssociatedClinicalTrials()) {
                ImmutableClinicalTrialInterpretation.Builder outputBuilderClinicalInterpretation =
                        ImmutableClinicalTrialInterpretation.builder();

                for (ClinicalTrial clinicalTrial : ckbEntry.clinicalTrials()) {
                    if (clinicalTrialInfo.nctId().equals(clinicalTrial.nctId())) {
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
                                .clinicalTrialVariantRequirementDetails(extractProfileName(clinicalTrial.variantRequirementDetails(),
                                        molecularProfile,
                                        ckbEntry))
                                .locations(extractClinicalTrialLocation(clinicalTrial.clinicalTrialLocations()))
                                .build());

                        for (IndicationInfo indicationInfo : clinicalTrial.indications()) {
                            for (Indication indication : ckbEntry.indications()) {
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
                            for (Therapy therapy : ckbEntry.therapies()) {
                                if (therapyInfo.id() == therapy.id()) {

                                    outputBuilderClinicalInterpretation.addTherapyInterpretations(extractTherapyInterpretation(therapy,
                                            ckbEntry));
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

    @NotNull
    private static TherapyInterpretation extractTherapyInterpretation(@NotNull Therapy therapy, @NotNull CkbJsonDatabase ckbEntry) {
        return ImmutableTherapyInterpretation.builder()
                .therapy(extractTherapy(therapy, ckbEntry))
                .addDrugs(extractDrugsInterpretation(therapy.drugs(), ckbEntry))
                .globalTherapyApprovalStatuses(extractGlobalApprovalStatus(therapy.globalApprovalStatuses()))
                .build();

    }

    @NotNull
    private static List<GlobalTherapyApprovalStatus> extractGlobalApprovalStatus(
            @NotNull List<GlobalApprovalStatusInfo> globalTherapyApprovalStatuses) {
        List<GlobalTherapyApprovalStatus> globalTherapyApprovalStatusesInterpretation = Lists.newArrayList();
        for (GlobalApprovalStatusInfo globalTherapyApprovalStatusInfo : globalTherapyApprovalStatuses) {
            globalTherapyApprovalStatusesInterpretation.add(ImmutableGlobalTherapyApprovalStatus.builder()
                    .id(globalTherapyApprovalStatusInfo.id())
                    .approvalStatus(globalTherapyApprovalStatusInfo.approvalStatus())
                    .approvalAuthority(globalTherapyApprovalStatusInfo.approvalAuthority())
                    .build());

        } return globalTherapyApprovalStatusesInterpretation;
    }

    @NotNull
    private static com.hartwig.hmftools.ckb.datamodelinterpretation.therapy.Therapy extractTherapy(@NotNull Therapy therapy,
            @NotNull CkbJsonDatabase ckbEntry) {
        return ImmutableTherapy.builder()
                .id(therapy.id())
                .therapyName(therapy.therapyName())
                .synonyms(therapy.synonyms())
                .descriptions(extractTherapyDescriptions(therapy.descriptions(), ckbEntry))
                .createDate(therapy.createDate())
                .updateDate(therapy.updateDate())
                .build();

    }

    @NotNull
    private static DrugsInterpretation extractDrugsInterpretation(@NotNull List<DrugInfo> drugs, @NotNull CkbJsonDatabase ckbEntry) {
        ImmutableDrugsInterpretation.Builder outputBuilderDrugInterpretation = ImmutableDrugsInterpretation.builder();

        for (DrugInfo drugInfo : drugs) {
            for (Drug drug : ckbEntry.drugs()) {
                if (drugInfo.id() == drug.id()) {
                    outputBuilderDrugInterpretation.drug(ImmutableDrug.builder()
                            .id(drug.id())
                            .drugName(drug.drugName())
                            .terms(drug.terms())
                            .synonyms(drug.synonyms())
                            .tradeName(drug.tradeName())
                            .drugDescriptions(createDrugDescriptions(drug.descriptions(), ckbEntry))
                            .casRegistryNum(drug.casRegistryNum())
                            .ncitId(drug.ncitId())
                            .createDate(drug.createDate())
                            .build());

                    for (DrugClassInfo drugClassInfo : drug.drugClasses()) {
                        for (DrugClass drugClass : ckbEntry.drugClasses()) {
                            if (drugClassInfo.id() == drugClass.id()) {
                                outputBuilderDrugInterpretation.drugClass(ImmutableDrugClass.builder()
                                        .id(drugClass.id())
                                        .drugClass(drugClass.drugClass())
                                        .createDate(drugClass.createDate())
                                        .build());

                            }
                        }

                    }

                }
            }

        }
        return outputBuilderDrugInterpretation.build();
    }

    @NotNull
    private static List<DrugDescription> createDrugDescriptions(@NotNull List<DescriptionInfo> descriptionInfos,
            @NotNull CkbJsonDatabase ckbEntry) {
        List<DrugDescription> drugDescriptions = Lists.newArrayList();

        for (DescriptionInfo descriptionInfo : descriptionInfos) {
            drugDescriptions.add(ImmutableDrugDescription.builder()
                    .description(descriptionInfo.description())
                    .references(extractReferences(descriptionInfo.references(), ckbEntry))
                    .build());
        }
        return drugDescriptions;
    }

    @NotNull
    private static List<ClinicalTrialLocation> extractClinicalTrialLocation(
            @NotNull List<com.hartwig.hmftools.ckb.json.clinicaltrial.ClinicalTrialLocation> clinicalTrialLocations) {
        List<ClinicalTrialLocation> clinicalTrialLocationsInterpretation = Lists.newArrayList();

        for (com.hartwig.hmftools.ckb.json.clinicaltrial.ClinicalTrialLocation location : clinicalTrialLocations) {
            clinicalTrialLocationsInterpretation.add(ImmutableClinicalTrialLocation.builder()
                    .nctId(location.nctId())
                    .facility(location.facility())
                    .city(location.city())
                    .country(location.country())
                    .status(location.status())
                    .state(location.state())
                    .zip(location.zip())
                    .contacts(extractClinicalTrialContacts(location.clinicalTrialContacts()))
                    .build());
        }
        return clinicalTrialLocationsInterpretation;
    }

    @NotNull
    private static List<com.hartwig.hmftools.ckb.datamodelinterpretation.clinicaltrial.ClinicalTrialContact> extractClinicalTrialContacts(
            @NotNull List<ClinicalTrialContact> contacts) {
        List<com.hartwig.hmftools.ckb.datamodelinterpretation.clinicaltrial.ClinicalTrialContact> clinicalTrialContacts =
                Lists.newArrayList();

        for (ClinicalTrialContact contact : contacts) {
            clinicalTrialContacts.add(ImmutableClinicalTrialContact.builder()
                    .name(contact.name())
                    .email(contact.email())
                    .phone(contact.phone())
                    .phoneExt(contact.phoneExt())
                    .role(contact.role())
                    .build());
        }
        return clinicalTrialContacts;
    }

    @NotNull
    private static List<TherapyDescription> extractTherapyDescriptions(@NotNull List<DescriptionInfo> descriptionInfos,
            @NotNull CkbJsonDatabase ckbEntry) {
        List<TherapyDescription> therapyDescriptions = Lists.newArrayList();

        for (DescriptionInfo descriptionInfo : descriptionInfos) {
            therapyDescriptions.add(ImmutableTherapyDescription.builder()
                    .description(descriptionInfo.description())
                    .references(extractReferences(descriptionInfo.references(), ckbEntry))
                    .build());
        }
        return therapyDescriptions;
    }

    @NotNull
    private static List<ReferenceExtend> extractReferences(@NotNull List<ReferenceInfo> referenceInfos, @NotNull CkbJsonDatabase ckbEntry) {
        List<ReferenceExtend> references = Lists.newArrayList();
        for (ReferenceInfo referenceInfo : referenceInfos) {
            for (Reference reference : ckbEntry.references()) {
                if (referenceInfo.id() == reference.id()) {
                    references.add(ImmutableReferenceExtend.builder()
                            .id(reference.id())
                            .pubMedId(reference.pubMedId())
                            .title(reference.title())
                            .url(reference.url())
                            .authors(reference.authors())
                            .journal(reference.journal())
                            .volume(reference.volume())
                            .issue(reference.issue())
                            .date(reference.date())
                            .abstractText(reference.abstractText())
                            .year(reference.year())
                            .build());
                }
            }
        }
        return references;
    }

    @NotNull
    private static List<ClinicalTrialVariantRequirementDetail> extractProfileName(
            @NotNull List<com.hartwig.hmftools.ckb.json.clinicaltrial.ClinicalTrialVariantRequirementDetail> molecularProfiles,
            @NotNull MolecularProfile molecularProfileDir, @NotNull CkbJsonDatabase ckbEntry) {

        List<ClinicalTrialVariantRequirementDetail> molecularProfileClinicalTrials = Lists.newArrayList();
        for (com.hartwig.hmftools.ckb.json.clinicaltrial.ClinicalTrialVariantRequirementDetail molecularProfile : molecularProfiles) {
            ImmutableVariantInterpretation.Builder outputBuilderVariantInterpretation = ImmutableVariantInterpretation.builder();
            if (molecularProfileDir.id() == molecularProfile.molecularProfile().id()) {
                for (VariantInfo variantInfo : molecularProfileDir.geneVariants()) {
                    for (Variant variant : ckbEntry.variants()) {
                        if (variantInfo.id() == variant.id()) {
                            outputBuilderVariantInterpretation.variant(ImmutableVariant.builder()
                                    .id(variant.id())
                                    .fullName(variant.fullName())
                                    .impact(variant.impact())
                                    .proteinEffect(variant.proteinEffect())
                                    .variantDescriptions(extractVariantDescriptions(variant.descriptions(), ckbEntry))
                                    .type(variant.type())
                                    .variant(variant.variant())
                                    .createDate(variant.createDate())
                                    .updateDate(variant.updateDate())
                                    .referenceTranscriptCoordinate(extractReferenceTranscriptCoordinates(variant.referenceTranscriptCoordinate()))
                                    .categoryVariantPaths(extractCategoryVariantPaths(variant.categoryVariantPaths()))
                                    .allTranscriptCoordinated(extractAllTranscriptCoordinates(variant.allTranscriptCoordinates()))
                                    .memberVariants(extractMemberVariants(variant.memberVariants(), ckbEntry))
                                    .build());

                            for (Gene gene : ckbEntry.genes()) {
                                if (variant.gene().id() == gene.id()) {
                                    outputBuilderVariantInterpretation.gene(ImmutableGene.builder()
                                            .id(gene.id())
                                            .geneSymbol(gene.geneSymbol())
                                            .terms(gene.terms())
                                            .entrezId(gene.entrezId())
                                            .synonyms(gene.synonyms())
                                            .chromosome(gene.chromosome())
                                            .mapLocation(gene.mapLocation())
                                            .geneDescriptions(extractGeneDescriptions(gene.descriptions(), ckbEntry))
                                            .canonicalTranscript(gene.canonicalTranscript())
                                            .geneRole(gene.geneRole())
                                            .createDate(gene.createDate())
                                            .updateDate(gene.updateDate())
                                            .build());
                                }
                            }
                        }
                    }
                }
            }
            molecularProfileClinicalTrials.add(ImmutableClinicalTrialVariantRequirementDetail.builder()
                    .id(molecularProfile.molecularProfile().id())
                    .profileName(molecularProfile.molecularProfile().profileName())
                    .requirementType(molecularProfile.requirementType())
                    .addVariantInterpretation(outputBuilderVariantInterpretation.build())
                    .build());

        }
        return molecularProfileClinicalTrials;
    }

    @NotNull
    private static List<GeneDescription> extractGeneDescriptions(@NotNull List<DescriptionInfo> descriptionInfos,
            @NotNull CkbJsonDatabase ckbEntry) {
        List<GeneDescription> geneDescriptions = Lists.newArrayList();

        for (DescriptionInfo descriptionInfo : descriptionInfos) {
            geneDescriptions.add(ImmutableGeneDescription.builder()
                    .description(descriptionInfo.description())
                    .references(extractReferences(descriptionInfo.references(), ckbEntry))
                    .build());
        }
        return geneDescriptions;
    }

    @NotNull
    private static List<CategoryVariantPath> extractCategoryVariantPaths(
            @NotNull List<VariantCategoryVariantPath> variantCategoryVariantPaths) {
        List<CategoryVariantPath> categoryVariantPaths = Lists.newArrayList();

        for (VariantCategoryVariantPath categoryVariantPath : variantCategoryVariantPaths) {
            categoryVariantPaths.add(ImmutableCategoryVariantPath.builder()
                    .variantPath(categoryVariantPath.variantPath())
                    .variantInfo(extractVariantInfo(categoryVariantPath.variants()))
                    .build());
        }
        return categoryVariantPaths;
    }

    @NotNull
    private static List<com.hartwig.hmftools.ckb.datamodelinterpretation.variant.VariantInfo> extractVariantInfo(
            @NotNull List<VariantInfo> variants) {
        List<com.hartwig.hmftools.ckb.datamodelinterpretation.variant.VariantInfo> variantInfos = Lists.newArrayList();

        for (VariantInfo variant : variants) {
            variantInfos.add(ImmutableVariantInfo.builder()
                    .id(variant.id())
                    .fullName(variant.fullName())
                    .impact(variant.impact())
                    .proteinEffect(variant.proteinEffect())
                    .build());
        }
        return variantInfos;
    }

    @NotNull
    private static List<MemberVariant> extractMemberVariants(@NotNull List<VariantInfo> variantsMembers,
            @NotNull CkbJsonDatabase ckbEntry) {
        List<MemberVariant> memberVariants = Lists.newArrayList();

        for (VariantInfo memberVariant : variantsMembers) {
            memberVariants.add(ImmutableMemberVariant.builder()
                    .id(memberVariant.id())
                    .fullName(memberVariant.fullName())
                    .impact(memberVariant.impact())
                    .proteinEffect(memberVariant.proteinEffect())
                    .variantDescriptions(extractVariantDescriptions(memberVariant.descriptions(), ckbEntry))
                    .build());
        }

        return memberVariants;

    }

    @NotNull
    private static List<VariantDescription> extractVariantDescriptions(@NotNull List<DescriptionInfo> descriptionInfos,
            @NotNull CkbJsonDatabase ckbEntry) {
        List<VariantDescription> variantDescriptions = Lists.newArrayList();

        for (DescriptionInfo descriptionInfo : descriptionInfos) {
            variantDescriptions.add(ImmutableVariantDescription.builder()
                    .description(descriptionInfo.description())
                    .references(extractReferences(descriptionInfo.references(), ckbEntry))
                    .build());
        }
        return variantDescriptions;
    }

    @Nullable
    private static ReferenceTranscriptCoordinate extractReferenceTranscriptCoordinates(@Nullable VariantTranscriptCoordinate coordinate) {
        if (coordinate != null) {
            return ImmutableReferenceTranscriptCoordinate.builder()
                    .id(coordinate.id())
                    .transcript(coordinate.transcript())
                    .gDna(coordinate.gDNA())
                    .cDna(coordinate.cDNA())
                    .protein(coordinate.protein())
                    .sourceDb(coordinate.sourceDB())
                    .refGenomeBuild(coordinate.refGenomeBuild())
                    .build();
        }
        return null;
    }

    @NotNull
    private static List<ReferenceTranscriptCoordinate> extractAllTranscriptCoordinates(
            @NotNull List<VariantTranscriptCoordinate> coordinates) {
        List<ReferenceTranscriptCoordinate> allTranscriptCoordinates = Lists.newArrayList();
        for (VariantTranscriptCoordinate coordinate : coordinates) {
            allTranscriptCoordinates.add(ImmutableReferenceTranscriptCoordinate.builder()
                    .id(coordinate.id())
                    .transcript(coordinate.transcript())
                    .gDna(coordinate.gDNA())
                    .cDna(coordinate.cDNA())
                    .protein(coordinate.protein())
                    .sourceDb(coordinate.sourceDB())
                    .refGenomeBuild(coordinate.refGenomeBuild())
                    .build());

        }
        return allTranscriptCoordinates;

    }

}