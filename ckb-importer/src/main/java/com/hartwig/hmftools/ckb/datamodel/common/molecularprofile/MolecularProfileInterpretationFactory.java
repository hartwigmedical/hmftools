package com.hartwig.hmftools.ckb.datamodel.common.molecularprofile;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ClinicalTrialVariantRequirementDetail;
import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ImmutableClinicalTrialVariantRequirementDetail;
import com.hartwig.hmftools.ckb.datamodel.common.CommonInterpretationFactory;
import com.hartwig.hmftools.ckb.datamodel.common.variant.CategoryVariantPath;
import com.hartwig.hmftools.ckb.datamodel.common.variant.GeneDescription;
import com.hartwig.hmftools.ckb.datamodel.common.variant.ImmutableCategoryVariantPath;
import com.hartwig.hmftools.ckb.datamodel.common.variant.ImmutableGene;
import com.hartwig.hmftools.ckb.datamodel.common.variant.ImmutableGeneDescription;
import com.hartwig.hmftools.ckb.datamodel.common.variant.ImmutableMemberVariant;
import com.hartwig.hmftools.ckb.datamodel.common.variant.ImmutableReferenceTranscriptCoordinate;
import com.hartwig.hmftools.ckb.datamodel.common.variant.ImmutableVariant;
import com.hartwig.hmftools.ckb.datamodel.common.variant.ImmutableVariantDescription;
import com.hartwig.hmftools.ckb.datamodel.common.variant.ImmutableVariantInfo;
import com.hartwig.hmftools.ckb.datamodel.common.variant.MemberVariant;
import com.hartwig.hmftools.ckb.datamodel.common.variant.ReferenceTranscriptCoordinate;
import com.hartwig.hmftools.ckb.datamodel.common.variant.VariantDescription;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.common.DescriptionInfo;
import com.hartwig.hmftools.ckb.json.common.VariantInfo;
import com.hartwig.hmftools.ckb.json.gene.Gene;
import com.hartwig.hmftools.ckb.json.molecularprofile.MolecularProfile;
import com.hartwig.hmftools.ckb.json.variant.Variant;
import com.hartwig.hmftools.ckb.json.variant.VariantCategoryVariantPath;
import com.hartwig.hmftools.ckb.json.variant.VariantTranscriptCoordinate;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class MolecularProfileInterpretationFactory {

    private MolecularProfileInterpretationFactory() {
    }

    @NotNull
    public static List<ClinicalTrialVariantRequirementDetail> extractProfileNameClinicalTrial(
            @NotNull List<com.hartwig.hmftools.ckb.json.clinicaltrial.ClinicalTrialVariantRequirementDetail> molecularProfiles,
            @NotNull MolecularProfile molecularProfileDir, @NotNull CkbJsonDatabase ckbEntry) {
        int countPartialRequirementTypes = 0;
        List<ClinicalTrialVariantRequirementDetail> molecularProfileClinicalTrials = Lists.newArrayList();
        for (com.hartwig.hmftools.ckb.json.clinicaltrial.ClinicalTrialVariantRequirementDetail molecularProfile : molecularProfiles) {
            if (molecularProfile.requirementType().equals("excluded")) { // variant is excluded from enrollment
                molecularProfileClinicalTrials.add(extractClinicalTrialVariantRequirementDetails(ckbEntry,
                        molecularProfile,
                        molecularProfileDir, countPartialRequirementTypes).build());
            }

            if (molecularProfile.requirementType().equals("required")) { // variant is requirement for enrollment
                molecularProfileClinicalTrials.add(extractClinicalTrialVariantRequirementDetails(ckbEntry,
                        molecularProfile,
                        molecularProfileDir, countPartialRequirementTypes).build());
            }

            if (molecularProfile.requirementType()
                    .equals("partial - required")) { // variant is required or excluded for a subset of the enrollment population
                ++countPartialRequirementTypes;
                if (molecularProfile.molecularProfile().id() == molecularProfileDir.id()) {
                    molecularProfileClinicalTrials.add(extractClinicalTrialVariantRequirementDetails(ckbEntry,
                            molecularProfile,
                            molecularProfileDir, countPartialRequirementTypes).build());
                }
            }
        }
        return molecularProfileClinicalTrials;
    }

    @NotNull
    private static ImmutableClinicalTrialVariantRequirementDetail.Builder extractClinicalTrialVariantRequirementDetails(
            @NotNull CkbJsonDatabase ckbEntry,
            @NotNull com.hartwig.hmftools.ckb.json.clinicaltrial.ClinicalTrialVariantRequirementDetail molecularProfile,
            @NotNull MolecularProfile molecularProfileDir, int countPartialRequirementTypes) {
        return ImmutableClinicalTrialVariantRequirementDetail.builder()
                .id(molecularProfile.molecularProfile().id())
                .profileName(molecularProfile.molecularProfile().profileName())
                .requirementType(molecularProfile.requirementType())
                .countPartialRequirementTypes(countPartialRequirementTypes)
                .molecularProfileInterpretation(extractVariantGeneInfo(ckbEntry, molecularProfileDir).build());
    }

    @NotNull
    public static ImmutableMolecularProfileInterpretation.Builder extractVariantGeneInfo(@NotNull CkbJsonDatabase ckbEntry,
            @NotNull MolecularProfile molecularProfileDir) {
        ImmutableMolecularProfileInterpretation.Builder outputBuilderVariantInterpretation =
                ImmutableMolecularProfileInterpretation.builder();
        for (VariantInfo variantInfo : molecularProfileDir.geneVariants()) {
            for (Variant variant : ckbEntry.variants()) {
                if (variantInfo.id() == variant.id()) {
                    outputBuilderVariantInterpretation.addVariants(ImmutableVariant.builder()
                            .id(variant.id())
                            .fullName(variant.fullName())
                            .impact(variant.impact())
                            .proteinEffect(variant.proteinEffect())
                            .variantDescriptions(extractVariantDescriptions(variant.descriptions(), ckbEntry))
                            .type(variant.type())
                            .variant(variant.variant())
                            .createDate(variant.createDate())
                            .updateDate(variant.updateDate())
                            .referenceTranscriptCoordinate(extractReferenceTranscriptCoordinate(variant.referenceTranscriptCoordinate()))
                            .categoryVariantPaths(extractCategoryVariantPaths(variant.categoryVariantPaths()))
                            .allTranscriptCoordinates(extractAllTranscriptCoordinates(variant.allTranscriptCoordinates()))
                            .memberVariants(extractMemberVariants(variant.memberVariants(), ckbEntry))
                            .gene(extractGene(ckbEntry, variant))
                            .build());
                }
            }
        }

        return outputBuilderVariantInterpretation;
    }

    @NotNull
    private static ImmutableGene extractGene(@NotNull CkbJsonDatabase ckbEntry, @NotNull Variant variant) {
        ImmutableGene.Builder outputBuilderGene = ImmutableGene.builder();
        for (Gene gene : ckbEntry.genes()) {
            if (variant.gene().id() == gene.id()) {
                outputBuilderGene.id(gene.id())
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
                        .updateDate(gene.updateDate());
            }
        }
        return outputBuilderGene.build();
    }

    @NotNull
    private static List<VariantDescription> extractVariantDescriptions(@NotNull List<DescriptionInfo> descriptionInfos,
            @NotNull CkbJsonDatabase ckbEntry) {
        List<VariantDescription> variantDescriptions = Lists.newArrayList();

        for (DescriptionInfo descriptionInfo : descriptionInfos) {
            variantDescriptions.add(ImmutableVariantDescription.builder()
                    .description(descriptionInfo.description())
                    .references(CommonInterpretationFactory.extractReferences(descriptionInfo.references(), ckbEntry))
                    .build());
        }
        return variantDescriptions;
    }

    @Nullable
    private static ReferenceTranscriptCoordinate extractReferenceTranscriptCoordinate(@Nullable VariantTranscriptCoordinate coordinate) {
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
    private static List<CategoryVariantPath> extractCategoryVariantPaths(
            @NotNull List<VariantCategoryVariantPath> variantCategoryVariantPaths) {
        List<CategoryVariantPath> categoryVariantPaths = Lists.newArrayList();

        for (VariantCategoryVariantPath categoryVariantPath : variantCategoryVariantPaths) {
            categoryVariantPaths.add(ImmutableCategoryVariantPath.builder()
                    .variantPath(categoryVariantPath.variantPath())
                    .variantInfos(extractVariantInfo(categoryVariantPath.variants()))
                    .build());
        }
        return categoryVariantPaths;
    }

    @NotNull
    private static List<com.hartwig.hmftools.ckb.datamodel.common.variant.VariantInfo> extractVariantInfo(
            @NotNull List<VariantInfo> variants) {
        List<com.hartwig.hmftools.ckb.datamodel.common.variant.VariantInfo> variantInfos = Lists.newArrayList();

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
    private static List<GeneDescription> extractGeneDescriptions(@NotNull List<DescriptionInfo> descriptionInfos,
            @NotNull CkbJsonDatabase ckbEntry) {
        List<GeneDescription> geneDescriptions = Lists.newArrayList();

        for (DescriptionInfo descriptionInfo : descriptionInfos) {
            geneDescriptions.add(ImmutableGeneDescription.builder()
                    .description(descriptionInfo.description())
                    .references(CommonInterpretationFactory.extractReferences(descriptionInfo.references(), ckbEntry))
                    .build());
        }
        return geneDescriptions;
    }
}