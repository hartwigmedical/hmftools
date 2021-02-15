package com.hartwig.hmftools.ckb.datamodel.common.molecularprofile;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ClinicalTrialVariantRequirementDetail;
import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ImmutableClinicalTrialVariantRequirementDetail;
import com.hartwig.hmftools.ckb.datamodel.common.CommonInterpretationFactory;
import com.hartwig.hmftools.ckb.datamodel.common.variant.CategoryVariantPath;
import com.hartwig.hmftools.ckb.datamodel.common.variant.Gene;
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
import com.hartwig.hmftools.ckb.datamodel.common.variant.Variant;
import com.hartwig.hmftools.ckb.datamodel.common.variant.VariantDescription;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.clinicaltrial.JsonClinicalTrialVariantRequirementDetail;
import com.hartwig.hmftools.ckb.json.common.DescriptionInfo;
import com.hartwig.hmftools.ckb.json.common.VariantInfo;
import com.hartwig.hmftools.ckb.json.gene.JsonGene;
import com.hartwig.hmftools.ckb.json.molecularprofile.JsonMolecularProfile;
import com.hartwig.hmftools.ckb.json.variant.JsonVariant;
import com.hartwig.hmftools.ckb.json.variant.JsonVariantCategoryVariantPath;
import com.hartwig.hmftools.ckb.json.variant.JsonVariantTranscriptCoordinate;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class MolecularProfileInterpretationFactory {

    private MolecularProfileInterpretationFactory() {
    }

    @NotNull
    public static List<ClinicalTrialVariantRequirementDetail> extractProfileNameClinicalTrial(@NotNull CkbJsonDatabase ckbJsonDatabase,
            @NotNull List<JsonClinicalTrialVariantRequirementDetail> molecularProfiles, @NotNull JsonMolecularProfile molecularProfileDir) {
        int countPartialRequirementTypes = 0;
        List<ClinicalTrialVariantRequirementDetail> molecularProfileClinicalTrials = Lists.newArrayList();
        for (JsonClinicalTrialVariantRequirementDetail molecularProfile : molecularProfiles) {
            if (molecularProfile.requirementType().equals("excluded")) { // variant is excluded from enrollment
                molecularProfileClinicalTrials.add(extractClinicalTrialVariantRequirementDetails(ckbJsonDatabase,
                        molecularProfile,
                        molecularProfileDir,
                        countPartialRequirementTypes).build());
            }

            if (molecularProfile.requirementType().equals("required")) { // variant is requirement for enrollment
                molecularProfileClinicalTrials.add(extractClinicalTrialVariantRequirementDetails(ckbJsonDatabase,
                        molecularProfile,
                        molecularProfileDir,
                        countPartialRequirementTypes).build());
            }

            if (molecularProfile.requirementType()
                    .equals("partial - required")) { // variant is required or excluded for a subset of the enrollment population
                ++countPartialRequirementTypes;
                if (molecularProfile.molecularProfile().id() == molecularProfileDir.id()) {
                    molecularProfileClinicalTrials.add(extractClinicalTrialVariantRequirementDetails(ckbJsonDatabase,
                            molecularProfile,
                            molecularProfileDir,
                            countPartialRequirementTypes).build());
                }
            }
        }
        return molecularProfileClinicalTrials;
    }

    @NotNull
    private static ImmutableClinicalTrialVariantRequirementDetail.Builder extractClinicalTrialVariantRequirementDetails(
            @NotNull CkbJsonDatabase ckbJsonDatabase, @NotNull JsonClinicalTrialVariantRequirementDetail molecularProfile,
            @NotNull JsonMolecularProfile molecularProfileDir, int countPartialRequirementTypes) {
        return ImmutableClinicalTrialVariantRequirementDetail.builder()
                .id(molecularProfile.molecularProfile().id())
                .profileName(molecularProfile.molecularProfile().profileName())
                .requirementType(molecularProfile.requirementType())
                .countPartialRequirementTypes(countPartialRequirementTypes)
                .variants(extractVariantGeneInfo(ckbJsonDatabase, molecularProfileDir));
    }

    @NotNull
    public static List<Variant> extractVariantGeneInfo(@NotNull CkbJsonDatabase ckbJsonDatabase,
            @NotNull JsonMolecularProfile molecularProfileDir) {
        List<Variant> variants = Lists.newArrayList();
        for (VariantInfo variantInfo : molecularProfileDir.geneVariants()) {
            for (JsonVariant variant : ckbJsonDatabase.variants()) {
                if (variantInfo.id() == variant.id()) {
                    variants.add(ImmutableVariant.builder()
                            .id(variant.id())
                            .fullName(variant.fullName())
                            .impact(variant.impact())
                            .proteinEffect(variant.proteinEffect())
                            .variantDescriptions(extractVariantDescriptions(ckbJsonDatabase, variant.descriptions()))
                            .type(variant.type())
                            .variant(variant.variant())
                            .createDate(variant.createDate())
                            .updateDate(variant.updateDate())
                            .referenceTranscriptCoordinate(extractReferenceTranscriptCoordinate(variant.referenceTranscriptCoordinate()))
                            .categoryVariantPaths(extractCategoryVariantPaths(variant.categoryVariantPaths()))
                            .allTranscriptCoordinates(extractAllTranscriptCoordinates(variant.allTranscriptCoordinates()))
                            .memberVariants(extractMemberVariants(ckbJsonDatabase, variant.memberVariants()))
                            .gene(extractGene(ckbJsonDatabase, variant))
                            .build());
                }
            }
        }

        return variants;
    }

    @NotNull
    private static Gene extractGene(@NotNull CkbJsonDatabase ckbJsonDatabase, @NotNull JsonVariant variant) {
        ImmutableGene.Builder outputBuilderGene = ImmutableGene.builder();
        for (JsonGene gene : ckbJsonDatabase.genes()) {
            if (variant.gene().id() == gene.id()) {
                outputBuilderGene.id(gene.id())
                        .geneSymbol(gene.geneSymbol())
                        .terms(gene.terms())
                        .entrezId(gene.entrezId())
                        .synonyms(gene.synonyms())
                        .chromosome(gene.chromosome())
                        .mapLocation(gene.mapLocation())
                        .geneDescriptions(extractGeneDescriptions(ckbJsonDatabase, gene.descriptions()))
                        .canonicalTranscript(gene.canonicalTranscript())
                        .geneRole(gene.geneRole())
                        .createDate(gene.createDate())
                        .updateDate(gene.updateDate());
            }
        }
        return outputBuilderGene.build();
    }

    @NotNull
    private static List<VariantDescription> extractVariantDescriptions(@NotNull CkbJsonDatabase ckbJsonDatabase,
            @NotNull List<DescriptionInfo> descriptionInfos) {
        List<VariantDescription> variantDescriptions = Lists.newArrayList();

        for (DescriptionInfo descriptionInfo : descriptionInfos) {
            variantDescriptions.add(ImmutableVariantDescription.builder()
                    .description(descriptionInfo.description())
                    .references(CommonInterpretationFactory.extractReferences(ckbJsonDatabase, descriptionInfo.references()))
                    .build());
        }
        return variantDescriptions;
    }

    @Nullable
    private static ReferenceTranscriptCoordinate extractReferenceTranscriptCoordinate(
            @Nullable JsonVariantTranscriptCoordinate coordinate) {
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
            @NotNull List<JsonVariantCategoryVariantPath> variantCategoryVariantPaths) {
        List<CategoryVariantPath> categoryVariantPaths = Lists.newArrayList();

        for (JsonVariantCategoryVariantPath categoryVariantPath : variantCategoryVariantPaths) {
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
            @NotNull List<JsonVariantTranscriptCoordinate> coordinates) {
        List<ReferenceTranscriptCoordinate> allTranscriptCoordinates = Lists.newArrayList();
        for (JsonVariantTranscriptCoordinate coordinate : coordinates) {
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
    private static List<MemberVariant> extractMemberVariants(@NotNull CkbJsonDatabase ckbJsonDatabase,
            @NotNull List<VariantInfo> variantsMembers) {
        List<MemberVariant> memberVariants = Lists.newArrayList();

        for (VariantInfo memberVariant : variantsMembers) {
            memberVariants.add(ImmutableMemberVariant.builder()
                    .id(memberVariant.id())
                    .fullName(memberVariant.fullName())
                    .impact(memberVariant.impact())
                    .proteinEffect(memberVariant.proteinEffect())
                    .variantDescriptions(extractVariantDescriptions(ckbJsonDatabase, memberVariant.descriptions()))
                    .build());
        }

        return memberVariants;
    }

    @NotNull
    private static List<GeneDescription> extractGeneDescriptions(@NotNull CkbJsonDatabase ckbJsonDatabase,
            @NotNull List<DescriptionInfo> descriptionInfos) {
        List<GeneDescription> geneDescriptions = Lists.newArrayList();

        for (DescriptionInfo descriptionInfo : descriptionInfos) {
            geneDescriptions.add(ImmutableGeneDescription.builder()
                    .description(descriptionInfo.description())
                    .references(CommonInterpretationFactory.extractReferences(ckbJsonDatabase, descriptionInfo.references()))
                    .build());
        }
        return geneDescriptions;
    }
}