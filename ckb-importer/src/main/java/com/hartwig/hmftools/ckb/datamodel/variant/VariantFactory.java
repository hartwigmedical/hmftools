package com.hartwig.hmftools.ckb.datamodel.variant;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodel.CommonInterpretationFactory;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.common.DescriptionInfo;
import com.hartwig.hmftools.ckb.json.common.VariantInfo;
import com.hartwig.hmftools.ckb.json.gene.JsonGene;
import com.hartwig.hmftools.ckb.json.variant.JsonVariant;
import com.hartwig.hmftools.ckb.json.variant.JsonVariantCategoryVariantPath;
import com.hartwig.hmftools.ckb.json.variant.JsonVariantTranscriptCoordinate;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class VariantFactory {

    private VariantFactory() {
    }

    @NotNull
    public static List<Variant> extractVariants(@NotNull CkbJsonDatabase ckbJsonDatabase, @NotNull List<VariantInfo> variantInfos) {
        List<Variant> variants = Lists.newArrayList();
        for (VariantInfo variantInfo : variantInfos) {
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
    private static List<com.hartwig.hmftools.ckb.datamodel.variant.VariantInfo> extractVariantInfo(
            @NotNull List<VariantInfo> variants) {
        List<com.hartwig.hmftools.ckb.datamodel.variant.VariantInfo> variantInfos = Lists.newArrayList();

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