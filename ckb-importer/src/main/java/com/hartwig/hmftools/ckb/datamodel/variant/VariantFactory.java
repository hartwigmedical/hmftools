package com.hartwig.hmftools.ckb.datamodel.variant;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodel.reference.ReferenceFactory;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.common.DescriptionInfo;
import com.hartwig.hmftools.ckb.json.common.GeneInfo;
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
            variants.add(resolveVariant(ckbJsonDatabase, variantInfo));
        }

        return variants;
    }

    @NotNull
    private static Variant resolveVariant(@NotNull CkbJsonDatabase ckbJsonDatabase, @NotNull VariantInfo variantInfo) {
        for (JsonVariant variant : ckbJsonDatabase.variants()) {
            if (variant.id() == variantInfo.id()) {
                return ImmutableVariant.builder()
                        .id(variant.id())
                        .gene(resolveGene(ckbJsonDatabase, variant.gene()))
                        .fullName(variant.fullName())
                        .impact(variant.impact())
                        .proteinEffect(variant.proteinEffect())
                        .descriptions(extractVariantDescriptions(ckbJsonDatabase, variant.descriptions()))
                        .type(variant.type())
                        .variant(variant.variant())
                        .createDate(variant.createDate())
                        .updateDate(variant.updateDate())
                        .referenceTranscriptCoordinate(convertReferenceTranscriptCoordinate(variant.referenceTranscriptCoordinate()))
                        .categoryVariantPaths(convertCategoryVariantPaths(variant.categoryVariantPaths()))
                        .allTranscriptCoordinates(convertAllTranscriptCoordinates(variant.allTranscriptCoordinates()))
                        .memberVariants(extractMemberVariants(ckbJsonDatabase, variant.memberVariants()))
                        .build();
            }
        }

        throw new IllegalStateException("Could not resolve CKB variant with id '" + variantInfo.id() + "'");
    }

    @NotNull
    private static Gene resolveGene(@NotNull CkbJsonDatabase ckbJsonDatabase, @NotNull GeneInfo geneInfo) {
        for (JsonGene gene : ckbJsonDatabase.genes()) {
            if (gene.id() == geneInfo.id()) {
                return ImmutableGene.builder()
                        .id(gene.id())
                        .geneSymbol(gene.geneSymbol())
                        .terms(gene.terms())
                        .entrezId(gene.entrezId())
                        .synonyms(gene.synonyms())
                        .chromosome(gene.chromosome())
                        .mapLocation(gene.mapLocation())
                        .descriptions(extractGeneDescriptions(ckbJsonDatabase, gene.descriptions()))
                        .canonicalTranscript(gene.canonicalTranscript())
                        .geneRole(gene.geneRole())
                        .createDate(gene.createDate())
                        .updateDate(gene.updateDate())
                        .build();
            }
        }

        throw new IllegalStateException("Could not resolve CKB gene with id '" + geneInfo.id() + "'");
    }

    @Nullable
    private static TranscriptCoordinate convertReferenceTranscriptCoordinate(@Nullable JsonVariantTranscriptCoordinate coordinate) {
        if (coordinate == null) {
            return null;
        }

        return ImmutableTranscriptCoordinate.builder()
                .id(coordinate.id())
                .transcript(coordinate.transcript())
                .gDna(coordinate.gDNA())
                .cDna(coordinate.cDNA())
                .protein(coordinate.protein())
                .sourceDb(coordinate.sourceDB())
                .refGenomeBuild(coordinate.refGenomeBuild())
                .build();
    }

    @NotNull
    private static List<CategoryVariantPath> convertCategoryVariantPaths(
            @NotNull List<JsonVariantCategoryVariantPath> variantCategoryVariantPaths) {
        List<CategoryVariantPath> categoryVariantPaths = Lists.newArrayList();

        for (JsonVariantCategoryVariantPath categoryVariantPath : variantCategoryVariantPaths) {
            categoryVariantPaths.add(ImmutableCategoryVariantPath.builder()
                    .variantPath(categoryVariantPath.variantPath())
                    .variants(convertCategoryVariants(categoryVariantPath.variants()))
                    .build());
        }

        return categoryVariantPaths;
    }

    @NotNull
    private static List<CategoryVariant> convertCategoryVariants(@NotNull List<VariantInfo> variants) {
        List<CategoryVariant> categoryVariants = Lists.newArrayList();

        for (VariantInfo variant : variants) {
            categoryVariants.add(ImmutableCategoryVariant.builder()
                    .id(variant.id())
                    .fullName(variant.fullName())
                    .impact(variant.impact())
                    .proteinEffect(variant.proteinEffect())
                    .build());
        }

        return categoryVariants;
    }

    @NotNull
    private static List<TranscriptCoordinate> convertAllTranscriptCoordinates(@NotNull List<JsonVariantTranscriptCoordinate> coordinates) {
        List<TranscriptCoordinate> allTranscriptCoordinates = Lists.newArrayList();

        for (JsonVariantTranscriptCoordinate coordinate : coordinates) {
            allTranscriptCoordinates.add(ImmutableTranscriptCoordinate.builder()
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
                    .descriptions(extractVariantDescriptions(ckbJsonDatabase, memberVariant.descriptions()))
                    .build());
        }

        return memberVariants;
    }

    @NotNull
    private static List<VariantDescription> extractVariantDescriptions(@NotNull CkbJsonDatabase ckbJsonDatabase,
            @NotNull List<DescriptionInfo> descriptionInfos) {
        List<VariantDescription> variantDescriptions = Lists.newArrayList();

        for (DescriptionInfo descriptionInfo : descriptionInfos) {
            variantDescriptions.add(ImmutableVariantDescription.builder()
                    .description(descriptionInfo.description())
                    .references(ReferenceFactory.extractReferences(ckbJsonDatabase, descriptionInfo.references()))
                    .build());
        }

        return variantDescriptions;
    }

    @NotNull
    private static List<GeneDescription> extractGeneDescriptions(@NotNull CkbJsonDatabase ckbJsonDatabase,
            @NotNull List<DescriptionInfo> descriptionInfos) {
        List<GeneDescription> geneDescriptions = Lists.newArrayList();

        for (DescriptionInfo descriptionInfo : descriptionInfos) {
            geneDescriptions.add(ImmutableGeneDescription.builder()
                    .description(descriptionInfo.description())
                    .references(ReferenceFactory.extractReferences(ckbJsonDatabase, descriptionInfo.references()))
                    .build());
        }

        return geneDescriptions;
    }
}