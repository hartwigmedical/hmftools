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
                        .createDate(variant.createDate())
                        .updateDate(variant.updateDate())
                        .fullName(variant.fullName())
                        .variant(variant.variant())
                        .impact(variant.impact())
                        .proteinEffect(variant.proteinEffect())
                        .type(variant.type())
                        .gene(resolveGene(ckbJsonDatabase, variant.gene()))
                        .referenceTranscriptCoordinate(convertReferenceTranscriptCoordinate(variant.referenceTranscriptCoordinate()))
                        .allTranscriptCoordinates(convertAllTranscriptCoordinates(variant.allTranscriptCoordinates()))
                        .categoryVariantPaths(convertCategoryVariantPaths(variant.categoryVariantPaths()))
                        .memberVariants(convertMemberVariants(variant.memberVariants()))
                        .descriptions(extractVariantDescriptions(ckbJsonDatabase, variant.descriptions()))
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
                        .createDate(gene.createDate())
                        .updateDate(gene.updateDate())
                        .geneSymbol(gene.geneSymbol())
                        .geneRole(gene.geneRole())
                        .entrezId(gene.entrezId())
                        .chromosome(gene.chromosome())
                        .mapLocation(gene.mapLocation())
                        .canonicalTranscript(gene.canonicalTranscript())
                        .terms(gene.terms())
                        .synonyms(gene.synonyms())
                        .descriptions(extractGeneDescriptions(ckbJsonDatabase, gene.descriptions()))
                        .build();
            }
        }

        throw new IllegalStateException("Could not resolve CKB gene with id '" + geneInfo.id() + "'");
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
    private static List<String> convertCategoryVariantPaths(@NotNull List<JsonVariantCategoryVariantPath> variantCategoryVariantPaths) {
        List<String> categoryVariantPaths = Lists.newArrayList();

        for (JsonVariantCategoryVariantPath categoryVariantPath : variantCategoryVariantPaths) {
            categoryVariantPaths.add(categoryVariantPath.variantPath());
        }

        return categoryVariantPaths;
    }

    @NotNull
    private static List<MemberVariant> convertMemberVariants(@NotNull List<VariantInfo> variantsMembers) {
        List<MemberVariant> memberVariants = Lists.newArrayList();

        for (VariantInfo memberVariant : variantsMembers) {
            memberVariants.add(ImmutableMemberVariant.builder()
                    .id(memberVariant.id())
                    .fullName(memberVariant.fullName())
                    .impact(memberVariant.impact())
                    .proteinEffect(memberVariant.proteinEffect())
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
}