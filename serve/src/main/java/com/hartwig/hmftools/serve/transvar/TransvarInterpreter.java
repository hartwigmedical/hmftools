package com.hartwig.hmftools.serve.transvar;

import static com.hartwig.hmftools.serve.util.AminoAcidFunctions.reverseAndFlip;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarAnnotation;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarComplexInsertDelete;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarDeletion;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarDuplication;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarFrameshift;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarInsertion;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarRecord;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarSnvMnv;
import com.hartwig.hmftools.serve.util.AminoAcidFunctions;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

class TransvarInterpreter {

    private static final Logger LOGGER = LogManager.getLogger(TransvarInterpreter.class);

    private static final List<String> BASES = Lists.newArrayList("G", "A", "T", "C");

    @NotNull
    private final RefGenomeVersion refGenomeVersion;
    @NotNull
    private final IndexedFastaSequenceFile refGenomeFasta;

    @NotNull
    static TransvarInterpreter withRefGenome(@NotNull RefGenomeVersion refGenomeVersion, @NotNull String refGenomeFastaFile)
            throws FileNotFoundException {
        return new TransvarInterpreter(refGenomeVersion, new IndexedFastaSequenceFile(new File(refGenomeFastaFile)));
    }

    private TransvarInterpreter(@NotNull final RefGenomeVersion refGenomeVersion, @NotNull final IndexedFastaSequenceFile refGenomeFasta) {
        this.refGenomeVersion = refGenomeVersion;
        this.refGenomeFasta = refGenomeFasta;
    }

    @NotNull
    List<VariantHotspot> convertRecordToHotspots(@NotNull TransvarRecord record, @NotNull Strand strand) {
        TransvarAnnotation annotation = record.annotation();
        if (annotation instanceof TransvarSnvMnv) {
            return convertSnvMnvToHotspots(record, (TransvarSnvMnv) annotation, strand);
        } else if (annotation instanceof TransvarDuplication) {
            return convertDuplicationToHotspots(record, (TransvarDuplication) annotation);
        } else if (annotation instanceof TransvarInsertion) {
            return convertInsertionToHotspots(record, (TransvarInsertion) annotation, strand);
        } else if (annotation instanceof TransvarDeletion) {
            return convertDeletionToHotspots(record, (TransvarDeletion) annotation);
        } else if (annotation instanceof TransvarComplexInsertDelete) {
            return convertComplexInsertDeleteToHotspots(record, (TransvarComplexInsertDelete) annotation, strand);
        } else if (annotation instanceof TransvarFrameshift) {
            return convertFrameshiftToHotspots(record, (TransvarFrameshift) annotation, strand);
        } else {
            LOGGER.warn("Unrecognized annotation type in transvar record: '{}'. Skipping interpretation.",
                    annotation.getClass().toString());
            return Lists.newArrayList();
        }
    }

    @NotNull
    private List<VariantHotspot> convertSnvMnvToHotspots(@NotNull TransvarRecord record, @NotNull TransvarSnvMnv snvMnv,
            @NotNull Strand strand) {
        List<VariantHotspot> hotspots = Lists.newArrayList();

        if (!record.variantSpanMultipleExons()) {
            // We need to look up which index of the ref codon is changed (0, 1 or 2) in case of SNV/MNV.
            int gdnaCodonIndex = findIndexInRefCodonForGdnaMatch(snvMnv, strand);

            for (String candidateCodon : snvMnv.candidateCodons()) {
                hotspots.add(fromCandidateCodon(record, snvMnv.referenceCodon(), candidateCodon, gdnaCodonIndex, strand));
            }
        } else {
            LOGGER.debug("SnvMnv spanning multiple exons. Ignoring '{}'", record);
        }

        return hotspots;
    }

    private static int findIndexInRefCodonForGdnaMatch(@NotNull TransvarSnvMnv snvMnv, @NotNull Strand strand) {
        String codonCompatibleRef = strand.equals(Strand.FORWARD) ? snvMnv.gdnaRef() : reverseAndFlip(snvMnv.gdnaRef());
        String codonCompatibleAlt = strand.equals(Strand.FORWARD) ? snvMnv.gdnaAlt() : reverseAndFlip(snvMnv.gdnaAlt());

        // Function only supports SNV and MNV
        assert codonCompatibleRef.length() == codonCompatibleAlt.length();

        // We look for the candidate codon where the mutation is exclusively the mutation implied by the ref>alt
        int mutLength = codonCompatibleRef.length();
        for (String candidateCodon : snvMnv.candidateCodons()) {
            for (int i = 0; i < 4 - mutLength; i++) {
                if (snvMnv.referenceCodon().substring(i, i + mutLength).equals(codonCompatibleRef) && candidateCodon.substring(i,
                        i + mutLength).equals(codonCompatibleAlt)) {
                    boolean match = true;
                    for (int j = 0; j < 3; j++) {
                        // No match if we find another mismatch outside of the ref->alt range.
                        if ((j - i >= mutLength || j - i < 0) && !snvMnv.referenceCodon()
                                .substring(j, j + 1)
                                .equals(candidateCodon.substring(j, j + 1))) {
                            match = false;
                        }
                    }
                    if (match) {
                        // For reverse strand searches we need to correct for the flipping of MNVs.
                        return strand == Strand.FORWARD ? i : i + mutLength - 1;
                    }
                }
            }
        }

        throw new IllegalStateException("Could not find codon index for gDNA match for " + snvMnv);
    }

    @NotNull
    private VariantHotspot fromCandidateCodon(@NotNull TransvarRecord record, @NotNull String referenceCodon,
            @NotNull String candidateCodon, int gdnaCodonIndex, @NotNull Strand strand) {
        String strandAdjustedRefCodon = strand == Strand.FORWARD ? referenceCodon : reverseAndFlip(referenceCodon);
        String strandAdjustedCandidateCodon = strand == Strand.FORWARD ? candidateCodon : reverseAndFlip(candidateCodon);
        int strandAdjustedGdnaCodingIndex = strand == Strand.FORWARD ? gdnaCodonIndex : 2 - gdnaCodonIndex;

        int firstMutatedPosition = -1;
        int lastMutatedPosition = -1;
        for (int i = 0; i < 3; i++) {
            if (!strandAdjustedRefCodon.substring(i, i + 1).equals(strandAdjustedCandidateCodon.substring(i, i + 1))) {
                if (firstMutatedPosition == -1) {
                    firstMutatedPosition = i;
                }
                lastMutatedPosition = i;
            }
        }

        String ref = strandAdjustedRefCodon.substring(firstMutatedPosition, lastMutatedPosition + 1);
        String alt = strandAdjustedCandidateCodon.substring(firstMutatedPosition, lastMutatedPosition + 1);

        return withRefBasedChromosome(record.chromosome()).position(
                record.gdnaPosition() - strandAdjustedGdnaCodingIndex + firstMutatedPosition).ref(ref).alt(alt).build();
    }

    @NotNull
    private List<VariantHotspot> convertDuplicationToHotspots(@NotNull TransvarRecord record, @NotNull TransvarDuplication dup) {
        List<VariantHotspot> hotspots = Lists.newArrayList();
        if (!record.variantSpanMultipleExons()) {
            // Dups don't have ref and alt information so need to look it up in ref genome.
            long position = record.gdnaPosition() - 1;
            String preMutatedSequence = refSequence(record.chromosome(), position, position);
            String dupBases = refSequence(record.chromosome(), position + 1, position + dup.duplicatedBaseCount());

            hotspots.add(withRefBasedChromosome(record.chromosome()).position(position)
                    .ref(preMutatedSequence)
                    .alt(preMutatedSequence + dupBases)
                    .build());
        } else {
            LOGGER.debug("Duplication spanning multiple exons. Ignoring '{}'", record);
        }

        return hotspots;
    }

    @NotNull
    private List<VariantHotspot> convertInsertionToHotspots(@NotNull TransvarRecord record, @NotNull TransvarInsertion insertion,
            @NotNull Strand strand) {
        List<VariantHotspot> hotspots = Lists.newArrayList();

        long position = record.gdnaPosition();
        String preMutatedSequence = refSequence(record.chromosome(), position, position);

        ImmutableVariantHotspotImpl.Builder hotspotBuilder =
                withRefBasedChromosome(record.chromosome()).position(position).ref(preMutatedSequence);

        // We assume inserts of length 3 are always (inframe) amino acid inserts,
        //  and we know the inserted bases match the proteins inserted.
        if (insertion.insertedBases().length() == 3) {
            for (String trinucleotide : AminoAcidFunctions.allTrinucleotidesForSameAminoAcid(insertion.insertedBases(), strand)) {
                hotspots.add(hotspotBuilder.alt(preMutatedSequence + trinucleotide).build());
            }
        } else {
            hotspots.add(hotspotBuilder.alt(preMutatedSequence + insertion.insertedBases()).build());
        }

        return hotspots;
    }

    @NotNull
    private List<VariantHotspot> convertDeletionToHotspots(@NotNull TransvarRecord record, @NotNull TransvarDeletion deletion) {
        List<VariantHotspot> hotspots = Lists.newArrayList();
        if (!record.variantSpanMultipleExons()) {
            ImmutableVariantHotspotImpl.Builder hotspotBuilder = withRefBasedChromosome(record.chromosome());
            for (long start = deletion.leftAlignedGDNAPosition(); start <= record.gdnaPosition(); start++) {
                long adjustedPosition = start - 1;
                String preMutatedSequence = refSequence(record.chromosome(), adjustedPosition, adjustedPosition);
                String deletedSequence =
                        refSequence(record.chromosome(), adjustedPosition + 1, adjustedPosition + deletion.deletedBaseCount());

                hotspots.add(hotspotBuilder.position(adjustedPosition)
                        .ref(preMutatedSequence + deletedSequence)
                        .alt(preMutatedSequence)
                        .build());
            }
        } else {
            LOGGER.debug("Deletion spanning multiple exons. Ignoring '{}'", record);
        }

        return hotspots;
    }

    @NotNull
    private List<VariantHotspot> convertComplexInsertDeleteToHotspots(@NotNull TransvarRecord record,
            @NotNull TransvarComplexInsertDelete insDel, @NotNull Strand strand) {
        List<VariantHotspot> hotspots = Lists.newArrayList();
        if (!record.variantSpanMultipleExons()) {
            long position = record.gdnaPosition();
            String deletedBases = refSequence(record.chromosome(), position, position + insDel.deletedBaseCount() - 1);

            ImmutableVariantHotspotImpl.Builder hotspotBuilder =
                    withRefBasedChromosome(record.chromosome()).position(position).ref(deletedBases);

            String insertedBases = insDel.insertedSequence();
            hotspots.add(reduceComplexityForComplexInsDel(hotspotBuilder.alt(insertedBases).build()));
            // For now we only add alternative sequences for insertions of one AA
            if (insertedBases.length() == 3) {
                for (String candidateAlternativeCodon : insDel.candidateAlternativeCodons()) {
                    assert candidateAlternativeCodon.length() == insertedBases.length();
                    String strandAdjustedAlternativeSequence =
                            strand == Strand.FORWARD ? candidateAlternativeCodon : reverseAndFlip(candidateAlternativeCodon);
                    if (!strandAdjustedAlternativeSequence.equals(insertedBases)) {
                        hotspots.add(reduceComplexityForComplexInsDel(hotspotBuilder.alt(strandAdjustedAlternativeSequence).build()));
                    }
                }
            }
        } else {
            LOGGER.debug("Complex insert/delete spanning multiple exons. Ignoring '{}'", record);
        }
        return hotspots;
    }

    @NotNull
    @VisibleForTesting
    VariantHotspot reduceComplexityForComplexInsDel(@NotNull VariantHotspot complexInsDel) {
        assert complexInsDel.ref().length() > 1 && complexInsDel.alt().length() > 1;

        String simplifiedRef = complexInsDel.ref();
        String simplifiedAlt = complexInsDel.alt();
        while (simplifiedRef.length() > 1 && simplifiedAlt.length() > 1 && simplifiedRef.charAt(simplifiedRef.length() - 1) == simplifiedAlt
                .charAt(simplifiedAlt.length() - 1)) {
            simplifiedRef = simplifiedRef.substring(0, simplifiedRef.length() - 1);
            simplifiedAlt = simplifiedAlt.substring(0, simplifiedAlt.length() - 1);
        }

        long adjustedPos = complexInsDel.position();
        while (simplifiedRef.length() > 2 && simplifiedAlt.length() > 2 && simplifiedRef.subSequence(0, 1)
                .equals(simplifiedAlt.substring(0, 1))) {
            adjustedPos++;
            simplifiedRef = simplifiedRef.substring(1);
            simplifiedAlt = simplifiedAlt.substring(1);
        }

        return withRefBasedChromosome(complexInsDel.chromosome())
                .position(adjustedPos)
                .ref(simplifiedRef)
                .alt(simplifiedAlt)
                .build();
    }

    @NotNull
    private List<VariantHotspot> convertFrameshiftToHotspots(@NotNull TransvarRecord record, @NotNull TransvarFrameshift frameshift,
            @NotNull Strand strand) {
        List<VariantHotspot> hotspots = Lists.newArrayList();
        if (!record.variantSpanMultipleExons()) {
            long posPriorToCodon = strand == Strand.FORWARD ? record.gdnaPosition() : record.gdnaPosition() - 3;
            // For frameshifts in start codons, transvar generates the start of the start codon rather than position prior.
            if (frameshift.isFrameshiftInsideStartCodon()) {
                posPriorToCodon = strand == Strand.FORWARD ? posPriorToCodon - 1 : posPriorToCodon + 1;
            }
            String referenceCodon = refSequence(record.chromosome(), posPriorToCodon + 1, posPriorToCodon + 3);
            String refAminoAcid =
                    AminoAcidFunctions.findAminoAcidForCodon(strand == Strand.FORWARD ? referenceCodon : reverseAndFlip(referenceCodon));

            if (refAminoAcid == null) {
                LOGGER.warn("Could not resolve a valid ref amino acid for '{}' based on reference codon {}. Skipping hotspot generation.",
                        record,
                        referenceCodon);
                return Lists.newArrayList();
            }

            hotspots.addAll(generateSingleBaseInserts(record, posPriorToCodon, strand, refAminoAcid));
            hotspots.addAll(generateSingleBaseDeletes(record, posPriorToCodon, strand, refAminoAcid));
            hotspots.addAll(generateDoubleBaseDeletes(record, posPriorToCodon, strand, refAminoAcid));
        } else {
            LOGGER.debug("Frameshift spanning multiple exons. Ignoring '{}'", record);
        }

        return hotspots;
    }

    @NotNull
    private List<VariantHotspot> generateSingleBaseInserts(@NotNull TransvarRecord record, long posPriorToCodon, @NotNull Strand strand,
            @NotNull String refAminoAcid) {
        List<VariantHotspot> hotspots = Lists.newArrayList();
        ImmutableVariantHotspotImpl.Builder builder = withRefBasedChromosome(record.chromosome());

        // Add 12 single base insertions in case they don't lead to synonymous impact in the impacted codon
        for (int i = 0; i < 3; i++) {
            // For reverse strand we need to move the position up by 1 for inserts
            int strandCorrection = strand == Strand.FORWARD ? 0 : 1;
            long pos = posPriorToCodon + i + strandCorrection;
            String ref = refSequence(record.chromosome(), pos, pos);
            builder.position(pos).ref(ref);

            String refBase1 =
                    refSequence(record.chromosome(), posPriorToCodon + 1 + strandCorrection, posPriorToCodon + 1 + strandCorrection);
            String refBase2 =
                    refSequence(record.chromosome(), posPriorToCodon + 2 + strandCorrection, posPriorToCodon + 2 + strandCorrection);
            for (String base : BASES) {
                String newRefCodon;
                if (i == 0) {
                    newRefCodon = base + refBase1 + refBase2;
                } else if (i == 1) {
                    newRefCodon = refBase1 + base + refBase2;
                } else {
                    newRefCodon = refBase1 + refBase2 + base;
                }
                String newAminoAcid =
                        AminoAcidFunctions.findAminoAcidForCodon(strand == Strand.FORWARD ? newRefCodon : reverseAndFlip(newRefCodon));

                if (newAminoAcid != null && !newAminoAcid.equals(refAminoAcid)) {
                    hotspots.add(builder.alt(ref + base).build());
                }
            }
        }
        return hotspots;
    }

    @NotNull
    private List<VariantHotspot> generateSingleBaseDeletes(@NotNull TransvarRecord record, long posPriorToCodon, @NotNull Strand strand,
            @NotNull String refAminoAcid) {
        List<VariantHotspot> hotspots = Lists.newArrayList();
        ImmutableVariantHotspotImpl.Builder builder = withRefBasedChromosome(record.chromosome());

        // Add the 3 single base deletes in case they don't lead to synonymous impact in the impacted codon
        for (int i = 0; i < 3; i++) {
            String newRefCodon;
            if (strand == Strand.FORWARD) {
                if (i == 0) {
                    newRefCodon = refSequence(record.chromosome(), posPriorToCodon + 2, posPriorToCodon + 4);
                } else if (i == 1) {
                    newRefCodon =
                            refSequence(record.chromosome(), posPriorToCodon + 1, posPriorToCodon + 1) + refSequence(record.chromosome(),
                                    posPriorToCodon + 3,
                                    posPriorToCodon + 4);
                } else {
                    newRefCodon =
                            refSequence(record.chromosome(), posPriorToCodon + 1, posPriorToCodon + 2) + refSequence(record.chromosome(),
                                    posPriorToCodon + 4,
                                    posPriorToCodon + 4);
                }
            } else {
                if (i == 0) {
                    newRefCodon = refSequence(record.chromosome(), posPriorToCodon, posPriorToCodon) + refSequence(record.chromosome(),
                            posPriorToCodon + 2,
                            posPriorToCodon + 3);
                } else if (i == 1) {
                    newRefCodon = refSequence(record.chromosome(), posPriorToCodon, posPriorToCodon + 1) + refSequence(record.chromosome(),
                            posPriorToCodon + 3,
                            posPriorToCodon + 3);
                } else {
                    newRefCodon = refSequence(record.chromosome(), posPriorToCodon, posPriorToCodon + 2);
                }
            }
            String newAminoAcid =
                    AminoAcidFunctions.findAminoAcidForCodon(strand == Strand.FORWARD ? newRefCodon : reverseAndFlip(newRefCodon));

            if (newAminoAcid != null && !newAminoAcid.equals(refAminoAcid)) {
                long pos = posPriorToCodon + i;
                String ref = refSequence(record.chromosome(), pos, pos + 1);
                String alt = refSequence(record.chromosome(), pos, pos);

                hotspots.add(builder.position(pos).ref(ref).alt(alt).build());
            }
        }
        return hotspots;
    }

    @NotNull
    private List<VariantHotspot> generateDoubleBaseDeletes(@NotNull TransvarRecord record, long posPriorToCodon, @NotNull Strand strand,
            @NotNull String refAminoAcid) {
        List<VariantHotspot> hotspots = Lists.newArrayList();
        ImmutableVariantHotspotImpl.Builder builder = withRefBasedChromosome(record.chromosome());

        // Add the 2 double base deletes in case they don't lead to synonymous impact in the impacted codon
        for (int i = 0; i < 2; i++) {
            String newRefCodon;
            if (strand == Strand.FORWARD) {
                if (i == 0) {
                    newRefCodon = refSequence(record.chromosome(), posPriorToCodon + 3, posPriorToCodon + 5);
                } else {
                    newRefCodon =
                            refSequence(record.chromosome(), posPriorToCodon + 1, posPriorToCodon + 1) + refSequence(record.chromosome(),
                                    posPriorToCodon + 4,
                                    posPriorToCodon + 5);
                }
            } else {
                if (i == 0) {
                    newRefCodon = refSequence(record.chromosome(), posPriorToCodon - 1, posPriorToCodon) + refSequence(record.chromosome(),
                            posPriorToCodon + 3,
                            posPriorToCodon + 3);
                } else {
                    newRefCodon = refSequence(record.chromosome(), posPriorToCodon - 1, posPriorToCodon + 1);
                }
            }

            String newAminoAcid =
                    AminoAcidFunctions.findAminoAcidForCodon(strand == Strand.FORWARD ? newRefCodon : reverseAndFlip(newRefCodon));

            if (newAminoAcid != null && !newAminoAcid.equals(refAminoAcid)) {
                long pos = posPriorToCodon + i;
                String ref = refSequence(record.chromosome(), pos, pos + 2);
                String alt = refSequence(record.chromosome(), pos, pos);

                hotspots.add(builder.position(pos).ref(ref).alt(alt).build());
            }
        }
        return hotspots;
    }

    @NotNull
    private String refSequence(@NotNull String chromosome, long start, long end) {
        String versionedChromosome = refGenomeVersion.versionedChromosome(chromosome);
        return refGenomeFasta.getSubsequenceAt(versionedChromosome, start, end).getBaseString();
    }

    @NotNull
    private ImmutableVariantHotspotImpl.Builder withRefBasedChromosome(@NotNull String chromosome) {
        return ImmutableVariantHotspotImpl.builder().chromosome(refGenomeVersion.versionedChromosome(chromosome));
    }
}
