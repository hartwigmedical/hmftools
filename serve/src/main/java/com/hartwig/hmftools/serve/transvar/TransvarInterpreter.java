package com.hartwig.hmftools.serve.transvar;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarDeletion;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarDuplication;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarInsertion;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarRecord;
import com.hartwig.hmftools.serve.transvar.datamodel.TransvarSnvMnv;
import com.hartwig.hmftools.serve.util.AminoAcidLookup;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

class TransvarInterpreter {

    private static final Logger LOGGER = LogManager.getLogger(TransvarInterpreter.class);

    @NotNull
    private final IndexedFastaSequenceFile refGenome;

    @NotNull
    static TransvarInterpreter fromRefGenomeFastaFile(@NotNull String refGenomeFastaFile) throws FileNotFoundException {
        return new TransvarInterpreter(new IndexedFastaSequenceFile(new File(refGenomeFastaFile)));
    }

    private TransvarInterpreter(@NotNull IndexedFastaSequenceFile refGenome) {
        this.refGenome = refGenome;
    }

    @NotNull
    List<VariantHotspot> convertRecordToHotspots(@NotNull TransvarRecord record, @NotNull Strand strand) {
        List<VariantHotspot> hotspots = Lists.newArrayList();

        if (record.annotation() instanceof TransvarSnvMnv) {
            // We need to look up which index of the ref codon is changed (0, 1 or 2) in case of SNV/MNV.
            TransvarSnvMnv snvMnv = (TransvarSnvMnv) record.annotation();
            int gdnaCodonIndex = findIndexInRefCodonForGdnaMatch(snvMnv, strand);

            for (String candidateCodon : snvMnv.candidateCodons()) {
                hotspots.add(fromCandidateCodon(record, snvMnv.referenceCodon(), candidateCodon, gdnaCodonIndex, strand));
            }
        } else {
            // For indels we assume we have to look up the base in front of the del/ins/dup and set the position 1 before the actual ref/alt
            long position = record.gdnaPosition() - 1;
            String preMutatedSequence = refGenome.getSubsequenceAt(record.chromosome(), position, position).getBaseString();

            ImmutableVariantHotspotImpl.Builder hotspotBuilder =
                    ImmutableVariantHotspotImpl.builder().chromosome(record.chromosome()).position(position);

            if (record.annotation() instanceof TransvarDuplication) {
                // Dups don't have ref and alt information so need to look it up in ref genome.
                TransvarDuplication dup = (TransvarDuplication) record.annotation();
                String dupBases = refGenome.getSubsequenceAt(record.chromosome(), position + 1, position + dup.duplicatedBaseCount())
                        .getBaseString();
                hotspotBuilder.ref(preMutatedSequence).alt(preMutatedSequence + dupBases);
                hotspots.add(hotspotBuilder.build());
            } else if (record.annotation() instanceof TransvarInsertion) {
                TransvarInsertion insertion = (TransvarInsertion) record.annotation();
                hotspotBuilder.ref(preMutatedSequence);

                // We assume inserts of length 3 are always (inframe) amino acid inserts.
                if (insertion.insertedBases().length() == 3) {
                    for (String trinucleotide : AminoAcidLookup.allTrinucleotidesForSameAminoAcid(insertion.insertedBases(), strand)) {
                        hotspots.add(hotspotBuilder.alt(preMutatedSequence + trinucleotide).build());
                    }
                } else {
                    hotspots.add(hotspotBuilder.alt(preMutatedSequence + insertion.insertedBases()).build());
                }
            } else {
                assert record.annotation() instanceof TransvarDeletion;
                TransvarDeletion deletion = (TransvarDeletion) record.annotation();
                hotspots.add(hotspotBuilder.ref(preMutatedSequence + deletion.deletedBases()).alt(preMutatedSequence).build());
            }
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
    private static VariantHotspot fromCandidateCodon(@NotNull TransvarRecord record, @NotNull String referenceCodon,
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

        return ImmutableVariantHotspotImpl.builder()
                .chromosome(record.chromosome())
                .position(record.gdnaPosition() - strandAdjustedGdnaCodingIndex + firstMutatedPosition)
                .ref(ref)
                .alt(alt)
                .build();
    }

    @NotNull
    private static String reverseAndFlip(@NotNull String string) {
        StringBuilder stringBuilder = new StringBuilder();
        for (int i = string.length() - 1; i >= 0; i--) {
            stringBuilder.append(flipBase(string.charAt(i)));
        }
        return stringBuilder.toString();
    }

    private static char flipBase(char base) {
        switch (base) {
            case 'A':
                return 'T';
            case 'T':
                return 'A';
            case 'G':
                return 'C';
            case 'C':
                return 'G';
        }

        throw new IllegalArgumentException("Cannot flip invalid base: " + base);
    }
}
