package com.hartwig.hmftools.knowledgebasegenerator.transvar;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

class TransvarInterpreter {

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

        // TODO Implement support for INDELs.
        if (record.gdnaRef().length() == record.gdnaAlt().length()) {
            int gdnaCodonIndex = findIndexInRefCodonForGdnaMatch(record, strand);

            for (String candidateCodon : record.candidateCodons()) {
                hotspots.add(fromCandidateCodon(record, candidateCodon, gdnaCodonIndex, strand));
            }
        }

        return hotspots;
    }

    private static int findIndexInRefCodonForGdnaMatch(@NotNull TransvarRecord record, @NotNull Strand strand) {
        String codonCompatibleRef = strand.equals(Strand.FORWARD) ? record.gdnaRef() : flipBases(record.gdnaRef());
        String codonCompatibleAlt = strand.equals(Strand.FORWARD) ? record.gdnaAlt() : flipBases(record.gdnaAlt());

        // Function only supports SNV and MNV
        assert codonCompatibleRef.length() == codonCompatibleAlt.length();

        // We look for the reference codon and candidate codon where the mutation is exclusively the mutation implied by the ref>alt
        int mutLength = codonCompatibleRef.length();
        for (String candidateCodon : record.candidateCodons()) {
            for (int i = 0; i < 4 - mutLength; i++) {
                if (record.referenceCodon().substring(i, i + mutLength).equals(codonCompatibleRef) && candidateCodon.substring(i,
                        i + mutLength).equals(codonCompatibleAlt)) {
                    boolean match = true;
                    for (int j = 0; j < 3; j++) {
                        // No match if we find another mismatch outside of the ref->alt range.
                        if ((j - i >= mutLength || j - i < 0) && !record.referenceCodon()
                                .substring(j, j + 1)
                                .equals(candidateCodon.substring(j, j + 1))) {
                            match = false;
                        }
                    }
                    if (match) {
                        return i;
                    }
                }
            }
        }

        throw new IllegalStateException("Could not find codon index for GDNA match for " + record);
    }

    @NotNull
    private static VariantHotspot fromCandidateCodon(@NotNull TransvarRecord record, @NotNull String candidateCodon, int gdnaCodonIndex,
            @NotNull Strand strand) {
        String strandAdjustedRefCodon = strand == Strand.FORWARD ? record.referenceCodon() : reverseAndFlip(record.referenceCodon());
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
            stringBuilder.append(flipBases(string.substring(i, i + 1)));
        }
        return stringBuilder.toString();
    }

    @NotNull
    private static String flipBases(@NotNull String bases) {
        StringBuilder flippedBases = new StringBuilder();
        for (char base : bases.toCharArray()) {
            flippedBases.append(flipBase(base));
        }

        return flippedBases.toString();
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
