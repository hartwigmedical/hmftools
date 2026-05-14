package com.hartwig.hmftools.redux.splice;

import static com.hartwig.hmftools.redux.splice.SpliceCommon.altContigName;

import java.util.ArrayList;
import java.util.List;

// packs a list of per-transcript sequences into a single per-chromosome alt contig. Lays transcripts down in
// input order with an N-spacer between consecutive transcripts so no real read can align across two.
// Output: the assembled sequence + one ContigEntry per packed transcript, each carrying its altStart/altEnd
// position on the alt contig. Pure function — no I/O — so it tests cleanly in isolation.
public final class AltContigPacker
{
    private final String mSpacer;

    public AltContigPacker(final int spacerLength)
    {
        mSpacer = "N".repeat(spacerLength);
    }

    public PackResult pack(final String chromosome, final List<TranscriptContigBuilder.TranscriptContigResult> transcripts)
    {
        String altContig = altContigName(chromosome);
        StringBuilder sequence = new StringBuilder();
        List<ContigEntry> entries = new ArrayList<>();

        for(TranscriptContigBuilder.TranscriptContigResult transcript : transcripts)
        {
            if(sequence.length() > 0)
                sequence.append(mSpacer);

            int altStart = sequence.length() + 1;
            sequence.append(transcript.Sequence);
            int altEnd = sequence.length();

            entries.add(new ContigEntry(
                    altContig, altStart, altEnd,
                    transcript.GeneId, transcript.GeneName, transcript.TransName, transcript.Chromosome,
                    transcript.ExonSpans));
        }

        return new PackResult(altContig, sequence.toString(), entries);
    }

    public static final class PackResult
    {
        public final String AltContig;
        public final String Sequence;
        public final List<ContigEntry> Entries;

        public PackResult(final String altContig, final String sequence, final List<ContigEntry> entries)
        {
            AltContig = altContig;
            Sequence = sequence;
            Entries = entries;
        }

        public boolean isEmpty() { return Entries.isEmpty(); }
    }
}
