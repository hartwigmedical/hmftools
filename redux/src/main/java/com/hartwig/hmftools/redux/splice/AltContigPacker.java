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
        final String altContig = altContigName(chromosome);
        final StringBuilder sequence = new StringBuilder();
        final List<ContigEntry> entries = new ArrayList<>();

        for(final TranscriptContigBuilder.TranscriptContigResult transcript : transcripts)
        {
            if(sequence.length() > 0)
                sequence.append(mSpacer);

            final int altStart = sequence.length() + 1;
            sequence.append(transcript.sequence());
            final int altEnd = sequence.length();

            entries.add(new ContigEntry(
                    altContig, altStart, altEnd,
                    transcript.geneId(), transcript.geneName(), transcript.transName(), transcript.chromosome(),
                    transcript.exonSpans()));
        }

        return new PackResult(altContig, sequence.toString(), entries);
    }

    public record PackResult(String altContig, String sequence, List<ContigEntry> entries)
    {
        public boolean isEmpty() { return entries.isEmpty(); }
    }
}
