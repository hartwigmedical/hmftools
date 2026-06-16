package com.hartwig.hmftools.tars.fasta;

import static com.hartwig.hmftools.tars.common.SpliceCommon.altContigName;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.tars.common.ContigEntry;

// Packs per-transcript sequences into a single per-chromosome alt contig with N-spacers between transcripts.
// Pure function - no I/O.
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
                    transcript.strand(), transcript.exonSpans()));
        }

        return new PackResult(altContig, sequence.toString(), entries);
    }

    public record PackResult(String altContig, String sequence, List<ContigEntry> entries)
    {
        public boolean isEmpty() { return entries.isEmpty(); }
    }
}
