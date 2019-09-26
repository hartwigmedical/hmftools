package com.hartwig.hmftools.sage.count;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.hotspot.VariantHotspot;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;

public class ReadContextConsumer implements Consumer<SAMRecord> {
    private final VariantHotspot hotspot;
    private final ReadContext readContext;

    private int full;
    private int partial;
    private int missaligned;

    public ReadContextConsumer(@NotNull final VariantHotspot hotspot, @NotNull final ReadContext readContext) {
        assert (readContext.isComplete());
        this.hotspot = hotspot;
        this.readContext = readContext;
    }

    public VariantHotspot hotspot() {
        return hotspot;
    }

    @NotNull
    public ReadContext readContext() {
        return readContext;
    }

    public int full() {
        return full;
    }

    public int partial() {
        return partial;
    }

    public int missaligned() {
        return missaligned;
    }

    @Override
    public void accept(final SAMRecord record) {

        if (record.getAlignmentStart() <= hotspot.position() && record.getAlignmentEnd() >= hotspot.position()) {
            for (AlignmentBlock alignmentBlock : record.getAlignmentBlocks()) {
                final int readBaseStart = alignmentBlock.getReadStart() - 1;

                for (int i = 0; i < alignmentBlock.getLength(); i++) {
                    long refPosition = alignmentBlock.getReferenceStart() + i;

                    final int readBasePosition = readBaseStart + i;
                    final ReadContext newContext = new ReadContext(readBasePosition, record.getReadBases());
                    ReadContext.ReadContextMatch match = readContext.match(newContext);
                    if (match.equals(ReadContext.ReadContextMatch.NONE)) {
                        if (refPosition == hotspot.position()) {
                            if (match.equals(ReadContext.ReadContextMatch.FULL)) {
                                full++;
                            } else {
                                partial++;
                            }
                        } else if (match.equals(ReadContext.ReadContextMatch.FULL)) {
                            missaligned++;
                        }
                    }
                }
            }
        }
    }
}
