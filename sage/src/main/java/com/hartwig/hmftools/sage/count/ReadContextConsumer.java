package com.hartwig.hmftools.sage.count;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.position.GenomePosition;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class ReadContextConsumer implements GenomePosition, Consumer<SAMRecord> {
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

    @NotNull
    @Override
    public String chromosome() {
        return hotspot.chromosome();
    }

    @Override
    public long position() {
        return hotspot.position();
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

            byte[] readBases = record.getReadBases();
            for (int readBasePosition = 0; readBasePosition < readBases.length; readBasePosition++) {
                long refPosition = record.getReferencePositionAtReadPosition(readBasePosition + 1); //TODO: Check + 1

                final ReadContext newContext = new ReadContext(readBasePosition, readBases);
                ReadContext.ReadContextMatch match = readContext.match(newContext);
                if (!match.equals(ReadContext.ReadContextMatch.NONE)) {
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
