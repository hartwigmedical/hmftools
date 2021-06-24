package com.hartwig.hmftools.sage.phase;

import java.util.Collection;
import java.util.List;
import java.util.Optional;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.genome.region.HmfExonRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.select.TranscriptRegionSelector;
import com.hartwig.hmftools.sage.variant.SageVariant;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class PhasedInframeIndel extends BufferedPostProcessor
{
    private static final int MAX_DISTANCE = 50;

    private int mPhase;
    private final TranscriptRegionSelector mSelector;

    public PhasedInframeIndel(final Consumer<SageVariant> consumer, final List<HmfTranscriptRegion> transcripts)
    {
        super(MAX_DISTANCE, consumer);
        mSelector = new TranscriptRegionSelector(transcripts);
    }

    @Override
    protected void processSageVariant(@NotNull final SageVariant newVariant, @NotNull final Collection<SageVariant> buffer)
    {
        final int lps = newVariant.localPhaseSet();
        if(lps == 0 || !newVariant.isPassing() || !newVariant.variant().isFrameshiftIndel())
        {
            return;
        }

        final Optional<HmfExonRegion> maybeExon = frameshiftExon(newVariant.variant());
        if(!maybeExon.isPresent())
        {
            return;
        }

        final int newLength = newVariant.variant().indelLength();

        for(final SageVariant other : buffer)
        {
            if(other.localPhaseSet() == lps && other.isPassing() && other.variant().isFrameshiftIndel())
            {
                final Optional<HmfExonRegion> maybeOtherExon = frameshiftExon(other.variant());
                if(!maybeOtherExon.filter(x -> x.equals(maybeExon.get())).isPresent())
                {
                    continue;
                }

                int otherLength = other.variant().indelLength();
                int combinedLength = newLength + otherLength;
                if(combinedLength % 3 == 0)
                {
                    if(other.phasedInframeIndel() != 0)
                    {
                        newVariant.phasedInframeIndel(other.phasedInframeIndel());
                    }
                    else
                    {
                        mPhase++;
                        newVariant.phasedInframeIndel(mPhase);
                        other.phasedInframeIndel(mPhase);
                    }
                }
            }
        }
    }

    @NotNull
    private Optional<HmfExonRegion> frameshiftExon(@Nullable final VariantHotspot hotspot)
    {
        return Optional.ofNullable(hotspot).filter(VariantHotspot::isFrameshiftIndel).flatMap(x -> selectExon(mSelector, x.position() + 1));
    }

    @NotNull
    static Optional<HmfExonRegion> selectExon(@NotNull final TranscriptRegionSelector selector, long position)
    {
        return selector.select(position).map(x -> selectExon(position, x));
    }

    @Nullable
    static HmfExonRegion selectExon(final long position, @NotNull final HmfTranscriptRegion transcript)
    {
        for(HmfExonRegion exon : transcript.exome())
        {
            if(position >= exon.start() && position <= exon.end())
            {
                return exon;
            }
        }
        return null;
    }
}
