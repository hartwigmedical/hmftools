package com.hartwig.hmftools.sage.phase;

import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Optional;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.HmfExonRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotComparator;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.select.RegionSelector;
import com.hartwig.hmftools.sage.variant.SageVariant;

import org.jetbrains.annotations.NotNull;

public class RightAlignMicrohomology implements Consumer<SageVariant>
{
    private final Consumer<SageVariant> mConsumer;
    private final RegionSelector<HmfTranscriptRegion> mSelector;
    private final List<SageVariant> mRightAlignedList = Lists.newArrayList();
    private final Comparator<VariantHotspot> mComparator = new VariantHotspotComparator();

    public RightAlignMicrohomology(final Consumer<SageVariant> consumer, final List<HmfTranscriptRegion> transcripts)
    {
        mSelector = new RegionSelector<>(transcripts);
        mConsumer = consumer;
    }

    @Override
    public void accept(final SageVariant variant)
    {
        if(!realign(variant))
        {
            Iterator<SageVariant> realignedIterator = mRightAlignedList.iterator();
            while(realignedIterator.hasNext())
            {
                final SageVariant rightAligned = realignedIterator.next();
                if(mComparator.compare(rightAligned.variant(), variant.variant()) < 0)
                {
                    mConsumer.accept(rightAligned);
                    realignedIterator.remove();
                }
            }

            mConsumer.accept(variant);
        }
    }

    boolean realign(@NotNull final SageVariant variant)
    {
        if(!variant.isPassing())
        {
            return false;
        }

        final boolean isInframeDel = variant.variant().isInframeIndel() && variant.alt().length() == 1;
        if(!variant.variant().isInframeIndel())
        {
            return false;
        }

        final String microhomology = variant.microhomology();
        final int microhomologyLength = microhomology.length();
        if(microhomologyLength == 0)
        {
            return false;
        }

        final Optional<HmfExonRegion> maybeLeftAlignedRegion = PhasedInframeIndel.selectExon(mSelector, variant.position() + 1);
        if(maybeLeftAlignedRegion.isPresent())
        {
            return false;
        }

        final VariantHotspot rightAligned =
                isInframeDel ? rightAlignDel(variant.variant(), microhomology) : rightAlignIns(variant.variant(), microhomology);
        final Optional<HmfTranscriptRegion> maybeRightAlignedTranscript = mSelector.select(rightAligned.position() + 1);
        if(!maybeRightAlignedTranscript.isPresent())
        {
            return false;
        }

        final HmfTranscriptRegion rightAlignedTranscript = maybeRightAlignedTranscript.get();
        if(rightAligned.position() + 1 < rightAlignedTranscript.codingStart()
                || rightAligned.position() + 1 > rightAlignedTranscript.codingEnd())
        {
            return false;
        }

        final Optional<HmfExonRegion> maybeRightAlignedExon =
                Optional.ofNullable(PhasedInframeIndel.selectExon(rightAligned.position() + 1, rightAlignedTranscript));
        if(!maybeRightAlignedExon.isPresent())
        {
            return false;
        }

        final Candidate oldCandidate = variant.candidate();
        final Candidate newCandidate =
                new Candidate(oldCandidate.tier(), rightAligned, oldCandidate.readContext(), oldCandidate.maxReadDepth(), oldCandidate.minNumberOfEvents());

        final SageVariant realignedVariant =
                new SageVariant(newCandidate, variant.filters(), variant.normalAltContexts(), variant.tumorAltContexts());
        realignedVariant.realigned(true);
        mRightAlignedList.add(realignedVariant);

        return true;
    }

    @NotNull
    static VariantHotspot rightAlignDel(@NotNull final VariantHotspot variant, @NotNull final String microhomology)
    {
        int microhomologyLength = microhomology.length();
        final String alt = microhomology.substring(microhomologyLength - 1);
        final String ref = variant.ref().substring(microhomologyLength) + microhomology;
        return ImmutableVariantHotspotImpl.builder()
                .from(variant)
                .position(variant.position() + microhomologyLength)
                .ref(ref)
                .alt(alt)
                .build();
    }

    @NotNull
    static VariantHotspot rightAlignIns(@NotNull final VariantHotspot variant, @NotNull final String microhomology)
    {
        int microhomologyLength = microhomology.length();
        final String alt = variant.alt().substring(microhomologyLength) + microhomology;
        final String ref = microhomology.substring(microhomologyLength - 1);
        return ImmutableVariantHotspotImpl.builder()
                .from(variant)
                .position(variant.position() + microhomologyLength)
                .ref(ref)
                .alt(alt)
                .build();
    }

    public void flush()
    {
        mRightAlignedList.sort((o1, o2) -> mComparator.compare(o1.variant(), o2.variant()));
        mRightAlignedList.forEach(mConsumer);
        mRightAlignedList.clear();
    }
}
