package com.hartwig.hmftools.sage.candidate;

import static java.lang.Math.min;

import static com.hartwig.hmftools.sage.SageConfig.isUltima;
import static com.hartwig.hmftools.sage.SageConstants.MIN_SECOND_CANDIDATE_FULL_READS;
import static com.hartwig.hmftools.sage.SageConstants.MIN_SECOND_CANDIDATE_FULL_READS_PERC;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.CORE;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.FULL;

import java.util.Collections;
import java.util.List;
import java.util.Set;

import javax.annotation.Nullable;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.google.common.primitives.Ints;
import com.hartwig.hmftools.common.sequencing.UltimaBamUtils;
import com.hartwig.hmftools.sage.common.ReadMatchInfo;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.common.VariantReadContextBuilder;
import com.hartwig.hmftools.sage.filter.FilterConfig;

import htsjdk.samtools.SAMRecord;

public class AltContext extends SimpleVariant
{
    public final RefContext RefContext;
    
    private final List<ReadContextCandidate> mReadContextCandidates;

    private boolean mAboveMinAltSupport;
    private final Set<String> mUniqueReadIds;
    private ReadContextCandidate mCandidate;
    private AltContext mSecondCandidate; // relevant if has a different read context and sufficient support

    public AltContext(final RefContext refContext, final String ref, final String alt)
    {
        super(refContext.Chromosome, refContext.Position, ref, alt);
        RefContext = refContext;

        mReadContextCandidates = Lists.newArrayList();
        mUniqueReadIds = Sets.newHashSet();
        mAboveMinAltSupport = false;
        mCandidate = null;
        mSecondCandidate = null;
    }

    private AltContext(final RefContext refContext, final String ref, final String alt, final ReadContextCandidate candidate)
    {
        super(refContext.Chromosome, refContext.Position, ref, alt);
        RefContext = refContext;

        mCandidate = candidate;
        mAboveMinAltSupport = true;
        mUniqueReadIds = null;
        mReadContextCandidates = null;
    }

    public VariantReadContext readContext() { return mCandidate.readContext(); }

    @Override
    public String ref() { return Ref; }

    @Override
    public String alt() { return Alt; }

    @Override
    public String chromosome() { return RefContext.chromosome(); }

    @Override
    public int position() { return RefContext.position(); }

    public int fullMatchSupport() { return mCandidate.FullMatch; }
    public int coreMatchSupport() { return mCandidate.CoreMatch; }
    public int lowQualInCoreCount() { return mCandidate.LowQualInCoreCount; }
    public int minNumberOfEvents() { return mCandidate.MinNumberOfEvents; }

    public boolean hasValidCandidate() { return mCandidate != null; }
    public boolean hasSecondCandidate() { return mSecondCandidate != null; }
    public AltContext secondCandidate() { return mSecondCandidate; }

    public void addReadContext(
            int numberOfEvents, final SAMRecord read, final int variantReadIndex,
            final VariantReadContextBuilder readContextBuilder, final RefSequence refSequence)
    {
        int coreMatch = 0;
        ReadContextCandidate fullMatchCandidate = null;

        for(ReadContextCandidate candidate : mReadContextCandidates)
        {
            // compare the core and flanks for the 2 contexts, not allowing for mismatches
            ReadMatchInfo matchInfo = candidate.matcher().determineReadMatchInfo(
                    read.getReadBases(), null, variantReadIndex, true);

            switch(matchInfo.MatchType)
            {
                case FULL:
                    candidate.incrementFull(1, numberOfEvents);
                    fullMatchCandidate = candidate;

                    break;

                case CORE:
                    candidate.CoreMatch++;
                    coreMatch++;
                    break;
            }

            if(isUltima())
            {
                if(matchInfo.MatchType == FULL || matchInfo.MatchType == CORE)
                {
                    if(lowQualInCore(candidate.readContext(), read, variantReadIndex))
                        ++candidate.LowQualInCoreCount;
                }
            }
        }

        if(fullMatchCandidate == null)
        {
            VariantReadContext readContext = readContextBuilder.createContext(this, read, variantReadIndex, refSequence);

            if(readContext == null)
                return;

            ReadContextCandidate candidate = new ReadContextCandidate(numberOfEvents, readContext);
            candidate.CoreMatch += coreMatch;

            if(isUltima() && lowQualInCore(candidate.readContext(), read, variantReadIndex))
                ++candidate.LowQualInCoreCount;

            mReadContextCandidates.add(candidate);
        }

        // keep enough reads to test the (unique) raw alt support limit
        if(!mAboveMinAltSupport && mUniqueReadIds.size() < FilterConfig.HardMinTumorRawAltSupport)
        {
            mUniqueReadIds.add(read.getReadName());

            if(mUniqueReadIds.size() >= FilterConfig.HardMinTumorRawAltSupport)
            {
                mAboveMinAltSupport = true;
                mUniqueReadIds.clear();
            }
        }
    }

    public void selectCandidates()
    {
        if(mReadContextCandidates.isEmpty())
            return;

        // sort by full, then partial then core read counts
        Collections.sort(mReadContextCandidates);

        if(mReadContextCandidates.isEmpty())
            return;

        mCandidate = mReadContextCandidates.get(0);

        if(mReadContextCandidates.size() > 1)
        {
            final String topCore = mReadContextCandidates.get(0).readContext().coreStr();
            double topCandidateRcThreshold = mReadContextCandidates.get(0).FullMatch * MIN_SECOND_CANDIDATE_FULL_READS_PERC;

            // add a second if its core is different and it has sufficient support
            for(int i = 0; i < mReadContextCandidates.size(); ++i)
            {
                ReadContextCandidate candidate = mReadContextCandidates.get(i);

                if(candidate.FullMatch < MIN_SECOND_CANDIDATE_FULL_READS || candidate.FullMatch < topCandidateRcThreshold)
                    break;

                String coreStr = candidate.readContext().coreStr();
                if(coreStr.contains(topCore) || topCore.contains(coreStr))
                    continue;

                mSecondCandidate = new AltContext(RefContext, Ref, Alt, candidate);
                break;
            }
        }

        mReadContextCandidates.clear();
    }

    private boolean lowQualInCore(final VariantReadContext readContext, final SAMRecord read, int variantReadIndex)
    {
        List<Integer> lowQualIndices = UltimaBamUtils.extractLowQualTag(read);

        if(lowQualIndices.isEmpty())
            return false;

        int coreIndexStart = variantReadIndex - readContext.leftCoreLength();
        int coreIndexEnd = variantReadIndex + readContext.rightCoreLength();

        for(int index = coreIndexStart; index <= coreIndexEnd; ++index)
        {
            if(lowQualIndices.contains(index))
                return true;
        }

        return false;
    }

    public boolean aboveMinAltSupport() { return mAboveMinAltSupport; }

    @Override
    public boolean equals(@Nullable Object another)
    {
        if(this == another)
            return true;

        return another instanceof SimpleVariant && equalTo((SimpleVariant) another);
    }

    private boolean equalTo(final SimpleVariant other)
    {
        return Ref.equals(other.Ref) && Alt.equals(other.alt()) && Chromosome.equals(other.Chromosome) && Position == other.Position;
    }

    @Override
    public int hashCode()
    {
        int h = 5381;
        h += (h << 5) + Ref.hashCode();
        h += (h << 5) + Alt.hashCode();
        h += (h << 5) + chromosome().hashCode();
        h += (h << 5) + Ints.hashCode(position());
        return h;
    }

    public String toString()
    {
        return String.format("var(%s:%d %s->%s) readCandidates(%d)",
                chromosome(), position(), Ref, Alt, mReadContextCandidates != null ? mReadContextCandidates.size() : 0);
    }

}
