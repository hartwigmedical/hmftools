package com.hartwig.hmftools.sage.candidate;

import static com.hartwig.hmftools.sage.SageConfig.isUltima;
import static com.hartwig.hmftools.sage.SageConstants.MIN_SECOND_CANDIDATE_FULL_READS;
import static com.hartwig.hmftools.sage.SageConstants.MIN_SECOND_CANDIDATE_FULL_READS_PERC;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.CORE;
import static com.hartwig.hmftools.sage.common.ReadContextMatch.FULL;
import static com.hartwig.hmftools.sage.filter.FilterConfig.ULTIMA_CANDIDATE_MIN_HIGH_BQ_THRESHOLD;
import static com.hartwig.hmftools.sage.seqtech.UltimaUtils.lowQualInReadContextCore;
import static com.hartwig.hmftools.sage.seqtech.UltimaUtils.lowQualInReadCore;

import java.util.Collections;
import java.util.List;
import java.util.Set;

import javax.annotation.Nullable;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.google.common.primitives.Ints;
import com.hartwig.hmftools.common.sequencing.UltimaBamUtils;
import com.hartwig.hmftools.sage.SageCallConfig;
import com.hartwig.hmftools.sage.common.ReadMatchInfo;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.common.VariantReadContextBuilder;
import com.hartwig.hmftools.sage.filter.FilterConfig;
import com.hartwig.hmftools.sage.seqtech.UltimaUtils;

import htsjdk.samtools.SAMRecord;

public class AltContext extends SimpleVariant
{
    public final RefContext RefContext;
    
    private final List<ReadContextCandidate> mReadContextCandidates;

    private boolean mAboveMinAltSupport;
    private final Set<String> mUniqueReadIds;

    private Boolean mCheckHighQualCoreSupport;
    private int mHighQualCoreSupport;
    private boolean mHasHighQualCoreSupport;

    private ReadContextCandidate mCandidate;
    private ReadContextCandidate mSecondCandidate; // relevant if has a different read context and sufficient support

    public AltContext(final RefContext refContext, final String ref, final String alt)
    {
        super(refContext.Chromosome, refContext.Position, ref, alt);
        RefContext = refContext;

        mReadContextCandidates = Lists.newArrayList();
        mUniqueReadIds = Sets.newHashSet();
        mAboveMinAltSupport = false;
        mHighQualCoreSupport = 0;
        mCandidate = null;
        mSecondCandidate = null;

        mCheckHighQualCoreSupport = isUltima() && ULTIMA_CANDIDATE_MIN_HIGH_BQ_THRESHOLD > 0 ? null : false;
        mHasHighQualCoreSupport = mCheckHighQualCoreSupport != null;
    }

    public VariantReadContext readContext() { return mCandidate.readContext(); }

    public String chromosome() { return RefContext.chromosome(); }
    public int position() { return RefContext.position(); }
    public String ref() { return Ref; }
    public String alt() { return Alt; }

    public boolean hasMinAltSupport() { return mAboveMinAltSupport && mHasHighQualCoreSupport; }

    public boolean hasValidCandidate() { return mCandidate != null; }
    public ReadContextCandidate candidate() { return mCandidate; }

    public boolean hasSecondCandidate() { return mSecondCandidate != null; }
    public ReadContextCandidate secondCandidate() { return mSecondCandidate; }

    public void addReadContext(
            int numberOfEvents, final SAMRecord read, final int variantReadIndex,
            final VariantReadContextBuilder readContextBuilder, final RefSequence refSequence)
    {
        int fullMatchCount = 0;
        ReadContextCandidate fullMatchCandidate = null;

        if(mCheckHighQualCoreSupport == null)
        {
            // check if this variant ought to check low-qual bases for this and subsequent reads
            mCheckHighQualCoreSupport = !UltimaUtils.hasLongHomopolymerInContext(this, read, variantReadIndex);
            mHasHighQualCoreSupport = !mCheckHighQualCoreSupport;
        }

        List<Integer> lowQualIndices = null;

        if(mCheckHighQualCoreSupport && (!mHasHighQualCoreSupport || SageCallConfig.LogCandidates))
        {
            // only extract Ultima low-qual bases if they are required
            lowQualIndices = UltimaBamUtils.extractLowQualIndices(read);
        }

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
                    fullMatchCount += candidate.FullMatch;
                    break;
            }

            if(lowQualIndices != null && (matchInfo.MatchType == FULL || matchInfo.MatchType == CORE))
            {
                // skip checking if not required for candidate logging
                if(UltimaUtils.lowQualInReadContextCore(variantReadIndex, lowQualIndices, candidate.readContext()))
                    ++candidate.LowQualInCoreCount;
                else
                    ++mHighQualCoreSupport;
            }
        }

        if(fullMatchCandidate == null)
        {
            VariantReadContext readContext = null;

            // skip creating candidates if this variant checks low-quals in core and finds them present
            if(!mCheckHighQualCoreSupport || !UltimaUtils.lowQualInReadCore(this, variantReadIndex, lowQualIndices))
            {
                readContext = readContextBuilder.createContext(this, read, variantReadIndex, refSequence);

                if(readContext == null)
                    return;
            }

            boolean hasValidCore = readContext != null;

            if(hasValidCore && lowQualIndices != null)
            {
                // now the core has been established, check again
                if(lowQualInReadContextCore(variantReadIndex, lowQualIndices, readContext))
                    hasValidCore = false;
                else
                    ++mHighQualCoreSupport;
            }

            if(hasValidCore)
            {
                ReadContextCandidate candidate = new ReadContextCandidate(numberOfEvents, readContext);
                candidate.CoreMatch = fullMatchCount;

                mReadContextCandidates.add(candidate);
            }
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

        if(!mHasHighQualCoreSupport)
            mHasHighQualCoreSupport = mHighQualCoreSupport >= ULTIMA_CANDIDATE_MIN_HIGH_BQ_THRESHOLD;
    }

    public void selectCandidates()
    {
        if(isUltima() && !mHasHighQualCoreSupport)
        {
            mReadContextCandidates.clear();
            return;
        }

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

                mSecondCandidate = candidate;
                break;
            }
        }

        mReadContextCandidates.clear();
    }

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
