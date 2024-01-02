package com.hartwig.hmftools.sage.candidate_;

import static com.hartwig.hmftools.sage.SageConstants.MIN_SECOND_CANDIDATE_FULL_READS;
import static com.hartwig.hmftools.sage.SageConstants.MIN_SECOND_CANDIDATE_FULL_READS_PERC;

import java.util.Collections;
import java.util.List;

import javax.annotation.Nullable;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.primitives.Ints;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.common.ReadContextMatch;

/**
 *  AltContext keeps track of interim ReadContextCandidates, and also the final selected candidate.
 *  Also keeps track of rawSupport (number of reads), and total rase base quality.
 * <p>
 * mRawAltSupport - number of alt reads.
 * mRawBaseQualityAlt - base quality sum of all alt reads.
 * mSecondCandidate - add a second candidate if its core is different and it has sufficient support
 */
public class AltContext_ implements VariantHotspot
{
    public final String Ref;
    public final String Alt;
    public final RefContext_ RefContext;

    private final List<ReadContextCandidate_> mInterimReadContextCandidates;

    private int mRawSupportAlt;
    private int mRawBaseQualityAlt;
    private ReadContextCandidate_ mCandidate;
    private AltContext_ mSecondCandidate; // relevant if has a different read context and sufficient support

    public AltContext_(final RefContext_ refContext, final String ref, final String alt)
    {
        RefContext = refContext;
        Ref = ref;
        Alt = alt;

        mInterimReadContextCandidates = Lists.newArrayList();

        mRawBaseQualityAlt = 0;
        mRawSupportAlt = 0;
        mCandidate = null;
        mSecondCandidate = null;
    }

    public AltContext_(
            final RefContext_ refContext, final String ref, final String alt, final ReadContextCandidate_ candidate,
            int rawSupportAlt, int rawBaseQualAlt)
    {
        RefContext = refContext;
        Ref = ref;
        Alt = alt;

        mCandidate = candidate;
        mRawSupportAlt = rawSupportAlt;
        mRawBaseQualityAlt = rawBaseQualAlt;

        mInterimReadContextCandidates = null;
    }

    /**
     * Register an alt read.
     */
    public void incrementAltRead(int baseQuality)
    {
        mRawSupportAlt++;
        mRawBaseQualityAlt += baseQuality;
    }

    /**
     * Adds a new read context.
     * <p>
     * Checks current interim candidates and increments counters based on a FULL/PARTIAL/CORE match for the newReadContext against these.
     * <p>
     * If no full match is found, the add a new interim candidate with either a core or partial match count bases on what was found above.
     * <p>
     * If full match is found and the newReadContext has a bigger buffer than the found fullMatchCandidate then swap these out, and absorb
     * the fullMatchCandidate counters into the new candidate.
     */
    public void addReadContext(int numberOfEvents, final ReadContext_ newReadContext)
    {
        int partialMatchCount = 0;
        int coreMatchCount = 0;
        ReadContextCandidate_ fullMatchCandidate = null;

        for (ReadContextCandidate_ candidate : mInterimReadContextCandidates)
        {
            ReadContextMatch matchType = candidate.ReadContext.indexedBases().matchAtPosition(newReadContext.indexedBases());
            switch (matchType)
            {
                case FULL:
                    candidate.incrementFull(1, numberOfEvents);
                    fullMatchCandidate = candidate;
                    break;
                case PARTIAL:
                    ++candidate.PartialMatch;
                    ++partialMatchCount;
                    break;
                case CORE:
                    ++candidate.CoreMatch;
                    ++coreMatchCount;
                    break;
            }
        }

        if (fullMatchCandidate == null)
        {
            ReadContextCandidate_ newCandidate = new ReadContextCandidate_(numberOfEvents, newReadContext);
            newCandidate.CoreMatch += coreMatchCount;
            newCandidate.PartialMatch += partialMatchCount;
            mInterimReadContextCandidates.add(newCandidate);
        }
        else if (newReadContext.maxFlankLength() > fullMatchCandidate.ReadContext.maxFlankLength())
        {
            // Full match and newReadCount has a bigger buffer around the core, so use this instead.
            mInterimReadContextCandidates.remove(fullMatchCandidate);
            ReadContextCandidate_ newCandidate = new ReadContextCandidate_(numberOfEvents, newReadContext);
            newCandidate.CoreMatch += fullMatchCandidate.CoreMatch;
            newCandidate.PartialMatch += fullMatchCandidate.PartialMatch;
            newCandidate.incrementFull(fullMatchCandidate.FullMatch, fullMatchCandidate.MinNumberOfEvents);
            mInterimReadContextCandidates.add(newCandidate);
        }
    }

    /**
     * Filter out interim read contexts with incomplete flanks, sorts interim candidates by DESC full, DESC partial, and then DESC core read
     * counts, sets the top one as the candidate. Add a second if its core is different (i.e. does not contain and is not contained within
     * the top interim candidate's core string) and it has sufficient support. Finally, clears the interim read contexts.
     */
    public void selectCandidates()
    {
        if (mInterimReadContextCandidates == null)
            return;

        // Remove all that has incomplete flanks.
        mInterimReadContextCandidates.removeIf(x -> x.ReadContext.hasIncompleteFlanks());

        // If there are no interim candidates then just return.
        if (mInterimReadContextCandidates.isEmpty())
            return;

        // Sort by DESC full, DESC partial, and then DESC core read counts.
        Collections.sort(mInterimReadContextCandidates);

        // Set candidate to the first one.
        mCandidate = mInterimReadContextCandidates.get(0);
        String topCore = mCandidate.ReadContext.coreString();
        double fullThreshold = MIN_SECOND_CANDIDATE_FULL_READS_PERC * mCandidate.FullMatch;

        // add a second if its core is different and it has sufficient support
        for (int i = 1; i < mInterimReadContextCandidates.size(); ++i)
        {
            ReadContextCandidate_ readContextCandidate = mInterimReadContextCandidates.get(i);
            if (readContextCandidate.FullMatch < fullThreshold || readContextCandidate.FullMatch < MIN_SECOND_CANDIDATE_FULL_READS)
                break;

            String candidateCore = readContextCandidate.ReadContext.coreString();
            if (candidateCore.contains(topCore) || topCore.contains(candidateCore))
                continue;

            mSecondCandidate = new AltContext_(RefContext, Ref, Alt, readContextCandidate, mRawSupportAlt, mRawBaseQualityAlt);
            break;
        }

        // Clear interim candidates
        mInterimReadContextCandidates.clear();
    }

    /**
     * Full + partial matches of the mCandidate.
     */
    public int readContextSupport()
    {
        return mCandidate.count();
    }

    /**
     * MinNumberOfEvents of the mCandidate.
     */
    public int minNumberOfEvents()
    {
        return mCandidate.MinNumberOfEvents;
    }

    @VisibleForTesting
    public List<ReadContextCandidate_> interimReadContexts()
    {
        return mInterimReadContextCandidates;
    }

    public boolean hasValidCandidate() { return mCandidate != null; }
    public boolean hasSecondCandidate() { return mSecondCandidate != null; }
    public AltContext_ secondCandidate() { return mSecondCandidate; }

    public ReadContext_ readContext() { return mCandidate.ReadContext; }

    @Override
    public String ref() { return Ref; }

    @Override
    public String alt() { return Alt; }

    @Override
    public String chromosome() { return RefContext.chromosome(); }

    @Override
    public int position() { return RefContext.position(); }

    public int rawAltSupport() { return mRawSupportAlt; }

    public int rawAltBaseQuality() { return mRawBaseQualityAlt; }

    private boolean equalTo(final VariantHotspot another)
    {
        return ref().equals(another.ref()) && alt().equals(another.alt()) && chromosome().equals(another.chromosome())
                && position() == another.position();
    }

    /**
     * Just based of ref, alt, chromosome, and position, and so just the variant.
     */
    @Override
    public boolean equals(@Nullable Object another)
    {
        if(this == another)
            return true;

        return another instanceof VariantHotspot && equalTo((VariantHotspot) another);
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

    @Override
    public String toString()
    {
        return String.format("var(%s:%d %s->%s) readCandidates(%d)",
                chromosome(), position(), Ref, Alt, mInterimReadContextCandidates != null ? mInterimReadContextCandidates.size() : 0);
    }
}
