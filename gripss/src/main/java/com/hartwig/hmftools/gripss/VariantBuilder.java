package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.createSingleBreakend;
import static com.hartwig.hmftools.gripss.filters.CommonFilters.isPolyATSequence;
import static com.hartwig.hmftools.gripss.common.SvData.hasLength;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantFactory;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;
import com.hartwig.hmftools.gripss.common.GenotypeIds;
import com.hartwig.hmftools.gripss.common.SvData;
import com.hartwig.hmftools.gripss.filters.FilterConstants;
import com.hartwig.hmftools.gripss.filters.HardFilters;
import com.hartwig.hmftools.gripss.filters.HotspotCache;

import htsjdk.variant.variantcontext.VariantContext;

public class VariantBuilder
{
    private final HardFilters mHardFilters;
    private final HotspotCache mHotspotCache;
    private StructuralVariantFactory mSvFactory;
    private final Set<String> mHotspotCandidateVcfIds;
    private final Set<String> mHardFilteredVcfIds;

    private int mHardFilteredCount;

    public VariantBuilder(final FilterConstants filterConstants, final HotspotCache hotspotCache)
    {
        mHardFilters = new HardFilters(filterConstants);
        mHotspotCache = hotspotCache;

        mSvFactory = new StructuralVariantFactory(new AlwaysPassFilter());
        mHardFilteredVcfIds = Sets.newHashSet();
        mHotspotCandidateVcfIds = Sets.newHashSet();
        mHardFilteredCount = 0;
    }

    public int hardFilteredCount() { return mHardFilteredCount; }
    public int incompleteSVs() { return mSvFactory.unmatched().size(); }

    public void clearState()
    {
        mSvFactory.clear();
        mHardFilteredCount = 0;
        mHardFilteredVcfIds.clear();
        mHotspotCandidateVcfIds.clear();
    }

    public SvData checkCreateVariant(final VariantContext variant, final GenotypeIds genotypeIds)
    {
        // each SV breakend can be a) not hard-filtered, b) hard-filtered but a hotspot candidate or c) neither
        // and if it's not a single then these 3 scenarios need to be considered together for the pair of breakends
        // if either are hard-filtered and not hotspot candidates, then drop them both
        boolean hardFiltered = mHardFilters.isFiltered(variant, genotypeIds);

        if(StructuralVariantFactory.isSingleBreakend(variant))
        {
            if(hardFiltered)
            {
                ++mHardFilteredCount;
                return null;
            }

            StructuralVariant sv = createSingleBreakend(variant);
            return new SvData(sv, genotypeIds);
        }

        String mateId = StructuralVariantFactory.mateId(variant);

        if(mateId == null)
            return null;

        if(mHardFilteredVcfIds.contains(mateId))
        {
            // drop since the mate was hard-filtered and not a hotspot candidate, so doesn't matter what this breakend is
            mHardFilteredVcfIds.remove(mateId);
            ++mHardFilteredCount;
            return null;
        }

        boolean mateHotspotCandidate = mHotspotCandidateVcfIds.contains(mateId);
        boolean hotspotCandidate = hardFiltered && mHotspotCache.matchesHotspotBreakend(variant.getContig(), variant.getStart());

        if(hardFiltered && !hotspotCandidate)
        {
            // clean up the first leg if it wasn't hard-filtered
            if(mateHotspotCandidate)
            {
                mHotspotCandidateVcfIds.remove(mateId);
                ++mHardFilteredCount;
            }
            else if(mSvFactory.hasUnmatchedVariant(mateId))
            {
                mSvFactory.removeUnmatchedVariant(mateId);
                ++mHardFilteredCount;
            }
            else
            {
                mHardFilteredVcfIds.add(variant.getID());
            }

            return null;
        }

        // remaining scenarios: neither leg hard-filtered, or one or both are hotspot candidates

        int currentSvCount = mSvFactory.results().size();
        mSvFactory.addVariantContext(variant);

        // check if both breakends have now been encountered
        if(currentSvCount == mSvFactory.results().size())
        {
            // this is the first breakend in the SV, so cache its state accordingly
            if(hotspotCandidate)
            {
                mHotspotCandidateVcfIds.add(variant.getID());
            }
            else if(hardFiltered)
            {
                mHardFilteredVcfIds.add(variant.getID());
            }

            return null;
        }

        final StructuralVariant sv = popLastSv(); // get and clear from storage

        if(sv == null)
            return null;

        if(mateHotspotCandidate)
            mHotspotCandidateVcfIds.remove(mateId);

        if(hotspotCandidate || mateHotspotCandidate)
        {
            // check whether this SV can be rescued as a hotspot
            if(!keepHotspotVariant(sv))
            {
                ++mHardFilteredCount;
                return null;
            }
        }

        return new SvData(sv, genotypeIds);
    }

    private boolean keepHotspotVariant(final StructuralVariant sv)
    {
        // check hotspot rescue
        if(hasLength(sv.type()) && SvData.length(sv) < FilterConstants.SHORT_RESCUE_LENGTH)
            return false;
        else if(isPolyATSequence(sv.startContext()) || (sv.endContext() != null && isPolyATSequence(sv.endContext())))
            return false;

        return mHotspotCache.matchesHotspot(sv);
    }

    private final StructuralVariant popLastSv()
    {
        if(mSvFactory.results().isEmpty())
            return null;

        StructuralVariant sv = mSvFactory.results().get(0);
        mSvFactory.results().remove(0);

        return sv;
    }

}
